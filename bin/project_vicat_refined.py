#!/usr/bin/env python3
"""Project precomputed viCAT loci onto provirus-refined sequence intervals."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from evidence_schema import TAXONOMY_COLUMNS, unclassified_taxonomy
from standardize_vicat import (
    OUTPUT_COLUMNS, aggregate_taxonomy, cluster_flanks, format_number, viral_clusters,
)


PROJECTION_COLUMNS = [
    "sample_id", "sequence_id", "parent_sequence_id", "record_type",
    "source_locus_id", "source_coordinates", "projected_coordinates",
    "locus_classification", "projection_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--loci", required=True, type=Path)
    parser.add_argument("--region-map", required=True, type=Path)
    parser.add_argument("--contig-taxonomy-support", required=True, type=float)
    parser.add_argument("--cluster-min-viral-loci", required=True, type=int)
    parser.add_argument("--cluster-max-neutral-gap", required=True, type=int)
    parser.add_argument("--output-evidence", required=True, type=Path)
    parser.add_argument("--output-projection", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return list(reader)


def coordinates(value: str) -> tuple[int, int]:
    start, end = value.split("-", 1)
    parsed = int(start), int(end)
    if parsed[0] < 1 or parsed[1] < parsed[0]:
        raise ValueError(f"Invalid coordinates: {value}")
    return parsed


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> None:
    locus_rows = read_tsv(
        args.loci,
        {"sample_id", "sequence_id", "locus_id", "coordinates", "locus_classification"}
        | set(TAXONOMY_COLUMNS),
    )
    region_rows = read_tsv(
        args.region_map,
        {
            "sample_id", "sequence_id", "parent_sequence_id", "record_type",
            "coordinates", "refined_length",
        },
    )
    by_parent: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in locus_rows:
        if row["sample_id"] != args.sample_id:
            raise ValueError(f"viCAT locus sample mismatch: {row['sample_id']}")
        row["_start"], row["_end"] = coordinates(row["coordinates"])
        by_parent[row["sequence_id"]].append(row)

    evidence_rows: list[dict[str, object]] = []
    projection_rows: list[dict[str, object]] = []
    for region in region_rows:
        if region["sample_id"] != args.sample_id:
            raise ValueError(f"Region-map sample mismatch: {region['sample_id']}")
        final_id = region["sequence_id"]
        parent_id = region["parent_sequence_id"].strip() or final_id
        if region["coordinates"].strip():
            region_start, region_end = coordinates(region["coordinates"])
        else:
            region_start, region_end = 1, int(region["refined_length"])

        projected: list[dict] = []
        for locus in sorted(by_parent.get(parent_id, []), key=lambda row: row["_start"]):
            contained = locus["_start"] >= region_start and locus["_end"] <= region_end
            overlaps = locus["_end"] >= region_start and locus["_start"] <= region_end
            status = "contained" if contained else "crosses_boundary" if overlaps else "outside_region"
            projected_coordinates = (
                f"{locus['_start'] - region_start + 1}-{locus['_end'] - region_start + 1}"
                if contained else ""
            )
            if overlaps:
                projection_rows.append({
                    "sample_id": args.sample_id,
                    "sequence_id": final_id,
                    "parent_sequence_id": parent_id if parent_id != final_id else "",
                    "record_type": region["record_type"],
                    "source_locus_id": locus["locus_id"],
                    "source_coordinates": locus["coordinates"],
                    "projected_coordinates": projected_coordinates,
                    "locus_classification": locus["locus_classification"],
                    "projection_status": status,
                })
            if contained:
                call = {
                    "taxonomy": {column: locus[column] for column in TAXONOMY_COLUMNS},
                    "aggregation_conflict": locus.get("taxonomy_conflict", "").lower() == "true",
                    "best_bitscore": float(locus.get("best_bitscore") or 0),
                }
                projected.append({
                    "row": locus,
                    "call": call,
                    "callers": locus.get("orf_callers", "").split(","),
                    "start": locus["_start"] - region_start + 1,
                    "end": locus["_end"] - region_start + 1,
                })

        counts = Counter(entry["row"]["locus_classification"] for entry in projected)
        viral_entries = [
            entry for entry in projected
            if entry["row"]["locus_classification"] == "viral_supported"
        ]
        clusters = viral_clusters(
            projected, args.cluster_min_viral_loci, args.cluster_max_neutral_gap
        )
        viral_count = counts["viral_supported"]
        cellular_count = counts["cellular_supported"]
        nonviral_class_counts = Counter(
            entry["row"].get("best_nonviral_reference_class", "")
            for entry in projected
            if entry["row"]["locus_classification"] == "cellular_supported"
            and entry["row"].get("best_nonviral_reference_class", "")
        )
        maximum_nonviral_count = max(nonviral_class_counts.values(), default=0)
        dominant_nonviral_candidates = sorted(
            label for label, count in nonviral_class_counts.items()
            if count == maximum_nonviral_count
        )
        if viral_count >= args.cluster_min_viral_loci and cellular_count == 0:
            classification, pattern = "virus", "predominantly_viral"
            reason = "refined_region_has_multiple_viral_loci"
        elif clusters and cellular_count:
            classification = "virus"
            if args.input_type == "dna" and region["record_type"] == "provirus":
                pattern = "resolved_provirus_with_vicat_support"
                reason = "refined_provirus_overlaps_spatial_viral_cluster"
            elif args.input_type == "dna":
                pattern = "localized_viral_cluster"
                reason = "unresolved_parent_retains_potential_provirus_pattern"
            else:
                pattern = "mixed_host_viral_signal"
                reason = "rna_region_contains_mixed_host_viral_signal"
        elif cellular_count:
            classification, pattern = "cellular", "predominantly_cellular"
            reason = "refined_region_lacks_qualified_vicat_viral_cluster"
        else:
            continue

        taxonomy_eligible = classification == "virus" and pattern in {
            "predominantly_viral", "resolved_provirus_with_vicat_support"
        }
        taxonomy_result = aggregate_taxonomy(
            viral_entries, args.contig_taxonomy_support
        ) if taxonomy_eligible else {
            "taxonomy": unclassified_taxonomy(""), "classification_rank": "",
            "taxonomy_support": None, "taxonomy_eligible_loci": 0,
            "taxonomy_supporting_loci": 0, "taxonomy_conflict": False,
        }
        best = max(viral_entries, key=lambda entry: entry["call"]["best_bitscore"]) if viral_entries else None
        cluster_coordinates = [f"{cluster[0]['start']}-{cluster[-1]['end']}" for cluster in clusters]
        flank_statuses = [cluster_flanks(cluster, projected)[0] for cluster in clusters]
        evidence_rows.append({
            "sample_id": args.sample_id,
            "sequence_id": final_id,
            "parent_sequence_id": parent_id if parent_id != final_id else "",
            "record_type": region["record_type"],
            "coordinates": region["coordinates"],
            "tool": "vicat",
            "classification": classification,
            "score": format_number(viral_count / len(projected) if projected else 0),
            "score_type": "projected_viral_supported_locus_fraction",
            "length": region["refined_length"],
            "topology": "",
            "n_genes": len(projected),
            "n_hallmarks": "",
            "evidence_strength": "qualified",
            "strength_basis": f"vicat_projected_{pattern}",
            **taxonomy_result["taxonomy"],
            "orf_loci": len(projected),
            "hit_loci": viral_count,
            "hit_locus_fraction": format_number(viral_count / len(projected) if projected else 0),
            "orf_callers": ",".join(sorted({caller for entry in projected for caller in entry["callers"] if caller})),
            "best_reference_id": best["row"].get("best_reference_id", "") if best else "",
            "best_bitscore": best["row"].get("best_bitscore", "") if best else "",
            "best_evalue": best["row"].get("best_evalue", "") if best else "",
            "best_identity": best["row"].get("best_identity", "") if best else "",
            "best_query_coverage": best["row"].get("best_query_coverage", "") if best else "",
            "best_subject_coverage": best["row"].get("best_subject_coverage", "") if best else "",
            "taxonomy_support": format_number(taxonomy_result["taxonomy_support"]),
            "taxonomy_eligible_loci": taxonomy_result["taxonomy_eligible_loci"],
            "taxonomy_supporting_loci": taxonomy_result["taxonomy_supporting_loci"],
            "classification_rank": taxonomy_result["classification_rank"],
            "taxonomy_conflict": "true" if taxonomy_result["taxonomy_conflict"] else "false",
            "reference_taxonomy_conflict_loci": sum(
                entry["row"].get("reference_taxonomy_conflict", "").lower() == "true"
                for entry in viral_entries
            ),
            "competitive_mode": (
                "true" if any(entry["row"].get("competitive_mode") == "true" for entry in projected)
                else "false"
            ),
            "evidence_scope": "refined_region",
            "origin_pattern": pattern,
            "potential_provirus": "true" if pattern == "localized_viral_cluster" else "false",
            "viral_supported_loci": viral_count,
            "cellular_supported_loci": cellular_count,
            "ambiguous_loci": counts["ambiguous"],
            "uninformative_loci": counts["uninformative"],
            "viral_cluster_count": len(clusters),
            "largest_viral_cluster_loci": max(
                (sum(item["row"]["locus_classification"] == "viral_supported" for item in cluster) for cluster in clusters),
                default=0,
            ),
            "viral_cluster_coordinates": ",".join(cluster_coordinates),
            "viral_cluster_flank_status": ",".join(flank_statuses),
            "competitive_decision_reason": reason,
            "nonviral_supported_classes": ",".join(
                f"{label}:{count}"
                for label, count in sorted(nonviral_class_counts.items())
            ),
            "dominant_nonviral_class": (
                dominant_nonviral_candidates[0]
                if len(dominant_nonviral_candidates) == 1
                else "TIED" if dominant_nonviral_candidates else ""
            ),
        })

    write_tsv(args.output_evidence, OUTPUT_COLUMNS, evidence_rows)
    write_tsv(args.output_projection, PROJECTION_COLUMNS, projection_rows)
    print(
        f"viCAT refinement projection sample={args.sample_id} "
        f"region_evidence={len(evidence_rows)} projected_loci={len(projection_rows)}"
    )


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, ValueError, csv.Error) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
