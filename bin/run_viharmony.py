#!/usr/bin/env python3

"""Combine viSUM evidence into final sequences, decisions, and taxonomy calls."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO


RANKS = ("domain", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species")
RANK_COLUMNS = {
    "domain": "d__Domain", "realm": "r__Realm", "kingdom": "k__Kingdom",
    "phylum": "p__Phylum", "class": "c__Class", "order": "o__Order",
    "family": "f__Family", "genus": "g__Genus", "species": "s__Species",
}
PREFIXES = {rank: column[:3] for rank, column in RANK_COLUMNS.items()}
UNCLASSIFIED = {rank: f"{PREFIXES[rank]}unclassified" for rank in RANKS}

ORIGIN_TOOLS = {
    "genomad", "virsorter2", "cenotetaker3", "deep6", "deepmicroclass2",
    "virbot", "gianthunter", "vicat", "checkv",
}
TAXONOMY_TOOLS = {
    "genomad", "cenotetaker3", "virbot", "gianthunter", "vicat", "vitap", "vcontact3",
}
METHOD_FAMILY = {
    "genomad": "marker_gene", "cenotetaker3": "hallmark_gene",
    "virbot": "protein_model", "gianthunter": "protein_model",
    "vicat": "protein_homology", "vitap": "network_homology",
    "vcontact3": "gene_sharing_network",
}
TIER_SCORE = {"provisional": 1, "moderate": 2, "high": 3}

METADATA_COLUMNS = [
    "original_contig_name", "normalized_name", "final_sequence_id", "record_type",
    "provirus_coordinates", "original_length", "refined_length", "viral_decision",
    "viral_confidence", "strong_tools", "qualified_tools", "cellular_conflict_tools",
    "plasmid_conflict_tools", "strict_taxonomy", "strict_taxonomy_rank",
    "analysis_taxonomy", "analysis_taxonomy_rank", "taxonomy_confidence",
    "taxonomy_supporting_tools", "taxonomy_conflict", "vcontact3_groups",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the viHARMONY decision engine.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", choices=("dna", "rna"), required=True)
    parser.add_argument("--normalized-fasta", required=True, type=Path)
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--discovery-gate", required=True, type=Path)
    parser.add_argument("--refined-fasta", required=True, type=Path)
    parser.add_argument("--region-map", required=True, type=Path)
    parser.add_argument("--ictv-msl", required=True, type=Path)
    parser.add_argument("--evidence", nargs="*", default=[], type=Path)
    parser.add_argument("--vcontact3-groups", nargs="*", default=[], type=Path)
    parser.add_argument("--audit-mode", choices=("none", "compact", "full"), default="compact")
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Table has no header: {path}")
        return list(reader.fieldnames), list(reader)


def require(columns: Iterable[str], required: set[str], path: Path) -> None:
    missing = required.difference(columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(parts)
                identifier = line[1:].split(maxsplit=1)[0]
                if not identifier or identifier in records:
                    raise ValueError(f"Invalid or duplicate FASTA identifier at {path}:{line_number}")
                parts = []
            elif identifier is None:
                raise ValueError(f"Sequence precedes FASTA header at {path}:{line_number}")
            else:
                parts.append("".join(line.split()))
    if identifier is not None:
        records[identifier] = "".join(parts)
    if any(not sequence for sequence in records.values()):
        raise ValueError(f"Empty FASTA sequence in {path}")
    return records


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_taxon(value: str, rank: str) -> str:
    value = (value or "").strip()
    prefix = PREFIXES[rank]
    if value.startswith(prefix):
        value = value[len(prefix):]
    if not value or value.lower() in {"unclassified", "na", "nan", "none", "-"}:
        return ""
    return value


def load_msl(path: Path) -> tuple[set[tuple[str, ...]], str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"ICTV MSL has no header: {path}")
        normalized = {name.strip().lower(): name for name in reader.fieldnames}
        required = [rank for rank in RANKS if rank != "domain"]
        missing = [rank for rank in required if rank not in normalized]
        if missing:
            raise ValueError(f"ICTV MSL lacks rank columns: {missing}")
        lineages = {
            tuple((row.get(normalized[rank]) or "").strip() for rank in required)
            for row in reader
        }
    if not lineages:
        raise ValueError(f"ICTV MSL contains no taxonomy rows: {path}")
    release = path.stem
    return lineages, release


def lineage_is_valid(lineage: dict[str, str], msl: set[tuple[str, ...]]) -> bool:
    named = [(rank, lineage.get(rank, "")) for rank in RANKS[1:] if lineage.get(rank)]
    if not named:
        return False
    rank_positions = {rank: index for index, rank in enumerate(RANKS[1:])}
    return any(all(msl_row[rank_positions[rank]] == taxon for rank, taxon in named) for msl_row in msl)


def vote_tier(row: dict[str, str]) -> str:
    tool = row.get("tool", "").strip().lower()
    if tool == "vitap":
        confidence = row.get("vitap_confidence_level", "").lower()
        if "high" in confidence:
            return "high"
        if "medium" in confidence or "uniref" in confidence:
            return "moderate"
        return "provisional"
    if tool == "vcontact3":
        return "moderate" if row.get("vcontact3_assignment_method") == "reference_supported" else "provisional"
    if tool == "vicat":
        try:
            loci = int(float(row.get("taxonomy_supporting_loci", "0") or 0))
        except ValueError:
            loci = 0
        return "moderate" if loci >= 2 and row.get("taxonomy_conflict", "").lower() not in {"true", "1", "yes"} else "provisional"
    if tool in {"cenotetaker3", "virbot"}:
        try:
            support = int(float(row.get("n_hallmarks", "0") or row.get("n_genes", "0") or 0))
        except ValueError:
            support = 0
        return "moderate" if support >= 2 else "provisional"
    if tool == "gianthunter":
        return "moderate" if "lca" in row.get("score_type", "").lower() else "provisional"
    if tool == "genomad":
        return "moderate"
    return "provisional"


def root_for_evidence(row: dict[str, str], final_to_parent: dict[str, str]) -> str:
    sequence_id = row.get("sequence_id", "").strip()
    parent = row.get("parent_sequence_id", "").strip()
    if sequence_id in final_to_parent:
        return sequence_id
    target_parent = parent or sequence_id
    matches = [final_id for final_id, final_parent in final_to_parent.items() if final_parent == target_parent]
    return matches[0] if len(matches) == 1 else target_parent


def taxonomy_string(lineage: dict[str, str]) -> str:
    return ";".join(f"{PREFIXES[rank]}{lineage.get(rank) or 'unclassified'}" for rank in RANKS)


def decide_taxonomy(rows: list[dict[str, str]], msl: set[tuple[str, ...]]) -> tuple[dict[str, str], str, str, list[str], list[dict[str, object]]]:
    votes: list[dict[str, object]] = []
    for row in rows:
        tool = row.get("tool", "").strip().lower()
        if tool not in TAXONOMY_TOOLS or row.get("classification") != "virus":
            continue
        lineage = {rank: clean_taxon(row.get(column, ""), rank) for rank, column in RANK_COLUMNS.items()}
        lineage["domain"] = "Viruses"
        valid = lineage_is_valid(lineage, msl)
        tier = vote_tier(row)
        for rank in RANKS[1:]:
            taxon = lineage.get(rank, "")
            if taxon:
                votes.append({"tool": tool, "method_family": METHOD_FAMILY[tool], "rank": rank, "taxon": taxon, "tier": tier, "valid": valid, "accepted": False, "reason": "pending" if valid else "invalid_configured_ictv_msl", "_lineage": lineage})

    accepted = {rank: "" for rank in RANKS}
    accepted["domain"] = "Viruses"
    supporting_tools: set[str] = set()
    conflict = "false"
    deepest = "domain"
    for rank in RANKS[1:]:
        earlier_ranks = RANKS[1:RANKS.index(rank)]
        rank_votes = [
            vote for vote in votes
            if vote["rank"] == rank
            and vote["valid"]
            and all(
                not accepted[ancestor]
                or vote["_lineage"].get(ancestor) == accepted[ancestor]
                for ancestor in earlier_ranks
            )
        ]
        if not rank_votes:
            break
        grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
        for vote in rank_votes:
            grouped[str(vote["taxon"])].append(vote)
        scores = {
            taxon: (
                max(TIER_SCORE[str(vote["tier"])] for vote in taxon_votes),
                len({str(vote["method_family"]) for vote in taxon_votes}),
                len({str(vote["tool"]) for vote in taxon_votes}),
            )
            for taxon, taxon_votes in grouped.items()
        }
        best_score = max(scores.values())
        winners = [taxon for taxon, score in scores.items() if score == best_score]
        if len(winners) != 1:
            conflict = "true"
            for vote in rank_votes:
                vote["reason"] = "rank_tie"
            break
        winner = winners[0]
        accepted[rank] = winner
        deepest = rank
        for vote in rank_votes:
            if vote["taxon"] == winner:
                vote["accepted"] = True
                vote["reason"] = "selected"
                supporting_tools.add(str(vote["tool"]))
            else:
                vote["reason"] = "lower_method_aware_support"

    deepest_votes = [vote for vote in votes if vote["rank"] == deepest and vote["accepted"]]
    families = {str(vote["method_family"]) for vote in deepest_votes}
    max_tier = max((TIER_SCORE[str(vote["tier"])] for vote in deepest_votes), default=0)
    if len(families) >= 2:
        confidence = "multi_method_consensus"
    elif max_tier == 3:
        confidence = "single_method_high"
    elif deepest == "domain":
        confidence = "unclassified"
    else:
        confidence = "provisional"
    return accepted, deepest, confidence, sorted(supporting_tools), votes


def original_output_id(original_id: str, record_type: str, coordinates: str) -> str:
    if record_type in {"provirus", "viral_region"}:
        return f"{original_id}|provirus_{coordinates.replace('-', '_')}"
    return original_id


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(args: argparse.Namespace) -> None:
    normalized = read_fasta(args.normalized_fasta)
    refined = read_fasta(args.refined_fasta)
    header_columns, header_rows = read_tsv(args.header_map)
    require(header_columns, {"sample_id", "sequence_id", "original_id", "original_header", "length"}, args.header_map)
    headers = {row["sequence_id"]: row for row in header_rows}
    gate_columns, gate_rows = read_tsv(args.discovery_gate)
    require(gate_columns, {"sequence_id", "advance_to_refinement", "discovery_status"}, args.discovery_gate)
    gate = {row["sequence_id"]: row for row in gate_rows}
    region_columns, region_rows = read_tsv(args.region_map)
    require(region_columns, {"sequence_id", "parent_sequence_id", "record_type", "coordinates", "original_length", "refined_length"}, args.region_map)
    regions = {row["sequence_id"]: row for row in region_rows}
    if set(refined) != set(regions):
        raise ValueError("Refined FASTA and provirus-region map contain different sequence IDs")

    final_to_parent = {sequence_id: row["parent_sequence_id"].strip() or sequence_id for sequence_id, row in regions.items()}
    evidence_rows: list[dict[str, str]] = []
    evidence_columns: list[str] = []
    for path in args.evidence:
        columns, rows = read_tsv(path)
        require(columns, {"sample_id", "sequence_id", "parent_sequence_id", "tool", "classification"}, path)
        evidence_columns.extend(column for column in columns if column not in evidence_columns)
        for row in rows:
            if row["sample_id"] != args.sample_id:
                raise ValueError(f"Evidence sample mismatch in {path}: {row['sample_id']}")
            row["_source_file"] = path.name
            evidence_rows.append(row)

    groups_by_sequence: dict[str, list[dict[str, str]]] = defaultdict(list)
    for path in args.vcontact3_groups:
        columns, rows = read_tsv(path)
        if "sequence_id" not in columns:
            continue
        for row in rows:
            groups_by_sequence[row["sequence_id"]].append(row)

    msl, msl_release = load_msl(args.ictv_msl)
    rows_by_final: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in evidence_rows:
        exact = row.get("sequence_id", "").strip()
        parent = row.get("parent_sequence_id", "").strip() or exact
        for final_id, final_parent in final_to_parent.items():
            # Region-specific refinement calls must never leak into a sibling
            # region from the same parent. Parent-contig discovery evidence is
            # inherited by every refined child from that parent.
            if exact in final_to_parent:
                applies = exact == final_id
            else:
                applies = parent == final_parent or exact == final_parent
            if applies:
                rows_by_final[final_id].append(row)

    metadata_rows: list[dict[str, object]] = []
    review_rows: list[dict[str, object]] = []
    evidence_audit_rows: list[dict[str, object]] = []
    taxonomy_audit_rows: list[dict[str, object]] = []
    normalized_fasta_rows: list[tuple[str, str]] = []
    original_fasta_rows: list[tuple[str, str]] = []
    used_original_output_ids: set[str] = set()
    map_rows: list[dict[str, object]] = []

    for final_id, sequence in refined.items():
        region = regions[final_id]
        parent_id = final_to_parent[final_id]
        header = headers[parent_id]
        rows = rows_by_final.get(final_id, [])
        tools_by_class: dict[str, set[str]] = defaultdict(set)
        strong_tools: set[str] = set()
        qualified_tools: set[str] = set()
        for row in rows:
            tool = row.get("tool", "").strip().lower()
            classification = row.get("classification", "").strip().lower()
            if tool in ORIGIN_TOOLS:
                tools_by_class[classification].add(tool)
                if classification == "virus":
                    if row.get("evidence_strength", "qualified").strip().lower() == "strong":
                        strong_tools.add(tool)
                    else:
                        qualified_tools.add(tool)
            audit = dict(row)
            audit["final_sequence_id"] = final_id
            audit["applied_scope"] = "region" if row.get("sequence_id") == final_id else "parent"
            evidence_audit_rows.append(audit)

        qualified_tools.difference_update(strong_tools)
        if len(strong_tools) >= 2:
            viral_confidence = "high"
        elif strong_tools or len(qualified_tools) >= 2:
            viral_confidence = "supported"
        else:
            viral_confidence = "provisional"

        strict, strict_rank, tax_confidence, tax_tools, vote_rows = decide_taxonomy(rows, msl)
        for vote in vote_rows:
            vote["sample_id"] = args.sample_id
            vote["final_sequence_id"] = final_id
            taxonomy_audit_rows.append(vote)

        analysis = dict(strict)
        group_labels: list[str] = []
        for group in groups_by_sequence.get(final_id, []):
            for rank in RANKS[1:-1]:
                prediction = (group.get(f"{rank}_prediction") or "").strip()
                if prediction.startswith(("novel_", "unplaced_")):
                    group_labels.append(f"{rank}:{prediction}")
                    if not analysis.get(rank):
                        analysis[rank] = f"viharmony_{args.sample_id}_{prediction}"

        coords = region["coordinates"].strip() or "NA"
        record_type = region["record_type"].strip()
        strict_taxonomy = taxonomy_string(strict)
        analysis_taxonomy = taxonomy_string(analysis)
        original_name = header["original_id"]
        cellular = sorted(tools_by_class.get("cellular", set()))
        plasmid = sorted(tools_by_class.get("plasmid", set()))
        metadata = {
            "original_contig_name": original_name,
            "normalized_name": parent_id,
            "final_sequence_id": final_id,
            "record_type": record_type,
            "provirus_coordinates": coords,
            "original_length": region["original_length"],
            "refined_length": len(sequence),
            "viral_decision": "retained_viral_candidate",
            "viral_confidence": viral_confidence,
            "strong_tools": ",".join(sorted(strong_tools)) or "NA",
            "qualified_tools": ",".join(sorted(qualified_tools)) or "NA",
            "cellular_conflict_tools": ",".join(cellular) or "NA",
            "plasmid_conflict_tools": ",".join(plasmid) or "NA",
            "strict_taxonomy": strict_taxonomy,
            "strict_taxonomy_rank": strict_rank,
            "analysis_taxonomy": analysis_taxonomy,
            "analysis_taxonomy_rank": next((rank for rank in reversed(RANKS) if analysis.get(rank)), "domain"),
            "taxonomy_confidence": tax_confidence,
            "taxonomy_supporting_tools": ",".join(tax_tools) or "NA",
            "taxonomy_conflict": "true" if any(vote["reason"] == "rank_tie" for vote in vote_rows) else "false",
            "vcontact3_groups": ",".join(sorted(set(group_labels))) or "NA",
        }
        metadata_rows.append(metadata)
        reasons = []
        if viral_confidence == "provisional": reasons.append("single_qualified_viral_tool")
        if cellular: reasons.append("cellular_conflict")
        if plasmid: reasons.append("plasmid_conflict")
        if metadata["taxonomy_conflict"] == "true": reasons.append("taxonomy_conflict")
        if strict_rank == "domain": reasons.append("taxonomy_unclassified")
        if reasons:
            review_rows.append({"final_sequence_id": final_id, "review_reasons": ",".join(reasons), **metadata})

        normalized_fasta_rows.append((final_id, sequence))
        original_output = original_output_id(original_name, record_type, coords)
        if original_output in used_original_output_ids:
            original_output = f"{original_output}|visum_{parent_id}"
        used_original_output_ids.add(original_output)
        original_fasta_rows.append((original_output, sequence))
        map_rows.append({
            "original_contig_name": original_name, "original_header": header["original_header"],
            "normalized_parent_id": parent_id, "final_sequence_id": final_id,
            "final_original_id": original_output, "record_type": record_type,
            "provirus_coordinates": coords,
        })

    prefix = args.output_prefix
    outputs = {
        "normalized_fasta": Path(f"{prefix}.final.normalized.fasta"),
        "original_fasta": Path(f"{prefix}.final.original_ids.fasta"),
        "metadata": Path(f"{prefix}.final_metadata.tsv"),
        "review": Path(f"{prefix}.review_queue.tsv"),
        "disposition": Path(f"{prefix}.sequence_disposition.tsv"),
        "map": Path(f"{prefix}.sequence_map.tsv"),
        "database_fasta": Path(f"{prefix}.database_candidates.fasta"),
        "database_tsv": Path(f"{prefix}.database_candidates.tsv"),
    }
    write_fasta(outputs["normalized_fasta"], normalized_fasta_rows)
    write_fasta(outputs["original_fasta"], original_fasta_rows)
    write_tsv(outputs["metadata"], METADATA_COLUMNS, metadata_rows)
    write_tsv(outputs["review"], ["review_reasons", *METADATA_COLUMNS], review_rows)
    write_tsv(outputs["map"], ["original_contig_name", "original_header", "normalized_parent_id", "final_sequence_id", "final_original_id", "record_type", "provirus_coordinates"], map_rows)
    write_fasta(outputs["database_fasta"], normalized_fasta_rows)
    write_tsv(outputs["database_tsv"], METADATA_COLUMNS, metadata_rows)

    final_by_parent: dict[str, list[str]] = defaultdict(list)
    for row in map_rows:
        final_by_parent[str(row["normalized_parent_id"])].append(str(row["final_sequence_id"]))
    disposition_rows = []
    for sequence_id in normalized:
        gate_row = gate[sequence_id]
        final_ids = final_by_parent.get(sequence_id, [])
        if final_ids:
            disposition = "retained_as_provirus" if any("|viral_region_" in item for item in final_ids) else "retained_contig"
        elif gate_row["advance_to_refinement"].strip().lower() == "true":
            disposition = "rejected_after_refinement"
        elif gate_row["discovery_status"] == "likely_nonviral":
            disposition = "discovery_noncandidate_conflicting_origin"
        else:
            disposition = "discovery_noncandidate_insufficient_support"
        disposition_rows.append({"original_contig_name": headers[sequence_id]["original_id"], "normalized_name": sequence_id, "disposition": disposition, "final_sequence_ids": ",".join(final_ids) or "NA"})
    write_tsv(outputs["disposition"], ["original_contig_name", "normalized_name", "disposition", "final_sequence_ids"], disposition_rows)

    if args.audit_mode == "full":
        evidence_path = Path(f"{prefix}.evidence_audit.tsv.gz")
        audit_columns = ["final_sequence_id", "applied_scope", "_source_file", *evidence_columns]
        with gzip.open(evidence_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=audit_columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader(); writer.writerows(evidence_audit_rows)
        taxonomy_path = Path(f"{prefix}.taxonomy_vote_audit.tsv.gz")
        tax_columns = ["sample_id", "final_sequence_id", "tool", "method_family", "rank", "taxon", "tier", "valid", "accepted", "reason"]
        with gzip.open(taxonomy_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=tax_columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader(); writer.writerows(taxonomy_audit_rows)
        outputs["evidence_audit"] = evidence_path
        outputs["taxonomy_audit"] = taxonomy_path

    manifest_path = Path(f"{prefix}.harmonizer_manifest.json")
    manifest = {
        "schema_version": "viharmony-0.1",
        "sample_id": args.sample_id,
        "input_type": args.input_type,
        "ictv_msl": args.ictv_msl.name,
        "ictv_msl_release": msl_release,
        "audit_mode": args.audit_mode,
        "tools_present": sorted({row.get("tool", "") for row in evidence_rows if row.get("tool")}),
        "policy": {"viral_confidence": "high=2+ strong; supported=1 strong or 2+ qualified; provisional=1 qualified", "taxonomy": "configured-ICTV-MSL-validated method-aware rank voting"},
        "outputs": {},
    }
    for label, path in outputs.items():
        manifest["outputs"][label] = {"file": path.name, "sha256": sha256(path)}
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"viHARMONY sample={args.sample_id} retained={len(metadata_rows)} review={len(review_rows)}")


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, ValueError, csv.Error) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
