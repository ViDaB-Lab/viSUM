#!/usr/bin/env python3
"""Convert viCAT ORF alignments into locus-level and contig-level evidence."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import duckdb

from evidence_schema import CORE_EVIDENCE_COLUMNS, TAXONOMY_COLUMNS, unclassified_taxonomy


RANK_NAMES = ("domain", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species")
EXTRA_COLUMNS = [
    "orf_loci", "hit_loci", "hit_locus_fraction", "orf_callers",
    "best_reference_id", "best_bitscore", "best_evalue", "best_identity",
    "best_query_coverage", "best_subject_coverage", "taxonomy_support",
    "classification_rank", "taxonomy_conflict", "reference_taxonomy_conflict_loci",
]
OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + EXTRA_COLUMNS
LOCUS_COLUMNS = [
    "sample_id", "sequence_id", "locus_id", "coordinates", "strand", "orf_ids",
    "orf_callers", "protein_length", "has_qualified_hit", "qualified_hit_count",
    "best_reference_id", "best_bitscore", "best_evalue", "best_identity",
    "best_query_coverage", "best_subject_coverage", "reference_votu_count",
    "reference_taxonomy_coverage", "reference_taxonomy_conflict", "taxonomy_support",
    "classification_rank", "caller_taxonomy_conflict", *TAXONOMY_COLUMNS,
]
AUDIT_COLUMNS = [
    "sample_id", "sequence_id", "locus_id", "orf_id", "orf_callers",
    "coordinates", "strand", "protein_length", "selected_orf_for_locus",
    "selected_best_reference", "reference_id", "identity", "alignment_length",
    "query_length", "subject_length", "query_start", "query_end", "subject_start",
    "subject_end", "evalue", "bitscore", "relative_bitscore", "query_coverage",
    "subject_coverage", "lineage_vote_representative", "lineage_vote_weight",
    "member_votu_count", "classified_votu_count", "reference_taxonomy_coverage",
    "reference_taxonomy_conflict", "reference_taxonomy_conflict_rank",
    "taxonomy_methods", "taxonomy_support_threshold", *TAXONOMY_COLUMNS,
]


@dataclass
class Orf:
    orf_id: str
    sequence_id: str
    start: int
    end: int
    strand: str
    protein_length: int
    callers: set[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--orf-map", required=True, type=Path)
    parser.add_argument("--diamond", required=True, type=Path)
    parser.add_argument("--taxonomy-lookup", required=True, type=Path)
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--taxonomy-support", required=True, type=float)
    parser.add_argument("--locus-overlap", required=True, type=float)
    parser.add_argument("--output-loci", required=True, type=Path)
    parser.add_argument("--output-audit", required=True, type=Path)
    parser.add_argument("--output-evidence", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        if required and not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} is missing columns: {sorted(required - set(reader.fieldnames))}")
        return list(reader)


def load_orfs(path: Path) -> dict[str, Orf]:
    required = {"orf_id", "sequence_id", "start", "end", "strand", "protein_length", "callers"}
    rows = read_tsv(path, required)
    output = {}
    for row in rows:
        orf = Orf(
            row["orf_id"], row["sequence_id"], int(row["start"]), int(row["end"]),
            row["strand"], int(row["protein_length"]), set(row["callers"].split(",")),
        )
        if orf.orf_id in output:
            raise ValueError(f"Duplicate ORF ID: {orf.orf_id}")
        output[orf.orf_id] = orf
    return output


def sql_path(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def diamond_has_rows(path: Path) -> bool:
    with path.open(encoding="utf-8") as handle:
        next(handle, None)
        return next(handle, None) is not None


def load_hits(diamond: Path, lookup: Path) -> dict[str, list[dict]]:
    if not diamond_has_rows(diamond):
        return {}
    rank_sql = ", ".join(f't."{column}" AS "{column}"' for column in TAXONOMY_COLUMNS)
    connection = duckdb.connect()
    missing = connection.execute(
        f"""
        SELECT count(*)
        FROM read_csv({sql_path(diamond)}, delim='\t', header=true, auto_detect=true) d
        LEFT JOIN read_parquet({sql_path(lookup)}) t
          ON d.sseqid = t.representative_protein_id
        WHERE t.representative_protein_id IS NULL
        """
    ).fetchone()[0]
    if missing:
        connection.close()
        raise ValueError(
            f"{missing} DIAMOND alignment(s) reference proteins absent from the viCAT taxonomy lookup"
        )
    query = f"""
        SELECT
            d.qseqid, d.sseqid, d.pident, d.length, d.qlen, d.slen,
            d.qstart, d.qend, d.sstart, d.send, d.evalue, d.bitscore,
            d.qcovhsp, d.scovhsp,
            t.member_votu_count, t.classified_votu_count, t.taxonomy_coverage,
            t.taxonomy_conflict AS reference_taxonomy_conflict,
            t.taxonomy_conflict_rank AS reference_taxonomy_conflict_rank,
            t.taxonomy_methods,
            {rank_sql}
        FROM read_csv(
            {sql_path(diamond)}, delim='\t', header=true, auto_detect=true
        ) d
        INNER JOIN read_parquet({sql_path(lookup)}) t
            ON d.sseqid = t.representative_protein_id
        ORDER BY d.qseqid, d.bitscore DESC, d.sseqid
    """
    relation = connection.execute(query)
    columns = [description[0] for description in relation.description]
    output: dict[str, list[dict]] = defaultdict(list)
    while True:
        batch = relation.fetchmany(10000)
        if not batch:
            break
        for values in batch:
            row = dict(zip(columns, values))
            output[row["qseqid"]].append(row)
    connection.close()
    return output


def is_unclassified(value: object) -> bool:
    if value is None:
        return True
    text = str(value).split("__", 1)[-1].strip().lower()
    return text in {"", "unclassified", "none", "null"}


def taxonomy_call(hits: list[dict], threshold: float) -> dict:
    if not hits:
        return {}
    best = max(hits, key=lambda row: (float(row["bitscore"]), float(row["qcovhsp"])))
    # One vote per distinct full lineage prevents repeated database entries
    # with identical taxonomy from manufacturing support.
    lineages: dict[tuple, dict] = {}
    for hit in hits:
        lineage = tuple(hit[column] for column in TAXONOMY_COLUMNS)
        previous = lineages.get(lineage)
        if previous is None or float(hit["bitscore"]) > float(previous["bitscore"]):
            lineages[lineage] = hit
    best_score = max(float(hit["bitscore"]) for hit in lineages.values())
    weighted = [
        (lineage, float(hit["bitscore"]) / best_score)
        for lineage, hit in lineages.items()
    ]
    total_weight = sum(weight for _, weight in weighted)
    taxonomy = unclassified_taxonomy("virus")
    deepest = "domain"
    deepest_support = 1.0
    aggregation_conflict = False
    accepted_prefix: list[str] = []
    for index, column in enumerate(TAXONOMY_COLUMNS[1:], start=1):
        support: dict[str, float] = defaultdict(float)
        for lineage, weight in weighted:
            if any(lineage[parent_index + 1] != accepted for parent_index, accepted in enumerate(accepted_prefix)):
                continue
            value = lineage[index]
            if not is_unclassified(value):
                support[str(value)] += weight
        if not support:
            break
        winner, winner_weight = max(support.items(), key=lambda item: (item[1], item[0]))
        fraction = winner_weight / total_weight
        if fraction + 1e-12 < threshold:
            aggregation_conflict = len(support) > 1
            break
        taxonomy[column] = winner
        accepted_prefix.append(winner)
        deepest = RANK_NAMES[index]
        deepest_support = fraction

    return {
        "taxonomy": taxonomy,
        "classification_rank": deepest,
        "taxonomy_support": deepest_support,
        "aggregation_conflict": aggregation_conflict,
        "qualified_hit_count": len(hits),
        "best_reference_id": best["sseqid"],
        "best_bitscore": float(best["bitscore"]),
        "best_evalue": float(best["evalue"]),
        "best_identity": float(best["pident"]),
        "best_query_coverage": float(best["qcovhsp"]),
        "best_subject_coverage": float(best["scovhsp"]),
        "reference_votu_count": int(best["member_votu_count"]),
        "reference_taxonomy_coverage": float(best["taxonomy_coverage"]),
        "reference_taxonomy_conflict": bool(best["reference_taxonomy_conflict"]),
    }


def lineage_vote_details(hits: list[dict]) -> tuple[set[int], float]:
    """Identify the best retained hit for each distinct full lineage."""
    lineages: dict[tuple, dict] = {}
    for hit in hits:
        lineage = tuple(hit[column] for column in TAXONOMY_COLUMNS)
        previous = lineages.get(lineage)
        if previous is None or float(hit["bitscore"]) > float(previous["bitscore"]):
            lineages[lineage] = hit
    if not lineages:
        return set(), 0.0
    return {id(hit) for hit in lineages.values()}, max(
        float(hit["bitscore"]) for hit in lineages.values()
    )


class UnionFind:
    def __init__(self, values):
        self.parent = {value: value for value in values}

    def find(self, value):
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left, right):
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def same_locus(left: Orf, right: Orf, threshold: float) -> bool:
    if left.sequence_id != right.sequence_id or left.strand != right.strand:
        return False
    if not left.callers.isdisjoint(right.callers):
        return False
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    if overlap == 0:
        return False
    reciprocal = min(overlap / (left.end - left.start + 1), overlap / (right.end - right.start + 1))
    left_stop = left.end if left.strand == "+" else left.start
    right_stop = right.end if right.strand == "+" else right.start
    return reciprocal >= threshold or left_stop == right_stop


def build_loci(orfs: dict[str, Orf], threshold: float) -> list[list[Orf]]:
    union = UnionFind(orfs)
    by_contig: dict[tuple[str, str], list[Orf]] = defaultdict(list)
    for orf in orfs.values():
        by_contig[(orf.sequence_id, orf.strand)].append(orf)
    for candidates in by_contig.values():
        for index, left in enumerate(candidates):
            for right in candidates[index + 1:]:
                if same_locus(left, right, threshold):
                    union.union(left.orf_id, right.orf_id)
    groups: dict[str, list[Orf]] = defaultdict(list)
    for orf in orfs.values():
        groups[union.find(orf.orf_id)].append(orf)
    return sorted(groups.values(), key=lambda group: (group[0].sequence_id, min(x.start for x in group)))


def first_taxonomy_conflict(left: dict[str, str], right: dict[str, str]) -> int | None:
    for index, column in enumerate(TAXONOMY_COLUMNS[1:], start=1):
        a, b = left[column], right[column]
        if is_unclassified(a) or is_unclassified(b):
            continue
        if a != b:
            return index
    return None


def truncate_taxonomy(taxonomy: dict[str, str], conflict_index: int) -> None:
    defaults = unclassified_taxonomy("virus")
    for column in TAXONOMY_COLUMNS[conflict_index:]:
        taxonomy[column] = defaults[column]


def format_number(value, digits=6):
    if value in (None, ""):
        return ""
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    return f"{float(value):.{digits}g}"


def main() -> None:
    args = parse_args()
    if not 0.5 <= args.taxonomy_support <= 1.0:
        raise SystemExit("--taxonomy-support must be between 0.5 and 1")
    if not 0.0 < args.locus_overlap <= 1.0:
        raise SystemExit("--locus-overlap must be greater than 0 and at most 1")

    orfs = load_orfs(args.orf_map)
    hits = load_hits(args.diamond, args.taxonomy_lookup)
    unknown = sorted(set(hits) - set(orfs))
    if unknown:
        raise SystemExit(f"DIAMOND contains ORF IDs absent from the ORF map: {unknown[:5]}")

    header_rows = read_tsv(args.header_map, {"sequence_id", "length"})
    sequence_lengths = {row["sequence_id"]: int(row["length"]) for row in header_rows}
    loci = build_loci(orfs, args.locus_overlap)
    locus_rows = []
    audit_rows = []
    contig_loci: dict[str, list[dict]] = defaultdict(list)
    for locus_index, alternatives in enumerate(loci, start=1):
        calls = [(orf, taxonomy_call(hits.get(orf.orf_id, []), args.taxonomy_support)) for orf in alternatives]
        hit_calls = [(orf, call) for orf, call in calls if call]
        selected_orf, selected = (None, {})
        if hit_calls:
            selected_orf, selected = max(
                hit_calls,
                key=lambda item: (item[1]["best_bitscore"], item[1]["best_query_coverage"]),
            )
            selected = {**selected, "taxonomy": dict(selected["taxonomy"])}
        caller_conflict = False
        if selected:
            for other_orf, other in hit_calls:
                if other_orf.orf_id == selected_orf.orf_id:
                    continue
                if other["best_bitscore"] < 0.95 * selected["best_bitscore"]:
                    continue
                conflict = first_taxonomy_conflict(selected["taxonomy"], other["taxonomy"])
                if conflict is not None:
                    caller_conflict = True
                    truncate_taxonomy(selected["taxonomy"], conflict)
                    selected["classification_rank"] = RANK_NAMES[conflict - 1]
                    selected["aggregation_conflict"] = True

        sequence_id = alternatives[0].sequence_id
        start, end = min(orf.start for orf in alternatives), max(orf.end for orf in alternatives)
        callers = sorted(set().union(*(orf.callers for orf in alternatives)))
        locus_id = f"{args.sample_id}_vicat_locus{locus_index:09d}"
        row = {
            "sample_id": args.sample_id,
            "sequence_id": sequence_id,
            "locus_id": locus_id,
            "coordinates": f"{start}-{end}",
            "strand": alternatives[0].strand,
            "orf_ids": ",".join(sorted(orf.orf_id for orf in alternatives)),
            "orf_callers": ",".join(callers),
            "protein_length": selected_orf.protein_length if selected_orf else max(orf.protein_length for orf in alternatives),
            "has_qualified_hit": "true" if selected else "false",
            "qualified_hit_count": selected.get("qualified_hit_count", 0),
            "best_reference_id": selected.get("best_reference_id", ""),
            "best_bitscore": format_number(selected.get("best_bitscore")),
            "best_evalue": format_number(selected.get("best_evalue")),
            "best_identity": format_number(selected.get("best_identity")),
            "best_query_coverage": format_number(selected.get("best_query_coverage")),
            "best_subject_coverage": format_number(selected.get("best_subject_coverage")),
            "reference_votu_count": selected.get("reference_votu_count", ""),
            "reference_taxonomy_coverage": format_number(selected.get("reference_taxonomy_coverage")),
            "reference_taxonomy_conflict": "true" if selected.get("reference_taxonomy_conflict") else "false",
            "taxonomy_support": format_number(selected.get("taxonomy_support")),
            "classification_rank": selected.get("classification_rank", ""),
            "caller_taxonomy_conflict": "true" if caller_conflict else "false",
            **(selected.get("taxonomy") or unclassified_taxonomy("")),
        }
        locus_rows.append(row)
        contig_loci[sequence_id].append({"row": row, "call": selected, "callers": callers})

        for orf, call in calls:
            orf_hits = hits.get(orf.orf_id, [])
            vote_representatives, best_lineage_score = lineage_vote_details(orf_hits)
            for hit in orf_hits:
                is_vote = id(hit) in vote_representatives
                bitscore = float(hit["bitscore"])
                audit_rows.append(
                    {
                        "sample_id": args.sample_id,
                        "sequence_id": sequence_id,
                        "locus_id": locus_id,
                        "orf_id": orf.orf_id,
                        "orf_callers": ",".join(sorted(orf.callers)),
                        "coordinates": f"{orf.start}-{orf.end}",
                        "strand": orf.strand,
                        "protein_length": orf.protein_length,
                        "selected_orf_for_locus": "true" if selected_orf is orf else "false",
                        "selected_best_reference": "true" if (
                            selected_orf is orf
                            and hit["sseqid"] == selected.get("best_reference_id")
                            and bitscore == selected.get("best_bitscore")
                        ) else "false",
                        "reference_id": hit["sseqid"],
                        "identity": format_number(hit["pident"]),
                        "alignment_length": hit["length"],
                        "query_length": hit["qlen"],
                        "subject_length": hit["slen"],
                        "query_start": hit["qstart"],
                        "query_end": hit["qend"],
                        "subject_start": hit["sstart"],
                        "subject_end": hit["send"],
                        "evalue": format_number(hit["evalue"]),
                        "bitscore": format_number(bitscore),
                        "relative_bitscore": format_number(
                            bitscore / best_lineage_score if best_lineage_score else None
                        ),
                        "query_coverage": format_number(hit["qcovhsp"]),
                        "subject_coverage": format_number(hit["scovhsp"]),
                        "lineage_vote_representative": "true" if is_vote else "false",
                        "lineage_vote_weight": format_number(
                            bitscore / best_lineage_score if is_vote and best_lineage_score else None
                        ),
                        "member_votu_count": hit["member_votu_count"],
                        "classified_votu_count": hit["classified_votu_count"],
                        "reference_taxonomy_coverage": format_number(hit["taxonomy_coverage"]),
                        "reference_taxonomy_conflict": (
                            "true" if hit["reference_taxonomy_conflict"] else "false"
                        ),
                        "reference_taxonomy_conflict_rank": (
                            hit["reference_taxonomy_conflict_rank"] or ""
                        ),
                        "taxonomy_methods": hit["taxonomy_methods"] or "",
                        "taxonomy_support_threshold": format_number(args.taxonomy_support),
                        **{column: hit[column] or "" for column in TAXONOMY_COLUMNS},
                    }
                )

    evidence_rows = []
    for sequence_id, entries in sorted(contig_loci.items()):
        hit_entries = [entry for entry in entries if entry["call"]]
        if not hit_entries:
            continue
        total_loci, hit_loci = len(entries), len(hit_entries)
        taxonomy = unclassified_taxonomy("virus")
        deepest, deepest_support = "domain", 1.0
        contig_conflict = any(
            entry["call"].get("aggregation_conflict") or
            entry["row"]["caller_taxonomy_conflict"] == "true"
            for entry in hit_entries
        )
        accepted_prefix: list[str] = []
        for index, column in enumerate(TAXONOMY_COLUMNS[1:], start=1):
            counts = Counter(
                entry["call"]["taxonomy"][column]
                for entry in hit_entries
                if all(
                    entry["call"]["taxonomy"][TAXONOMY_COLUMNS[parent_index + 1]] == accepted
                    for parent_index, accepted in enumerate(accepted_prefix)
                )
                if not is_unclassified(entry["call"]["taxonomy"][column])
            )
            if not counts:
                break
            winner, count = max(counts.items(), key=lambda item: (item[1], item[0]))
            support = count / hit_loci
            if support + 1e-12 < args.taxonomy_support:
                contig_conflict = len(counts) > 1 or count < hit_loci
                break
            taxonomy[column] = winner
            accepted_prefix.append(winner)
            deepest, deepest_support = RANK_NAMES[index], support

        best_entry = max(hit_entries, key=lambda entry: entry["call"]["best_bitscore"])
        all_callers = sorted(set().union(*(set(entry["callers"]) for entry in entries)))
        reference_conflicts = sum(
            bool(entry["call"].get("reference_taxonomy_conflict")) for entry in hit_entries
        )
        fraction = hit_loci / total_loci
        evidence_rows.append(
            {
                "sample_id": args.sample_id,
                "sequence_id": sequence_id,
                "parent_sequence_id": "",
                "record_type": "input_contig",
                "coordinates": "",
                "tool": "vicat",
                "classification": "virus",
                "score": format_number(fraction),
                "score_type": "qualified_orf_locus_fraction",
                "length": sequence_lengths.get(sequence_id, ""),
                "topology": "",
                "n_genes": total_loci,
                "n_hallmarks": "",
                **taxonomy,
                "orf_loci": total_loci,
                "hit_loci": hit_loci,
                "hit_locus_fraction": format_number(fraction),
                "orf_callers": ",".join(all_callers),
                "best_reference_id": best_entry["call"]["best_reference_id"],
                "best_bitscore": format_number(best_entry["call"]["best_bitscore"]),
                "best_evalue": format_number(best_entry["call"]["best_evalue"]),
                "best_identity": format_number(best_entry["call"]["best_identity"]),
                "best_query_coverage": format_number(best_entry["call"]["best_query_coverage"]),
                "best_subject_coverage": format_number(best_entry["call"]["best_subject_coverage"]),
                "taxonomy_support": format_number(deepest_support),
                "classification_rank": deepest,
                "taxonomy_conflict": "true" if contig_conflict else "false",
                "reference_taxonomy_conflict_loci": reference_conflicts,
            }
        )

    args.output_loci.parent.mkdir(parents=True, exist_ok=True)
    args.output_audit.parent.mkdir(parents=True, exist_ok=True)
    args.output_evidence.parent.mkdir(parents=True, exist_ok=True)
    with args.output_loci.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOCUS_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(locus_rows)
    with args.output_audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)
    with args.output_evidence.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(evidence_rows)
    print(
        f"viCAT standardized: loci={len(locus_rows)} reference_hits={len(audit_rows)} "
        f"contigs_with_hits={len(evidence_rows)}"
    )


if __name__ == "__main__":
    main()
