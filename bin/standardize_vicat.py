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

# A single DNA ORF is allowed to enter discovery only when the viral alignment
# is both extensive and decisively better than any nonviral competitor. These
# values are deliberately stricter than the ordinary multi-locus filters.
STRICT_SINGLE_LOCUS_MIN_BITSCORE = 100.0
STRICT_SINGLE_LOCUS_MIN_QUERY_COVERAGE = 70.0
STRICT_SINGLE_LOCUS_MIN_SUBJECT_COVERAGE = 70.0
STRICT_SINGLE_LOCUS_MIN_COMPETITIVE_MARGIN = 0.30
EXTRA_COLUMNS = [
    "orf_loci", "hit_loci", "hit_locus_fraction", "orf_callers",
    "best_reference_id", "best_bitscore", "best_evalue", "best_identity",
    "best_query_coverage", "best_subject_coverage", "taxonomy_support",
    "taxonomy_eligible_loci", "taxonomy_supporting_loci",
    "classification_rank", "taxonomy_conflict", "reference_taxonomy_conflict_loci",
    "competitive_mode", "evidence_scope", "origin_pattern", "potential_provirus",
    "viral_supported_loci", "cellular_supported_loci", "ambiguous_loci",
    "uninformative_loci", "viral_cluster_count", "largest_viral_cluster_loci",
    "viral_cluster_coordinates", "viral_cluster_flank_status",
    "competitive_decision_reason", "nonviral_supported_classes",
    "dominant_nonviral_class", "cluster_min_viral_loci",
    "embedded_context_min_viral_loci",
]
OUTPUT_COLUMNS = (
    CORE_EVIDENCE_COLUMNS
    + ["evidence_strength", "strength_basis"]
    + EXTRA_COLUMNS
)
LOCUS_COLUMNS = [
    "sample_id", "sequence_id", "locus_id", "coordinates", "strand", "orf_ids",
    "orf_callers", "protein_length", "has_qualified_hit", "qualified_hit_count",
    "best_reference_id", "best_bitscore", "best_evalue", "best_identity",
    "best_query_coverage", "best_subject_coverage", "reference_votu_count",
    "reference_taxonomy_coverage", "reference_taxonomy_conflict", "taxonomy_support",
    "classification_rank", "caller_taxonomy_conflict", "locus_classification",
    "best_viral_bitscore", "best_cellular_bitscore", "competitive_score_margin",
    "competitive_decision_reason", "best_nonviral_reference_class",
    "competitive_mode", *TAXONOMY_COLUMNS,
]
AUDIT_COLUMNS = [
    "sample_id", "sequence_id", "locus_id", "orf_id", "orf_callers",
    "coordinates", "strand", "protein_length", "selected_orf_for_locus",
    "selected_best_reference", "reference_id", "reference_class",
    "nonviral_reference_class", "source_classes", "replicon_accessions",
    "replicon_types", "provirus_flank_eligible", "classification_note",
    "cluster_member_count",
    "source_protein_id", "cellular_group", "source_accession", "organism_name",
    "taxid", "identity", "alignment_length",
    "query_length", "subject_length", "query_start", "query_end", "subject_start",
    "subject_end", "evalue", "bitscore", "relative_bitscore", "query_coverage",
    "subject_coverage", "lineage_vote_representative", "lineage_vote_weight",
    "member_votu_count", "classified_votu_count", "reference_taxonomy_coverage",
    "reference_taxonomy_conflict", "reference_taxonomy_conflict_rank",
    "taxonomy_methods", "orf_taxonomy_support_threshold", "locus_classification",
    "best_viral_bitscore", "best_cellular_bitscore", "competitive_score_margin",
    "competitive_decision_reason", *TAXONOMY_COLUMNS,
]
CLUSTER_COLUMNS = [
    "sample_id", "sequence_id", "cluster_id", "coordinates",
    "viral_supported_loci", "intervening_neutral_loci", "flank_status",
    "left_cellular_loci", "right_cellular_loci", "classification_rank",
    "taxonomy_support", "taxonomy_eligible_loci", "taxonomy_supporting_loci",
    "taxonomy_conflict", *TAXONOMY_COLUMNS,
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
    parser.add_argument("--diamond", type=Path)
    parser.add_argument("--taxonomy-lookup", type=Path)
    parser.add_argument("--reference-manifest", type=Path)
    parser.add_argument("--prepared-hits", type=Path)
    parser.add_argument("--prepared-metadata", type=Path)
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--orf-taxonomy-support", required=True, type=float)
    parser.add_argument("--contig-taxonomy-support", required=True, type=float)
    parser.add_argument("--locus-overlap", required=True, type=float)
    parser.add_argument("--competitive-min-margin", type=float, default=0.05)
    parser.add_argument("--cluster-min-viral-loci", type=int, default=2)
    parser.add_argument("--cluster-max-neutral-gap", type=int, default=1)
    parser.add_argument(
        "--dna-single-locus-rescue",
        choices=("off", "strict"),
        default="off",
    )
    parser.add_argument("--output-loci", required=True, type=Path)
    parser.add_argument("--output-clusters", type=Path)
    parser.add_argument("--output-context", type=Path)
    parser.add_argument("--output-provirus-evidence", type=Path)
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


def load_hits(
    diamond: Path,
    lookup: Path,
    reference_manifest: Path | None = None,
    threads: int = 1,
) -> dict[str, list[dict]]:
    if not diamond_has_rows(diamond):
        return {}
    rank_sql = ", ".join(f't."{column}" AS "{column}"' for column in TAXONOMY_COLUMNS)
    connection = duckdb.connect()
    connection.execute(f"SET threads = {threads}")
    if reference_manifest is not None:
        missing = connection.execute(
            f"""
            SELECT count(*)
            FROM read_csv({sql_path(diamond)}, delim='\t', header=true, auto_detect=true) d
            LEFT JOIN read_parquet({sql_path(reference_manifest)}) m
              ON d.sseqid = m.reference_id
            WHERE m.reference_id IS NULL
            """
        ).fetchone()[0]
        if missing:
            connection.close()
            raise ValueError(
                f"{missing} DIAMOND alignment(s) absent from the competitive reference manifest"
            )
        query = f"""
            SELECT
                d.qseqid, d.sseqid, d.pident, d.length, d.qlen, d.slen,
                d.qstart, d.qend, d.sstart, d.send, d.evalue, d.bitscore,
                d.qcovhsp, d.scovhsp,
                m.reference_class, m.source_protein_id, m.cellular_group,
                m.source_accession, m.organism_name, m.taxid,
                coalesce(t.member_votu_count, 0) AS member_votu_count,
                coalesce(t.classified_votu_count, 0) AS classified_votu_count,
                coalesce(t.taxonomy_coverage, 0) AS taxonomy_coverage,
                coalesce(t.taxonomy_conflict, false) AS reference_taxonomy_conflict,
                coalesce(t.taxonomy_conflict_rank, '') AS reference_taxonomy_conflict_rank,
                coalesce(t.taxonomy_methods, '') AS taxonomy_methods,
                {rank_sql}
            FROM read_csv(
                {sql_path(diamond)}, delim='\t', header=true, auto_detect=true
            ) d
            INNER JOIN read_parquet({sql_path(reference_manifest)}) m
              ON d.sseqid = m.reference_id
            LEFT JOIN read_parquet({sql_path(lookup)}) t
              ON m.reference_class = 'viral'
             AND m.source_protein_id = t.representative_protein_id
            ORDER BY d.qseqid, d.bitscore DESC, d.sseqid
        """
    else:
        query = f"""
            SELECT
                d.qseqid, d.sseqid, d.pident, d.length, d.qlen, d.slen,
                d.qstart, d.qend, d.sstart, d.send, d.evalue, d.bitscore,
                d.qcovhsp, d.scovhsp,
                'viral' AS reference_class, d.sseqid AS source_protein_id,
                '' AS cellular_group, '' AS source_accession,
                '' AS organism_name, '' AS taxid,
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
    if reference_manifest is not None:
        missing_viral = connection.execute(
            f"""
            SELECT count(*)
            FROM (
                SELECT DISTINCT m.source_protein_id
                FROM read_csv(
                    {sql_path(diamond)}, delim='\t', header=true, auto_detect=true
                ) d
                INNER JOIN read_parquet({sql_path(reference_manifest)}) m
                  ON d.sseqid = m.reference_id
                WHERE m.reference_class = 'viral'
            ) m
            LEFT JOIN read_parquet({sql_path(lookup)}) t
              ON m.source_protein_id = t.representative_protein_id
            WHERE t.representative_protein_id IS NULL
            """
        ).fetchone()[0]
        if missing_viral:
            connection.close()
            raise ValueError(
                f"{missing_viral} viral manifest reference(s) absent from the viCAT taxonomy lookup"
            )
    else:
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


def load_prepared_hits(
    prepared_hits: Path, threads: int
) -> dict[str, list[dict]]:
    connection = duckdb.connect()
    connection.execute(f"SET threads = {threads}")
    relation = connection.execute(
        f"""
        SELECT * EXCLUDE (manifest_matched, taxonomy_matched, competitive_mode)
        FROM read_parquet({sql_path(prepared_hits)})
        """
    )
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
    for rows in output.values():
        rows.sort(
            key=lambda row: (
                -float(row["bitscore"]),
                str(row["sseqid"]),
            )
        )
    return output


def prepared_competitive_mode(path: Path) -> bool:
    rows = read_tsv(path, {"competitive_mode"})
    if len(rows) != 1 or rows[0]["competitive_mode"] not in {"true", "false"}:
        raise ValueError(f"Invalid viCAT hit-preparation metadata: {path}")
    return rows[0]["competitive_mode"] == "true"


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
    taxonomy = unclassified_taxonomy("virus")
    deepest = "domain"
    deepest_support = 1.0
    aggregation_conflict = False
    accepted_prefix: list[str] = []
    for index, column in enumerate(TAXONOMY_COLUMNS[1:], start=1):
        support: dict[str, float] = defaultdict(float)
        eligible_weight = 0.0
        for lineage, weight in weighted:
            if any(lineage[parent_index + 1] != accepted for parent_index, accepted in enumerate(accepted_prefix)):
                continue
            value = lineage[index]
            if not is_unclassified(value):
                support[str(value)] += weight
                eligible_weight += weight
        if not support:
            break
        winner, winner_weight = max(support.items(), key=lambda item: (item[1], item[0]))
        # A reference whose database taxonomy was truncated at this rank is an
        # abstention here. It remains a qualified viral-homology hit, but must
        # contribute neither support nor opposition below its last safe rank.
        fraction = winner_weight / eligible_weight
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


def best_class_hit(alternatives: list[Orf], hits: dict[str, list[dict]], label: str):
    candidates = [
        (orf, hit)
        for orf in alternatives
        for hit in hits.get(orf.orf_id, [])
        if hit["reference_class"] == label
    ]
    if not candidates:
        return None, None
    return max(
        candidates,
        key=lambda item: (float(item[1]["bitscore"]), float(item[1]["qcovhsp"])),
    )


def classify_locus(
    alternatives: list[Orf],
    hits: dict[str, list[dict]],
    competitive_mode: bool,
    minimum_margin: float,
) -> dict[str, object]:
    viral_orf, viral_hit = best_class_hit(alternatives, hits, "viral")
    cellular_orf, cellular_hit = best_class_hit(alternatives, hits, "cellular")
    viral_score = float(viral_hit["bitscore"]) if viral_hit else None
    cellular_score = float(cellular_hit["bitscore"]) if cellular_hit else None

    if not competitive_mode:
        classification = "viral_supported" if viral_hit else "uninformative"
        reason = "viral_only_database_hit" if viral_hit else "no_qualified_hit"
        margin = None
    elif viral_score is None and cellular_score is None:
        classification, reason, margin = "uninformative", "no_qualified_hit", None
    elif cellular_score is None:
        classification, reason, margin = "viral_supported", "viral_hit_without_returned_cellular_competitor", 1.0
    elif viral_score is None:
        classification, reason, margin = "cellular_supported", "cellular_hit_without_returned_viral_competitor", -1.0
    else:
        signed_margin = (viral_score - cellular_score) / max(viral_score, cellular_score)
        margin = signed_margin
        if abs(signed_margin) + 1e-12 < minimum_margin:
            classification, reason = "ambiguous", "viral_cellular_scores_within_margin"
        elif signed_margin > 0:
            classification, reason = "viral_supported", "viral_score_exceeds_cellular_margin"
        else:
            classification, reason = "cellular_supported", "cellular_score_exceeds_viral_margin"

    return {
        "classification": classification,
        "reason": reason,
        "margin": margin,
        "viral_orf": viral_orf,
        "viral_hit": viral_hit,
        "viral_score": viral_score,
        "cellular_orf": cellular_orf,
        "cellular_hit": cellular_hit,
        "cellular_score": cellular_score,
        "nonviral_reference_class": (
            str(cellular_hit.get("nonviral_reference_class") or "CELLULAR")
            if cellular_hit else ""
        ),
    }


def qualifies_strict_single_locus_rescue(entry: dict) -> bool:
    """Return whether one viral-supported DNA locus merits guarded discovery."""
    row = entry["row"]
    try:
        bitscore = float(row.get("best_viral_bitscore", ""))
        query_coverage = float(row.get("best_query_coverage", ""))
        subject_coverage = float(row.get("best_subject_coverage", ""))
        margin = float(row.get("competitive_score_margin", ""))
    except (TypeError, ValueError):
        return False
    return (
        bitscore >= STRICT_SINGLE_LOCUS_MIN_BITSCORE
        and query_coverage >= STRICT_SINGLE_LOCUS_MIN_QUERY_COVERAGE
        and subject_coverage >= STRICT_SINGLE_LOCUS_MIN_SUBJECT_COVERAGE
        and margin >= STRICT_SINGLE_LOCUS_MIN_COMPETITIVE_MARGIN
    )


def viral_clusters(entries: list[dict], minimum_loci: int, maximum_neutral_gap: int) -> list[list[dict]]:
    """Return viral-locus clusters; cellular-supported loci always break a cluster."""
    clusters: list[list[dict]] = []
    current: list[dict] = []
    neutral: list[dict] = []

    def finish() -> None:
        nonlocal current, neutral
        viral_count = sum(
            entry["row"]["locus_classification"] == "viral_supported"
            for entry in current
        )
        if viral_count >= minimum_loci:
            clusters.append(current)
        current, neutral = [], []

    for entry in sorted(entries, key=lambda item: item["start"]):
        classification = entry["row"]["locus_classification"]
        if classification == "cellular_supported":
            finish()
        elif classification == "viral_supported":
            if current and neutral:
                current.extend(neutral)
            current.append(entry)
            neutral = []
        elif current:
            neutral.append(entry)
            if len(neutral) > maximum_neutral_gap:
                finish()
    finish()
    return clusters


def cluster_flanks(cluster: list[dict], entries: list[dict]) -> tuple[str, int, int]:
    start, end = cluster[0]["start"], cluster[-1]["end"]
    def flank_eligible(entry: dict) -> bool:
        # Legacy competitive databases had no class-aware field and their
        # generic cellular references remain eligible for compatibility.
        return entry["row"].get("best_nonviral_reference_class", "") in {
            "", "CELLULAR", "CELLULAR_CHROMOSOME"
        }

    left = sum(
        entry["row"]["locus_classification"] == "cellular_supported"
        and flank_eligible(entry)
        and entry["end"] < start
        for entry in entries
    )
    right = sum(
        entry["row"]["locus_classification"] == "cellular_supported"
        and flank_eligible(entry)
        and entry["start"] > end
        for entry in entries
    )
    if left and right:
        status = "both_sides"
    elif left:
        status = "left_only"
    elif right:
        status = "right_only"
    else:
        status = "none"
    return status, left, right


def aggregate_taxonomy(entries: list[dict], threshold: float) -> dict[str, object]:
    taxonomy = unclassified_taxonomy("virus")
    deepest, deepest_support = "domain", 1.0
    deepest_eligible_loci = len(entries)
    deepest_supporting_loci = len(entries)
    conflict = any(
        entry["call"].get("aggregation_conflict")
        or entry["row"]["caller_taxonomy_conflict"] == "true"
        for entry in entries
    )
    accepted_prefix: list[str] = []
    for index, column in enumerate(TAXONOMY_COLUMNS[1:], start=1):
        eligible = [
            entry for entry in entries
            if all(
                entry["call"]["taxonomy"][TAXONOMY_COLUMNS[parent + 1]] == accepted
                for parent, accepted in enumerate(accepted_prefix)
            )
            if not is_unclassified(entry["call"]["taxonomy"][column])
        ]
        counts = Counter(entry["call"]["taxonomy"][column] for entry in eligible)
        if not counts:
            break
        winner, count = max(counts.items(), key=lambda item: (item[1], item[0]))
        support = count / len(eligible)
        if support + 1e-12 < threshold:
            conflict = len(counts) > 1
            break
        taxonomy[column] = winner
        accepted_prefix.append(winner)
        deepest, deepest_support = RANK_NAMES[index], support
        deepest_eligible_loci, deepest_supporting_loci = len(eligible), count
    return {
        "taxonomy": taxonomy,
        "classification_rank": deepest,
        "taxonomy_support": deepest_support,
        "taxonomy_eligible_loci": deepest_eligible_loci,
        "taxonomy_supporting_loci": deepest_supporting_loci,
        "taxonomy_conflict": conflict,
    }


def main() -> None:
    args = parse_args()
    if not 0.5 <= args.orf_taxonomy_support <= 1.0:
        raise SystemExit("--orf-taxonomy-support must be between 0.5 and 1")
    if not 0.5 <= args.contig_taxonomy_support <= 1.0:
        raise SystemExit("--contig-taxonomy-support must be between 0.5 and 1")
    if not 0.0 < args.locus_overlap <= 1.0:
        raise SystemExit("--locus-overlap must be greater than 0 and at most 1")
    if not 0.0 <= args.competitive_min_margin < 1.0:
        raise SystemExit("--competitive-min-margin must be between 0 and 1")
    if args.cluster_min_viral_loci < 1:
        raise SystemExit("--cluster-min-viral-loci must be at least 1")
    if args.cluster_max_neutral_gap < 0:
        raise SystemExit("--cluster-max-neutral-gap must be zero or greater")
    if args.threads < 1:
        raise SystemExit("--threads must be a positive integer")

    prepared_mode = args.prepared_hits is not None or args.prepared_metadata is not None
    if prepared_mode:
        if args.prepared_hits is None or args.prepared_metadata is None:
            raise SystemExit(
                "--prepared-hits and --prepared-metadata must be provided together"
            )
        if any((args.diamond, args.taxonomy_lookup, args.reference_manifest)):
            raise SystemExit(
                "Prepared viCAT hits cannot be combined with raw DIAMOND/database inputs"
            )
    elif args.diamond is None or args.taxonomy_lookup is None:
        raise SystemExit(
            "Provide --prepared-hits/--prepared-metadata or "
            "--diamond/--taxonomy-lookup"
        )

    orfs = load_orfs(args.orf_map)
    if prepared_mode:
        hits = load_prepared_hits(args.prepared_hits, args.threads)
        competitive_mode = prepared_competitive_mode(args.prepared_metadata)
    else:
        hits = load_hits(
            args.diamond, args.taxonomy_lookup, args.reference_manifest, args.threads
        )
        competitive_mode = args.reference_manifest is not None
    unknown = sorted(set(hits) - set(orfs))
    if unknown:
        raise SystemExit(f"DIAMOND contains ORF IDs absent from the ORF map: {unknown[:5]}")

    header_rows = read_tsv(args.header_map, {"sequence_id", "length"})
    sequence_lengths = {row["sequence_id"]: int(row["length"]) for row in header_rows}
    loci = build_loci(orfs, args.locus_overlap)
    locus_rows = []
    args.output_audit.parent.mkdir(parents=True, exist_ok=True)
    audit_handle = args.output_audit.open("w", encoding="utf-8", newline="")
    audit_writer = csv.DictWriter(
        audit_handle, fieldnames=AUDIT_COLUMNS, delimiter="\t", lineterminator="\n"
    )
    audit_writer.writeheader()
    audit_count = 0
    contig_loci: dict[str, list[dict]] = defaultdict(list)
    for locus_index, alternatives in enumerate(loci, start=1):
        competition = classify_locus(
            alternatives, hits, competitive_mode, args.competitive_min_margin
        )
        calls = [
            (
                orf,
                taxonomy_call(
                    [
                        hit for hit in hits.get(orf.orf_id, [])
                        if hit["reference_class"] == "viral"
                    ],
                    args.orf_taxonomy_support,
                ),
            )
            for orf in alternatives
        ]
        hit_calls = [(orf, call) for orf, call in calls if call]
        selected_orf, selected = (None, {})
        if competition["classification"] == "viral_supported":
            selected_orf = competition["viral_orf"]
            selected = next(
                (call for orf, call in hit_calls if orf is selected_orf), {}
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
            "locus_classification": competition["classification"],
            "best_viral_bitscore": format_number(competition["viral_score"]),
            "best_cellular_bitscore": format_number(competition["cellular_score"]),
            "competitive_score_margin": format_number(competition["margin"]),
            "competitive_decision_reason": competition["reason"],
            "best_nonviral_reference_class": competition["nonviral_reference_class"],
            "competitive_mode": "true" if competitive_mode else "false",
            **(selected.get("taxonomy") or unclassified_taxonomy("")),
        }
        locus_rows.append(row)
        contig_loci[sequence_id].append(
            {
                "row": row,
                "call": selected,
                "callers": callers,
                "start": start,
                "end": end,
            }
        )

        for orf, call in calls:
            orf_hits = hits.get(orf.orf_id, [])
            viral_orf_hits = [
                hit for hit in orf_hits if hit["reference_class"] == "viral"
            ]
            vote_representatives, best_lineage_score = lineage_vote_details(
                viral_orf_hits
            )
            for hit in orf_hits:
                is_vote = id(hit) in vote_representatives
                bitscore = float(hit["bitscore"])
                audit_writer.writerow(
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
                        "reference_class": hit["reference_class"],
                        "nonviral_reference_class": hit.get("nonviral_reference_class", "") or "",
                        "source_classes": hit.get("source_classes", "") or "",
                        "replicon_accessions": hit.get("replicon_accessions", "") or "",
                        "replicon_types": hit.get("replicon_types", "") or "",
                        "provirus_flank_eligible": str(
                            bool(hit.get("provirus_flank_eligible", False))
                        ).lower(),
                        "classification_note": hit.get("classification_note", "") or "",
                        "cluster_member_count": hit.get("cluster_member_count", "") or "",
                        "source_protein_id": hit["source_protein_id"] or "",
                        "cellular_group": hit["cellular_group"] or "",
                        "source_accession": hit["source_accession"] or "",
                        "organism_name": hit["organism_name"] or "",
                        "taxid": hit["taxid"] or "",
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
                        "orf_taxonomy_support_threshold": format_number(args.orf_taxonomy_support),
                        "locus_classification": competition["classification"],
                        "best_viral_bitscore": format_number(competition["viral_score"]),
                        "best_cellular_bitscore": format_number(competition["cellular_score"]),
                        "competitive_score_margin": format_number(competition["margin"]),
                        "competitive_decision_reason": competition["reason"],
                        **{column: hit[column] or "" for column in TAXONOMY_COLUMNS},
                    }
                )
                audit_count += 1

    audit_handle.close()

    evidence_rows = []
    context_rows = []
    cluster_rows = []
    provirus_rows = []
    for sequence_id, entries in sorted(contig_loci.items()):
        total_loci = len(entries)
        counts = Counter(entry["row"]["locus_classification"] for entry in entries)
        viral_entries = [
            entry for entry in entries
            if entry["row"]["locus_classification"] == "viral_supported"
        ]
        clusters = viral_clusters(
            entries, args.cluster_min_viral_loci, args.cluster_max_neutral_gap
        )
        embedded_context_min = max(2, args.cluster_min_viral_loci)
        embedded_clusters = [
            cluster for cluster in clusters
            if sum(
                item["row"]["locus_classification"] == "viral_supported"
                for item in cluster
            ) >= embedded_context_min
        ]
        viral_count = counts["viral_supported"]
        cellular_count = counts["cellular_supported"]
        ambiguous_count = counts["ambiguous"]
        uninformative_count = counts["uninformative"]
        nonviral_class_counts = Counter(
            entry["row"].get("best_nonviral_reference_class", "")
            for entry in entries
            if entry["row"]["locus_classification"] == "cellular_supported"
            and entry["row"].get("best_nonviral_reference_class", "")
        )
        nonviral_supported_classes = ",".join(
            f"{label}:{count}"
            for label, count in sorted(nonviral_class_counts.items())
        )
        if nonviral_class_counts:
            maximum_nonviral_count = max(nonviral_class_counts.values())
            dominant_candidates = sorted(
                label for label, count in nonviral_class_counts.items()
                if count == maximum_nonviral_count
            )
            dominant_nonviral_class = (
                dominant_candidates[0] if len(dominant_candidates) == 1 else "TIED"
            )
        else:
            dominant_nonviral_class = ""

        strict_single_locus_rescue = (
            args.input_type == "dna"
            and args.dna_single_locus_rescue == "strict"
            and viral_count == 1
            and cellular_count == 0
            and qualifies_strict_single_locus_rescue(viral_entries[0])
        )

        if viral_count >= args.cluster_min_viral_loci and cellular_count == 0:
            pattern = "predominantly_viral"
            classification = "virus"
            decision_reason = "viral_locus_threshold_met_without_cellular_supported_loci"
        elif strict_single_locus_rescue:
            pattern = "single_locus_viral_rescue"
            classification = "virus"
            decision_reason = (
                "strict_single_locus_alignment_thresholds_met_without_"
                "vicat_cellular_supported_loci"
            )
        elif embedded_clusters and cellular_count:
            pattern = (
                "localized_viral_cluster"
                if args.input_type == "dna"
                else "mixed_host_viral_signal"
            )
            classification = "virus"
            decision_reason = "spatial_viral_cluster_with_cellular_context"
        elif viral_count == 1:
            pattern = "isolated_viral_locus"
            classification = "cellular" if cellular_count else ""
            decision_reason = "isolated_viral_locus_is_not_qualified_discovery_evidence"
        elif viral_count >= 2:
            pattern = "mixed_distributed"
            classification = "cellular" if cellular_count else ""
            decision_reason = "viral_loci_do_not_form_a_qualified_cluster"
        elif cellular_count:
            pattern = "predominantly_cellular"
            classification = "cellular"
            decision_reason = "cellular_supported_loci_without_qualified_viral_pattern"
        else:
            pattern = "unresolved"
            classification = ""
            decision_reason = "no_qualified_competitive_origin_pattern"

        cluster_coordinates = [
            f"{cluster[0]['start']}-{cluster[-1]['end']}" for cluster in clusters
        ]
        flank_statuses = [cluster_flanks(cluster, entries)[0] for cluster in clusters]

        # Weak viCAT patterns do not enter the Discovery Gate, but a separate
        # audit-only context table allows viHARMONY to preserve their history.
        if not classification:
            context_rows.append({
                "sample_id": args.sample_id,
                "sequence_id": sequence_id,
                "parent_sequence_id": "",
                "record_type": "input_contig",
                "coordinates": "",
                "tool": "vicat_context",
                "classification": "ambiguous",
                "score": format_number(viral_count / total_loci),
                "score_type": "viral_supported_locus_fraction",
                "length": sequence_lengths.get(sequence_id, ""),
                "topology": "",
                "n_genes": total_loci,
                "n_hallmarks": "",
                "evidence_strength": "weak",
                "strength_basis": f"vicat_audit_only_{pattern}",
                **unclassified_taxonomy(""),
                "orf_loci": total_loci,
                "hit_loci": viral_count,
                "hit_locus_fraction": format_number(viral_count / total_loci),
                "competitive_mode": "true" if competitive_mode else "false",
                "evidence_scope": "parent_discovery",
                "origin_pattern": pattern,
                "potential_provirus": "false",
                "viral_supported_loci": viral_count,
                "cellular_supported_loci": cellular_count,
                "ambiguous_loci": ambiguous_count,
                "uninformative_loci": uninformative_count,
                "viral_cluster_count": len(clusters),
                "largest_viral_cluster_loci": 0,
                "viral_cluster_coordinates": ",".join(cluster_coordinates),
                "viral_cluster_flank_status": ",".join(flank_statuses),
                "competitive_decision_reason": decision_reason,
                "nonviral_supported_classes": nonviral_supported_classes,
                "dominant_nonviral_class": dominant_nonviral_class,
                "cluster_min_viral_loci": args.cluster_min_viral_loci,
                "embedded_context_min_viral_loci": embedded_context_min,
            })
            continue

        taxonomy_entries = (
            viral_entries
            if pattern in {"predominantly_viral", "single_locus_viral_rescue"}
            else []
        )
        taxonomy_result = aggregate_taxonomy(
            taxonomy_entries, args.contig_taxonomy_support
        ) if taxonomy_entries else {
            "taxonomy": unclassified_taxonomy("virus" if classification == "virus" else ""),
            "classification_rank": "domain" if classification == "virus" else "",
            "taxonomy_support": 1.0 if classification == "virus" else None,
            "taxonomy_eligible_loci": 0,
            "taxonomy_supporting_loci": 0,
            "taxonomy_conflict": False,
        }

        best_entry = (
            max(viral_entries, key=lambda entry: entry["call"]["best_bitscore"])
            if viral_entries else None
        )
        all_callers = sorted(set().union(*(set(entry["callers"]) for entry in entries)))
        reference_conflicts = sum(
            bool(entry["call"].get("reference_taxonomy_conflict")) for entry in viral_entries
        )
        fraction = viral_count / total_loci
        flank_statuses = []
        for cluster_index, cluster in enumerate(clusters, start=1):
            cluster_viral_entries = [
                entry for entry in cluster
                if entry["row"]["locus_classification"] == "viral_supported"
            ]
            cluster_taxonomy = aggregate_taxonomy(
                cluster_viral_entries, args.contig_taxonomy_support
            )
            flank_status, left_cellular, right_cellular = cluster_flanks(
                cluster, entries
            )
            flank_statuses.append(flank_status)
            cluster_rows.append({
                "sample_id": args.sample_id,
                "sequence_id": sequence_id,
                "cluster_id": f"{sequence_id}|vicat_cluster_{cluster_index}",
                "coordinates": f"{cluster[0]['start']}-{cluster[-1]['end']}",
                "viral_supported_loci": len(cluster_viral_entries),
                "intervening_neutral_loci": len(cluster) - len(cluster_viral_entries),
                "flank_status": flank_status,
                "left_cellular_loci": left_cellular,
                "right_cellular_loci": right_cellular,
                "classification_rank": cluster_taxonomy["classification_rank"],
                "taxonomy_support": format_number(cluster_taxonomy["taxonomy_support"]),
                "taxonomy_eligible_loci": cluster_taxonomy["taxonomy_eligible_loci"],
                "taxonomy_supporting_loci": cluster_taxonomy["taxonomy_supporting_loci"],
                "taxonomy_conflict": "true" if cluster_taxonomy["taxonomy_conflict"] else "false",
                **cluster_taxonomy["taxonomy"],
            })
            if (
                args.input_type == "dna"
                and cellular_count
                and cluster in embedded_clusters
            ):
                cluster_start, cluster_end = cluster[0]["start"], cluster[-1]["end"]
                provirus_rows.append({
                    "sample_id": args.sample_id,
                    "sequence_id": f"{sequence_id}|vicat_provirus_{cluster_start}_{cluster_end}",
                    "parent_sequence_id": sequence_id,
                    "record_type": "provirus",
                    "coordinates": f"{cluster_start}-{cluster_end}",
                    "tool": "vicat",
                    "classification": "virus",
                    "score": format_number(len(cluster_viral_entries) / len(cluster)),
                    "score_type": "viral_supported_locus_fraction",
                    "length": cluster_end - cluster_start + 1,
                    "topology": "",
                    "n_genes": len(cluster),
                    "n_hallmarks": "",
                    "evidence_strength": "qualified",
                    "strength_basis": "vicat_advisory_provirus_boundary",
                    **cluster_taxonomy["taxonomy"],
                    "orf_loci": len(cluster),
                    "hit_loci": len(cluster_viral_entries),
                    "hit_locus_fraction": format_number(
                        len(cluster_viral_entries) / len(cluster)
                    ),
                    "taxonomy_support": format_number(cluster_taxonomy["taxonomy_support"]),
                    "taxonomy_eligible_loci": cluster_taxonomy["taxonomy_eligible_loci"],
                    "taxonomy_supporting_loci": cluster_taxonomy["taxonomy_supporting_loci"],
                    "classification_rank": cluster_taxonomy["classification_rank"],
                    "taxonomy_conflict": (
                        "true" if cluster_taxonomy["taxonomy_conflict"] else "false"
                    ),
                    "competitive_mode": "true" if competitive_mode else "false",
                    "evidence_scope": "provirus_boundary_advisory",
                    "origin_pattern": "localized_viral_cluster",
                    "potential_provirus": "true",
                    "viral_supported_loci": len(cluster_viral_entries),
                    "cellular_supported_loci": cellular_count,
                    "ambiguous_loci": sum(
                        item["row"]["locus_classification"] == "ambiguous"
                        for item in cluster
                    ),
                    "uninformative_loci": sum(
                        item["row"]["locus_classification"] == "uninformative"
                        for item in cluster
                    ),
                    "viral_cluster_count": 1,
                    "largest_viral_cluster_loci": len(cluster_viral_entries),
                    "viral_cluster_coordinates": f"{cluster_start}-{cluster_end}",
                    "viral_cluster_flank_status": flank_status,
                    "competitive_decision_reason": "spatial_viral_cluster_with_nonviral_context",
                    "nonviral_supported_classes": nonviral_supported_classes,
                    "dominant_nonviral_class": dominant_nonviral_class,
                    "cluster_min_viral_loci": args.cluster_min_viral_loci,
                    "embedded_context_min_viral_loci": embedded_context_min,
                })
        largest_cluster = max(
            (
                sum(item["row"]["locus_classification"] == "viral_supported" for item in cluster)
                for cluster in clusters
            ),
            default=0,
        )
        evidence_rows.append(
            {
                "sample_id": args.sample_id,
                "sequence_id": sequence_id,
                "parent_sequence_id": "",
                "record_type": "input_contig",
                "coordinates": "",
                "tool": "vicat",
                "classification": classification,
                "score": format_number(fraction),
                "score_type": "viral_supported_locus_fraction",
                "length": sequence_lengths.get(sequence_id, ""),
                "topology": "",
                "n_genes": total_loci,
                "n_hallmarks": "",
                "evidence_strength": "qualified",
                "strength_basis": (
                    "vicat_strict_single_locus_rescue"
                    if pattern == "single_locus_viral_rescue"
                    else (
                        "vicat_viral_protein_homology"
                        if not competitive_mode
                        else f"vicat_competitive_{pattern}"
                    )
                ),
                **taxonomy_result["taxonomy"],
                "orf_loci": total_loci,
                "hit_loci": viral_count,
                "hit_locus_fraction": format_number(fraction),
                "orf_callers": ",".join(all_callers),
                "best_reference_id": best_entry["call"]["best_reference_id"] if best_entry else "",
                "best_bitscore": format_number(best_entry["call"]["best_bitscore"] if best_entry else None),
                "best_evalue": format_number(best_entry["call"]["best_evalue"] if best_entry else None),
                "best_identity": format_number(best_entry["call"]["best_identity"] if best_entry else None),
                "best_query_coverage": format_number(best_entry["call"]["best_query_coverage"] if best_entry else None),
                "best_subject_coverage": format_number(best_entry["call"]["best_subject_coverage"] if best_entry else None),
                "taxonomy_support": format_number(taxonomy_result["taxonomy_support"]),
                "taxonomy_eligible_loci": taxonomy_result["taxonomy_eligible_loci"],
                "taxonomy_supporting_loci": taxonomy_result["taxonomy_supporting_loci"],
                "classification_rank": taxonomy_result["classification_rank"],
                "taxonomy_conflict": "true" if taxonomy_result["taxonomy_conflict"] else "false",
                "reference_taxonomy_conflict_loci": reference_conflicts,
                "competitive_mode": "true" if competitive_mode else "false",
                "evidence_scope": "parent_discovery",
                "origin_pattern": pattern,
                "potential_provirus": "true" if pattern == "localized_viral_cluster" else "false",
                "viral_supported_loci": viral_count,
                "cellular_supported_loci": cellular_count,
                "ambiguous_loci": ambiguous_count,
                "uninformative_loci": uninformative_count,
                "viral_cluster_count": len(clusters),
                "largest_viral_cluster_loci": largest_cluster,
                "viral_cluster_coordinates": ",".join(cluster_coordinates),
                "viral_cluster_flank_status": ",".join(flank_statuses),
                "competitive_decision_reason": decision_reason,
                "nonviral_supported_classes": nonviral_supported_classes,
                "dominant_nonviral_class": dominant_nonviral_class,
                "cluster_min_viral_loci": args.cluster_min_viral_loci,
                "embedded_context_min_viral_loci": embedded_context_min,
            }
        )

    args.output_loci.parent.mkdir(parents=True, exist_ok=True)
    args.output_evidence.parent.mkdir(parents=True, exist_ok=True)
    with args.output_loci.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOCUS_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(locus_rows)
    if args.output_clusters is not None:
        args.output_clusters.parent.mkdir(parents=True, exist_ok=True)
        with args.output_clusters.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=CLUSTER_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(cluster_rows)
    if args.output_context is not None:
        args.output_context.parent.mkdir(parents=True, exist_ok=True)
        with args.output_context.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(context_rows)
    if args.output_provirus_evidence is not None:
        args.output_provirus_evidence.parent.mkdir(parents=True, exist_ok=True)
        with args.output_provirus_evidence.open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(provirus_rows)
    with args.output_evidence.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(evidence_rows)
    print(
        f"viCAT standardized: loci={len(locus_rows)} reference_hits={audit_count} "
        f"contigs_with_hits={len(evidence_rows)}"
    )


if __name__ == "__main__":
    main()
