#!/usr/bin/env python3

"""Convert selected VITAP assignments into sparse viSUM taxonomy evidence."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

from evidence_schema import CORE_EVIDENCE_COLUMNS, unclassified_taxonomy


VITAP_LINEAGE_RANKS = (
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "realm",
)
RANK_COLUMNS = {
    "realm": "r__Realm",
    "kingdom": "k__Kingdom",
    "phylum": "p__Phylum",
    "class": "c__Class",
    "order": "o__Order",
    "family": "f__Family",
    "genus": "g__Genus",
    "species": "s__Species",
}
RANK_PREFIXES = {
    "realm": "r__",
    "kingdom": "k__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}

OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "vitap_confidence_level",
    "vitap_assignment_method",
    "classification_rank",
]

AUDIT_COLUMNS = [
    "sample_id",
    "sequence_id",
    "in_refinement_map",
    "accepted",
    "decision",
    "raw_lineage",
    "raw_score",
    "vitap_confidence_level",
    "vitap_assignment_method",
    "classification",
    "classification_rank",
    "alternative_lineage_count",
    "uniref90_fallback_taxon",
    "uniref90_fallback_rank",
    "uniref90_fallback_score",
]

MAP_REQUIRED = {
    "sample_id",
    "input_type",
    "sequence_id",
    "parent_sequence_id",
    "record_type",
    "coordinates",
    "refined_length",
}
BEST_REQUIRED = {
    "Genome_ID",
    "lineage",
    "lineage_score/participation_index",
    "Confidence_level",
}
ALL_REQUIRED = {"Genome_ID", "lineage", "lineage_score/participation_index"}
FALLBACK_REQUIRED = {"genome_id", "taxa_name", "participation_index", "taxon_level"}
METADATA_REQUIRED = {
    "sample_id",
    "input_type",
    "refined_sequence_count",
    "best_lineage_row_count",
    "all_lineage_row_count",
    "uniref90_fallback_row_count",
    "database_release",
    "run_status",
}
VMR_REQUIRED = {rank.title() for rank in VITAP_LINEAGE_RANKS}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate and standardize VITAP best-lineage assignments."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--region-map", required=True, type=Path)
    parser.add_argument("--best-lineages", required=True, type=Path)
    parser.add_argument("--all-lineages", required=True, type=Path)
    parser.add_argument("--uniref-fallback", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--vmr", required=True, type=Path)
    parser.add_argument("--output-evidence", required=True, type=Path)
    parser.add_argument("--output-audit", required=True, type=Path)
    return parser.parse_args()


def read_table(path: Path, delimiter: str, required: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Table has no header: {path}")
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return [
            {key.strip(): (value or "").strip() for key, value in row.items() if key}
            for row in reader
        ]


def unique_rows(
    rows: list[dict[str, str]], key: str, label: str
) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = row[key]
        if not identifier:
            raise ValueError(f"{label} contains an empty {key}")
        if identifier in records:
            raise ValueError(f"Duplicate {label} row: {identifier}")
        records[identifier] = row
    return records


def parse_nonnegative_integer(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise ValueError(f"{label} is not an integer: {value}") from error
    if parsed < 0:
        raise ValueError(f"{label} must be nonnegative: {value}")
    return parsed


def parse_score(value: str, label: str) -> float:
    try:
        score = float(value)
    except ValueError as error:
        raise ValueError(f"{label} is not numeric: {value}") from error
    if not math.isfinite(score) or score < 0:
        raise ValueError(f"{label} must be finite and nonnegative: {value}")
    return score


def load_metadata(path: Path, sample_id: str, input_type: str) -> dict[str, str]:
    rows = read_table(path, "\t", METADATA_REQUIRED)
    if len(rows) != 1:
        raise ValueError(f"Expected one VITAP metadata row in {path}; found {len(rows)}")
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError("VITAP metadata does not match the requested sample and input type")
    if row["run_status"] not in {"completed", "skipped_no_refined_candidates"}:
        raise ValueError(f"Unrecognized VITAP run status: {row['run_status']}")
    return row


def load_vmr_taxa(path: Path) -> dict[str, set[str]]:
    rows = read_table(path, ",", VMR_REQUIRED)
    taxa: dict[str, set[str]] = {rank: set() for rank in VITAP_LINEAGE_RANKS}
    for row in rows:
        for rank in VITAP_LINEAGE_RANKS:
            value = row[rank.title()]
            if value and value != "-":
                taxa[rank].add(value)
    if not any(taxa.values()):
        raise ValueError(f"VMR contains no populated taxonomy: {path}")
    return taxa


def parse_lineage(raw: str, vmr_taxa: dict[str, set[str]]) -> dict[str, str]:
    values = [value.strip() for value in raw.split(";")]
    if len(values) != len(VITAP_LINEAGE_RANKS):
        raise ValueError(
            f"VITAP lineage must contain eight ranks (species through realm): {raw}"
        )
    lineage = dict(zip(VITAP_LINEAGE_RANKS, values))
    if all(value in {"", "-"} for value in values):
        raise ValueError(f"VITAP selected an entirely unclassified lineage: {raw}")
    for rank, value in lineage.items():
        if value not in {"", "-"} and value not in vmr_taxa[rank]:
            raise ValueError(f"VITAP taxon is absent from VMR {rank}: {value}")
    return lineage


def assignment_method(confidence: str) -> str:
    return (
        "uniref90_fallback"
        if confidence.casefold() == "uniref90-based".casefold()
        else "vitap_graph_lineage"
    )


def deepest_rank(lineage: dict[str, str]) -> str:
    for rank in VITAP_LINEAGE_RANKS:
        if lineage[rank] not in {"", "-"}:
            return rank
    raise ValueError("Cannot determine a classification rank from an empty lineage")


def is_viriform(lineage: dict[str, str]) -> bool:
    return any(
        "viriform" in value.casefold()
        for value in lineage.values()
        if value not in {"", "-"}
    )


def standardized_taxonomy(
    lineage: dict[str, str], classification: str
) -> dict[str, str]:
    taxonomy = unclassified_taxonomy(classification)
    for rank, column in RANK_COLUMNS.items():
        value = lineage[rank]
        taxonomy[column] = (
            f"{RANK_PREFIXES[rank]}{value}"
            if value not in {"", "-"}
            else f"{RANK_PREFIXES[rank]}unclassified"
        )
    return taxonomy


def format_score(score: float) -> str:
    return f"{score:.15g}"


def validate_fallback(
    sequence_id: str,
    best: dict[str, str],
    lineage: dict[str, str],
    fallback: dict[str, str] | None,
) -> None:
    if fallback is None:
        raise ValueError(f"UniRef90-based assignment lacks a fallback row: {sequence_id}")
    rank = fallback["taxon_level"].casefold()
    if rank not in RANK_COLUMNS:
        raise ValueError(f"Invalid UniRef90 fallback rank for {sequence_id}: {rank}")
    if lineage[rank] != fallback["taxa_name"]:
        raise ValueError(f"UniRef90 fallback taxon disagrees for {sequence_id}")
    best_score = parse_score(best["lineage_score/participation_index"], "VITAP score")
    fallback_score = parse_score(fallback["participation_index"], "fallback score")
    if not math.isclose(best_score, fallback_score, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError(f"UniRef90 fallback score disagrees for {sequence_id}")


def build_audit_row(
    sample_id: str,
    sequence_id: str,
    in_map: bool,
    best: dict[str, str] | None,
    decision: str,
    alternative_count: int,
    fallback: dict[str, str] | None,
    classification: str = "",
    rank: str = "",
) -> dict[str, str]:
    confidence = best["Confidence_level"] if best else ""
    return {
        "sample_id": sample_id,
        "sequence_id": sequence_id,
        "in_refinement_map": str(in_map).lower(),
        "accepted": str(decision == "accepted_best_assignment").lower(),
        "decision": decision,
        "raw_lineage": best["lineage"] if best else "",
        "raw_score": best["lineage_score/participation_index"] if best else "",
        "vitap_confidence_level": confidence,
        "vitap_assignment_method": assignment_method(confidence) if confidence else "",
        "classification": classification,
        "classification_rank": rank,
        "alternative_lineage_count": str(alternative_count),
        "uniref90_fallback_taxon": fallback["taxa_name"] if fallback else "",
        "uniref90_fallback_rank": fallback["taxon_level"] if fallback else "",
        "uniref90_fallback_score": fallback["participation_index"] if fallback else "",
    }


def run(args: argparse.Namespace) -> None:
    metadata = load_metadata(args.run_metadata, args.sample_id, args.input_type)
    map_rows = unique_rows(
        read_table(args.region_map, "\t", MAP_REQUIRED), "sequence_id", "refinement map"
    )
    for sequence_id, row in map_rows.items():
        if row["sample_id"] != args.sample_id or row["input_type"] != args.input_type:
            raise ValueError(f"Refinement-map sample/type mismatch: {sequence_id}")
        parse_nonnegative_integer(row["refined_length"], f"refined length for {sequence_id}")

    best_list = read_table(args.best_lineages, "\t", BEST_REQUIRED)
    best_rows = unique_rows(best_list, "Genome_ID", "VITAP best-lineage")
    all_rows = read_table(args.all_lineages, "\t", ALL_REQUIRED)
    fallback_list = read_table(args.uniref_fallback, "\t", FALLBACK_REQUIRED)
    fallback_rows = unique_rows(fallback_list, "genome_id", "VITAP fallback")
    vmr_taxa = load_vmr_taxa(args.vmr)

    expected_counts = {
        "refined_sequence_count": len(map_rows),
        "best_lineage_row_count": len(best_list),
        "all_lineage_row_count": len(all_rows),
        "uniref90_fallback_row_count": len(fallback_list),
    }
    for field, observed in expected_counts.items():
        if parse_nonnegative_integer(metadata[field], field.replace("_", " ")) != observed:
            raise ValueError(f"VITAP metadata count disagrees for {field}")
    expected_status = "skipped_no_refined_candidates" if not map_rows else "completed"
    if metadata["run_status"] != expected_status:
        raise ValueError("VITAP run status disagrees with the refinement map")

    alternative_counts = Counter(row["Genome_ID"] for row in all_rows)
    evidence_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []

    for sequence_id, map_row in map_rows.items():
        best = best_rows.get(sequence_id)
        fallback = fallback_rows.get(sequence_id)
        if best is None:
            audit_rows.append(
                build_audit_row(
                    args.sample_id,
                    sequence_id,
                    True,
                    None,
                    "no_best_assignment",
                    alternative_counts[sequence_id],
                    fallback,
                )
            )
            continue

        score = parse_score(best["lineage_score/participation_index"], "VITAP score")
        lineage = parse_lineage(best["lineage"], vmr_taxa)
        method = assignment_method(best["Confidence_level"])
        if method == "uniref90_fallback":
            validate_fallback(sequence_id, best, lineage, fallback)
        classification = "viriform" if is_viriform(lineage) else "virus"
        rank = deepest_rank(lineage)
        evidence_rows.append(
            {
                "sample_id": args.sample_id,
                "sequence_id": sequence_id,
                "parent_sequence_id": map_row["parent_sequence_id"],
                "record_type": map_row["record_type"],
                "coordinates": map_row["coordinates"],
                "tool": "vitap",
                "classification": classification,
                "score": format_score(score),
                "score_type": (
                    "vitap_uniref90_participation_index"
                    if method == "uniref90_fallback"
                    else "vitap_lineage_score"
                ),
                "length": map_row["refined_length"],
                "topology": "",
                "n_genes": "",
                "n_hallmarks": "",
                **standardized_taxonomy(lineage, classification),
                "vitap_confidence_level": best["Confidence_level"],
                "vitap_assignment_method": method,
                "classification_rank": rank,
            }
        )
        audit_rows.append(
            build_audit_row(
                args.sample_id,
                sequence_id,
                True,
                best,
                "accepted_best_assignment",
                alternative_counts[sequence_id],
                fallback,
                classification,
                rank,
            )
        )

    reference_ids = (
        set(best_rows)
        | {row["Genome_ID"] for row in all_rows}
        | set(fallback_rows)
    ).difference(map_rows)
    for sequence_id in sorted(reference_ids):
        audit_rows.append(
            build_audit_row(
                args.sample_id,
                sequence_id,
                False,
                best_rows.get(sequence_id),
                "excluded_non_input_reference",
                alternative_counts[sequence_id],
                fallback_rows.get(sequence_id),
            )
        )

    args.output_evidence.parent.mkdir(parents=True, exist_ok=True)
    args.output_audit.parent.mkdir(parents=True, exist_ok=True)
    with args.output_evidence.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(evidence_rows)
    with args.output_audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=AUDIT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(audit_rows)

    print(
        f"Standardized VITAP sample={args.sample_id} evidence={len(evidence_rows)} "
        f"unassigned={len(map_rows) - len(evidence_rows)} "
        f"excluded_references={len(reference_ids)}"
    )


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
