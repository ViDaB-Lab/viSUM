#!/usr/bin/env python3

"""Convert VirSorter2 calls into sparse, standardized viSUM evidence."""

import argparse
import csv
import re
from pathlib import Path

from evidence_schema import CORE_EVIDENCE_COLUMNS, unclassified_taxonomy


OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "call_type",
    "parent_length",
    "viral_gene_percent",
    "cellular_gene_percent",
    "max_score_group",
    "evidence_strength",
    "strength_basis",
]

SCORE_REQUIRED = {
    "seqname",
    "max_score",
    "max_score_group",
    "length",
    "hallmark",
    "viral",
    "cellular",
}

BOUNDARY_REQUIRED = {
    "seqname",
    "trim_orf_index_start",
    "trim_orf_index_end",
    "trim_bp_start",
    "trim_bp_end",
    "partial",
    "hallmark_cnt",
    "shape",
    "seqname_new",
    "final_max_score",
    "final_max_score_group",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert VirSorter2 outputs into standardized viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--score-table", required=True, type=Path)
    parser.add_argument("--boundary-table", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path, required_columns: set[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        missing = required_columns.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return reader.fieldnames, list(reader)


def clean_missing(value: str | None) -> str:
    value = "" if value is None else value.strip()
    return "" if value.lower() in {"", "na", "nan", "none"} else value


def load_header_map(path: Path, sample_id: str) -> dict[str, dict[str, str]]:
    _, rows = read_tsv(path, {"sample_id", "sequence_id", "length"})
    records: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["sample_id"] != sample_id:
            raise ValueError(
                f"Header-map sample '{row['sample_id']}' does not match '{sample_id}'"
            )
        sequence_id = row["sequence_id"].strip()
        if not sequence_id:
            raise ValueError("Header map contains an empty sequence_id")
        if sequence_id in records:
            raise ValueError(f"Duplicate sequence_id in header map: {sequence_id}")
        parse_positive_int(row["length"], "header-map length", sequence_id)
        records[sequence_id] = row
    if not records:
        raise ValueError("Header map contains no sequence records")
    return records


def load_run_metadata(path: Path, sample_id: str, input_type: str) -> dict[str, str]:
    _, rows = read_tsv(
        path,
        {
            "sample_id",
            "input_type",
            "virsorter2_version",
            "classifier_groups",
            "min_length",
            "min_score",
            "run_status",
            "virus_call_count",
        },
    )
    if len(rows) != 1:
        raise ValueError(
            f"Expected one VirSorter2 metadata row in {path}; found {len(rows)}"
        )
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError(
            "VirSorter2 metadata does not match the requested sample and input type"
        )
    if not clean_missing(row["virsorter2_version"]):
        raise ValueError("VirSorter2 metadata contains an empty version")
    parse_positive_int(row["min_length"], "minimum length", sample_id)
    parse_score(row["min_score"], "minimum score", sample_id)
    if not clean_missing(row["classifier_groups"]):
        raise ValueError("VirSorter2 metadata contains no classifier groups")
    if row["run_status"] not in {
        "completed_with_virus_calls",
        "completed_no_viruses_detected",
    }:
        raise ValueError(f"Unrecognized VirSorter2 run status: {row['run_status']}")
    parse_nonnegative_int(row["virus_call_count"], "virus call count", sample_id)
    return row


def unique_rows(
    rows: list[dict[str, str]], key: str, label: str
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        value = row[key].strip()
        if not value:
            raise ValueError(f"{label} contains an empty {key}")
        if value in indexed:
            raise ValueError(f"Duplicate {key} in {label}: {value}")
        indexed[value] = row
    return indexed


def parse_score(value: str, label: str, sequence_id: str) -> tuple[float, str]:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing {label} for {sequence_id}")
    numeric = float(cleaned)
    if not 0.0 <= numeric <= 1.0:
        raise ValueError(f"{label} outside 0-1 for {sequence_id}: {cleaned}")
    return numeric, cleaned


def parse_percentage(value: str, label: str, sequence_id: str) -> str:
    cleaned = clean_missing(value)
    if not cleaned:
        return ""
    numeric = float(cleaned)
    if not 0.0 <= numeric <= 100.0:
        raise ValueError(f"{label} outside 0-100 for {sequence_id}: {cleaned}")
    return cleaned


def classify_evidence_strength(score: float) -> tuple[str, str] | None:
    """Apply the score cutoffs recommended by the VirSorter2 authors."""
    if score >= 0.90:
        return "strong", "virsorter2_high_confidence_cutoff"
    if score >= 0.50:
        return "qualified", "virsorter2_default_cutoff"
    return None


def parse_positive_int(value: str, label: str, sequence_id: str) -> int:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing {label} for {sequence_id}")
    numeric = int(cleaned)
    if numeric < 1:
        raise ValueError(f"{label} must be positive for {sequence_id}: {cleaned}")
    return numeric


def parse_nonnegative_int(value: str, label: str, sequence_id: str) -> int:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing {label} for {sequence_id}")
    numeric = int(cleaned)
    if numeric < 0:
        raise ValueError(f"{label} must be nonnegative for {sequence_id}: {cleaned}")
    return numeric


def parse_prediction_id(sequence_id: str) -> tuple[str, str, str]:
    if "||" not in sequence_id:
        raise ValueError(f"VirSorter2 prediction lacks a call suffix: {sequence_id}")
    parent_id, suffix = sequence_id.split("||", 1)
    if suffix == "full":
        return parent_id, "full", "input_contig"
    if suffix == "lt2gene":
        return parent_id, "lt2gene", "short_hallmark_region"
    if re.fullmatch(r"\d+_partial", suffix):
        return parent_id, "partial", "provirus"
    raise ValueError(f"Unrecognized VirSorter2 call suffix in {sequence_id}: {suffix}")


def boundary_fields(
    sequence_id: str,
    parent_id: str,
    call_type: str,
    region_length: int,
    score: float,
    max_group: str,
    hallmark_count: int,
    boundary: dict[str, str] | None,
    parent_length: int,
) -> dict[str, str]:
    if boundary is None:
        if call_type != "lt2gene":
            raise ValueError(f"Missing VirSorter2 boundary row for {sequence_id}")
        return {"coordinates": "", "topology": "", "n_genes": ""}

    if boundary["seqname"].strip() != parent_id:
        raise ValueError(f"Boundary parent mismatch for {sequence_id}")

    start = parse_positive_int(boundary["trim_bp_start"], "boundary start", sequence_id)
    end = parse_positive_int(boundary["trim_bp_end"], "boundary end", sequence_id)
    if end < start or end > parent_length:
        raise ValueError(f"Invalid boundary coordinates for {sequence_id}: {start}-{end}")
    if end - start + 1 != region_length:
        raise ValueError(f"Score and boundary lengths disagree for {sequence_id}")

    orf_start = parse_positive_int(
        boundary["trim_orf_index_start"], "ORF start", sequence_id
    )
    orf_end = parse_positive_int(boundary["trim_orf_index_end"], "ORF end", sequence_id)
    if orf_end < orf_start:
        raise ValueError(f"Invalid ORF range for {sequence_id}: {orf_start}-{orf_end}")

    boundary_score, _ = parse_score(
        boundary["final_max_score"], "boundary final score", sequence_id
    )
    if abs(boundary_score - score) > 1e-9:
        raise ValueError(f"Score and boundary maximum scores disagree for {sequence_id}")
    if boundary["final_max_score_group"].strip() != max_group:
        raise ValueError(f"Score and boundary maximum groups disagree for {sequence_id}")

    boundary_hallmarks = parse_nonnegative_int(
        boundary["hallmark_cnt"], "boundary hallmark count", sequence_id
    )
    if boundary_hallmarks != hallmark_count:
        raise ValueError(f"Score and boundary hallmark counts disagree for {sequence_id}")

    partial_flag = boundary["partial"].strip()
    expected_partial = "1" if call_type == "partial" else "0"
    if partial_flag != expected_partial:
        raise ValueError(f"Call suffix and boundary partial flag disagree for {sequence_id}")

    return {
        "coordinates": f"{start}-{end}",
        "topology": clean_missing(boundary["shape"]),
        "n_genes": str(orf_end - orf_start + 1),
    }


def main() -> None:
    args = parse_args()
    header_records = load_header_map(args.header_map, args.sample_id)
    metadata = load_run_metadata(
        args.run_metadata, args.sample_id, args.input_type
    )
    groups = [group.strip() for group in metadata["classifier_groups"].split(",")]
    min_score, _ = parse_score(metadata["min_score"], "minimum score", args.sample_id)

    score_columns, score_rows = read_tsv(args.score_table, SCORE_REQUIRED)
    missing_group_columns = set(groups).difference(score_columns)
    if missing_group_columns:
        raise ValueError(
            "VirSorter2 score table lacks configured classifier columns: "
            f"{sorted(missing_group_columns)}"
        )

    _, boundary_rows = read_tsv(args.boundary_table, BOUNDARY_REQUIRED)
    score_records = unique_rows(score_rows, "seqname", "VirSorter2 score table")
    boundary_records = unique_rows(
        boundary_rows, "seqname_new", "VirSorter2 boundary table"
    )

    metadata_call_count = parse_nonnegative_int(
        metadata["virus_call_count"], "virus call count", args.sample_id
    )
    if metadata_call_count != len(score_records):
        raise ValueError("VirSorter2 metadata and score-table call counts disagree")
    expected_status = (
        "completed_with_virus_calls"
        if score_records
        else "completed_no_viruses_detected"
    )
    if metadata["run_status"] != expected_status:
        raise ValueError("VirSorter2 metadata status disagrees with the score table")

    evidence_rows: list[dict[str, str]] = []
    for sequence_id, row in score_records.items():
        parent_id, call_type, record_type = parse_prediction_id(sequence_id)
        if parent_id not in header_records:
            raise ValueError(f"VirSorter2 parent is absent from header map: {parent_id}")

        parent_length = parse_positive_int(
            header_records[parent_id]["length"], "parent length", sequence_id
        )
        region_length = parse_positive_int(row["length"], "region length", sequence_id)
        if region_length > parent_length:
            raise ValueError(f"VirSorter2 region exceeds its parent for {sequence_id}")

        # VirSorter2 keeps ``lt2gene`` entries for short regions that could not
        # be assigned a classifier score. They are native diagnostic rows, not
        # viral calls, so they should not enter viSUM evidence.
        if call_type == "lt2gene" and not clean_missing(row["max_score"]):
            continue

        score, score_text = parse_score(row["max_score"], "maximum score", sequence_id)
        if score + 1e-9 < min_score:
            raise ValueError(f"VirSorter2 call is below the configured cutoff: {sequence_id}")
        strength = classify_evidence_strength(score)
        if strength is None:
            # A user may run VirSorter2 below its documented default cutoff for
            # exploratory purposes. Such rows remain in the native score table
            # but are not formal evidence for the viSUM discovery gate.
            continue
        evidence_strength, strength_basis = strength
        max_group = row["max_score_group"].strip()
        if max_group not in groups:
            raise ValueError(
                f"Maximum group for {sequence_id} was not configured: {max_group}"
            )

        hallmark_count = parse_nonnegative_int(
            row["hallmark"], "hallmark count", sequence_id
        )
        region = boundary_fields(
            sequence_id,
            parent_id,
            call_type,
            region_length,
            score,
            max_group,
            hallmark_count,
            boundary_records.get(sequence_id),
            parent_length,
        )

        output = {
            "sample_id": args.sample_id,
            "sequence_id": sequence_id,
            "parent_sequence_id": parent_id,
            "record_type": record_type,
            "coordinates": region["coordinates"],
            "tool": "virsorter2",
            "classification": "virus",
            "score": score_text,
            "score_type": "virsorter2_max_score",
            "length": str(region_length),
            "topology": region["topology"],
            "n_genes": region["n_genes"],
            "n_hallmarks": str(hallmark_count),
            "call_type": call_type,
            "parent_length": str(parent_length),
            "viral_gene_percent": parse_percentage(
                row["viral"], "viral gene percentage", sequence_id
            ),
            "cellular_gene_percent": parse_percentage(
                row["cellular"], "cellular gene percentage", sequence_id
            ),
            "max_score_group": max_group,
            "evidence_strength": evidence_strength,
            "strength_basis": strength_basis,
        }
        output.update(unclassified_taxonomy("virus"))
        evidence_rows.append(output)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(evidence_rows)

    print(f"Wrote {args.output}: {len(evidence_rows)} virus calls")


if __name__ == "__main__":
    main()
