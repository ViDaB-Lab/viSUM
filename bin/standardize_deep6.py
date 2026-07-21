#!/usr/bin/env python3

"""Convert Deep6 predictions into sparse, standardized viSUM evidence."""

import argparse
import csv
import math
import statistics
from pathlib import Path

from evidence_schema import CORE_EVIDENCE_COLUMNS, unclassified_taxonomy


DEEP6_CLASSES = ("duplo", "euk", "mono", "pro", "ribo", "vari")
SCORE_REQUIRED = {"name", "length", *DEEP6_CLASSES}
OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + ["deep6_class"]

CLASSIFICATION_MAP = {
    "duplo": "virus",
    "euk": "cellular",
    "mono": "virus",
    "pro": "cellular",
    "ribo": "virus",
    "vari": "virus",
}

REALM_MAP = {
    "duplo": "Duplodnaviria",
    "mono": "Monodnaviria",
    "ribo": "Riboviria",
    "vari": "Varidnaviria",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Deep6 scores into standardized viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("rna",))
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--score-table", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--minimum-score", required=True, type=float)
    parser.add_argument("--median-multiplier", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path, required_columns: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        missing = required_columns.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return list(reader)


def clean_missing(value: str | None) -> str:
    value = "" if value is None else value.strip()
    return "" if value.lower() in {"", "na", "nan", "none"} else value


def parse_integer(
    value: str | None,
    label: str,
    sequence_id: str,
    minimum: int = 0,
) -> int:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing {label} for {sequence_id}")
    numeric = float(cleaned)
    if not math.isfinite(numeric) or not numeric.is_integer():
        raise ValueError(f"{label} is not an integer for {sequence_id}: {cleaned}")
    integer = int(numeric)
    if integer < minimum:
        raise ValueError(
            f"{label} must be at least {minimum} for {sequence_id}: {cleaned}"
        )
    return integer


def parse_score(value: str | None, label: str, sequence_id: str) -> tuple[float, str]:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing {label} for {sequence_id}")
    numeric = float(cleaned)
    if not math.isfinite(numeric) or not 0.0 <= numeric <= 1.0:
        raise ValueError(f"{label} outside 0-1 for {sequence_id}: {cleaned}")
    return numeric, cleaned


def validate_thresholds(minimum_score: float, median_multiplier: float) -> None:
    if not math.isfinite(minimum_score) or not 0.0 <= minimum_score <= 1.0:
        raise ValueError("Deep6 minimum score must be between 0 and 1")
    if not math.isfinite(median_multiplier) or median_multiplier < 1.0:
        raise ValueError("Deep6 median multiplier must be at least 1")


def load_header_map(path: Path, sample_id: str) -> dict[str, dict[str, str]]:
    rows = read_tsv(path, {"sample_id", "sequence_id", "record_type", "length"})
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
        parse_integer(row["length"], "header-map length", sequence_id, minimum=1)
        records[sequence_id] = row
    if not records:
        raise ValueError("Header map contains no sequence records")
    return records


def load_run_metadata(
    path: Path,
    sample_id: str,
    input_type: str,
) -> dict[str, str]:
    rows = read_tsv(
        path,
        {
            "sample_id",
            "input_type",
            "minimum_length",
            "prediction_count",
            "deep6_version",
            "deep6_revision",
            "raw_score_file",
        },
    )
    if len(rows) != 1:
        raise ValueError(f"Expected one Deep6 metadata row in {path}; found {len(rows)}")
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError("Deep6 metadata does not match the requested sample and input type")
    parse_integer(row["minimum_length"], "minimum length", sample_id, minimum=1)
    parse_integer(row["prediction_count"], "prediction count", sample_id)
    if not clean_missing(row["deep6_version"]):
        raise ValueError("Deep6 metadata contains an empty version")
    if not clean_missing(row["deep6_revision"]):
        raise ValueError("Deep6 metadata contains an empty revision")
    if not clean_missing(row["raw_score_file"]):
        raise ValueError("Deep6 metadata contains an empty raw score filename")
    return row


def taxonomy_fields(deep6_class: str, classification: str) -> dict[str, str]:
    taxonomy = unclassified_taxonomy(classification)
    if deep6_class in REALM_MAP:
        taxonomy["r__Realm"] = f"r__{REALM_MAP[deep6_class]}"
    elif deep6_class == "euk":
        taxonomy["d__Domain"] = "d__Eukaryota"
    elif deep6_class == "pro":
        taxonomy["d__Domain"] = "d__Prokaryota"
    return taxonomy


def main() -> None:
    args = parse_args()
    validate_thresholds(args.minimum_score, args.median_multiplier)
    header_records = load_header_map(args.header_map, args.sample_id)
    metadata = load_run_metadata(args.run_metadata, args.sample_id, args.input_type)
    score_rows = read_tsv(args.score_table, SCORE_REQUIRED)

    metadata_prediction_count = parse_integer(
        metadata["prediction_count"], "prediction count", args.sample_id
    )
    if metadata_prediction_count != len(score_rows):
        raise ValueError("Deep6 metadata and score-table prediction counts disagree")

    minimum_length = parse_integer(
        metadata["minimum_length"], "minimum length", args.sample_id, minimum=1
    )
    expected_ids = {
        sequence_id
        for sequence_id, row in header_records.items()
        if parse_integer(
            row["length"], "header-map length", sequence_id, minimum=1
        )
        >= minimum_length
    }

    score_records: dict[str, dict[str, str]] = {}
    for row in score_rows:
        sequence_id = row["name"].strip()
        if not sequence_id:
            raise ValueError("Deep6 score table contains an empty name")
        if sequence_id in score_records:
            raise ValueError(f"Duplicate Deep6 sequence name: {sequence_id}")
        if sequence_id not in header_records:
            raise ValueError(f"Deep6 sequence is absent from the header map: {sequence_id}")
        score_records[sequence_id] = row

    if set(score_records) != expected_ids:
        missing = sorted(expected_ids.difference(score_records))
        unexpected = sorted(set(score_records).difference(expected_ids))
        details = []
        if missing:
            details.append(f"missing eligible sequences: {missing[:5]}")
        if unexpected:
            details.append(f"unexpected ineligible sequences: {unexpected[:5]}")
        raise ValueError("Deep6 score-table coverage mismatch; " + "; ".join(details))

    evidence_rows: list[dict[str, str]] = []
    for sequence_id, row in score_records.items():
        reported_length = parse_integer(
            row["length"], "Deep6 sequence length", sequence_id, minimum=1
        )
        header_length = parse_integer(
            header_records[sequence_id]["length"],
            "header-map length",
            sequence_id,
            minimum=1,
        )
        if reported_length != header_length:
            raise ValueError(
                f"Deep6 and header-map lengths disagree for {sequence_id}"
            )

        numeric_scores: list[float] = []
        score_text: dict[str, str] = {}
        for deep6_class in DEEP6_CLASSES:
            numeric, cleaned = parse_score(
                row[deep6_class], f"{deep6_class} score", sequence_id
            )
            numeric_scores.append(numeric)
            score_text[deep6_class] = cleaned
        if abs(sum(numeric_scores) - 1.0) > 1e-4:
            raise ValueError(f"Deep6 class scores do not sum to 1 for {sequence_id}")

        top_score = max(numeric_scores)
        top_indexes = [
            index
            for index, value in enumerate(numeric_scores)
            if math.isclose(value, top_score, rel_tol=0.0, abs_tol=1e-12)
        ]
        median_score = statistics.median(numeric_scores)
        if (
            len(top_indexes) != 1
            or top_score < args.minimum_score
            or top_score < args.median_multiplier * median_score
        ):
            continue

        deep6_class = DEEP6_CLASSES[top_indexes[0]]
        classification = CLASSIFICATION_MAP[deep6_class]
        output = {
            "sample_id": args.sample_id,
            "sequence_id": sequence_id,
            "parent_sequence_id": "",
            "record_type": header_records[sequence_id]["record_type"].strip(),
            "coordinates": "",
            "tool": "deep6",
            "classification": classification,
            "score": score_text[deep6_class],
            "score_type": "deep6_top_class_score",
            "length": str(reported_length),
            "topology": "",
            "n_genes": "",
            "n_hallmarks": "",
            "deep6_class": deep6_class,
        }
        output.update(taxonomy_fields(deep6_class, classification))
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

    viral_count = sum(row["classification"] == "virus" for row in evidence_rows)
    cellular_count = sum(
        row["classification"] == "cellular" for row in evidence_rows
    )
    print(
        f"Wrote {args.output}: {len(score_rows)} predictions, "
        f"{len(evidence_rows)} confident calls "
        f"({viral_count} viral, {cellular_count} cellular)"
    )


if __name__ == "__main__":
    main()
