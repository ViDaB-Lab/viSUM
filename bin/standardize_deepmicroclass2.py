#!/usr/bin/env python3

"""Convert DeepMicroClass2 probabilities into sparse viSUM evidence."""

import argparse
import csv
import math
from pathlib import Path

from evidence_schema import CORE_EVIDENCE_COLUMNS, unclassified_taxonomy


DEEPMICROCLASS2_CLASSES = (
    "arc",
    "bac",
    "chlor",
    "euk",
    "eukvir",
    "mit",
    "pls",
    "prokvir",
)

# Defaults in DeepMicroClass2 predict.py at the viSUM-pinned revision.
CLASS_THRESHOLDS = {
    "arc": 0.625000,
    "bac": 0.400000,
    "chlor": 0.390625,
    "euk": 0.435547,
    "eukvir": 0.951172,
    "mit": 0.261719,
    "pls": 0.832031,
    "prokvir": 0.997925,
}
THRESHOLD_METADATA = ",".join(
    f"{label}={CLASS_THRESHOLDS[label]:.6f}"
    for label in DEEPMICROCLASS2_CLASSES
)

MODEL_MINIMUM_LENGTH = {
    "8class": 500,
    "high_precision": 300,
    "300bp": 300,
}

CLASSIFICATION_MAP = {
    "arc": "cellular",
    "bac": "cellular",
    "chlor": "cellular",
    "euk": "cellular",
    "eukvir": "virus",
    "mit": "cellular",
    "pls": "plasmid",
    "prokvir": "virus",
}

# Chloroplast and mitochondrial calls are cellular evidence, but assigning
# either organelle to an organismal domain would overstate the model output.
DOMAIN_MAP = {
    "arc": "Archaea",
    "bac": "Bacteria",
    "euk": "Eukaryota",
}

SCORE_REQUIRED = {
    "contig",
    "label",
    "confidence",
    *DEEPMICROCLASS2_CLASSES,
}
OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "deepmicroclass2_class",
    "evidence_strength",
    "strength_basis",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert DeepMicroClass2 scores into standardized viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna",))
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--score-table", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
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


def parse_probability(
    value: str | None,
    label: str,
    sequence_id: str,
) -> tuple[float, str]:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing {label} for {sequence_id}")
    numeric = float(cleaned)
    if not math.isfinite(numeric) or not 0.0 <= numeric <= 1.0:
        raise ValueError(f"{label} outside 0-1 for {sequence_id}: {cleaned}")
    return numeric, cleaned


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
            "model_mode",
            "minimum_length",
            "input_sequence_count",
            "eligible_sequence_count",
            "prediction_count",
            "class_thresholds",
            "deepmicroclass2_revision",
            "run_status",
            "score_file",
        },
    )
    if len(rows) != 1:
        raise ValueError(
            f"Expected one DeepMicroClass2 metadata row in {path}; found {len(rows)}"
        )
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError(
            "DeepMicroClass2 metadata does not match the requested sample and input type"
        )

    model_mode = clean_missing(row["model_mode"])
    if model_mode not in MODEL_MINIMUM_LENGTH:
        raise ValueError(f"Unrecognized DeepMicroClass2 model mode: {model_mode}")
    minimum_length = parse_integer(
        row["minimum_length"], "minimum length", sample_id, minimum=1
    )
    if minimum_length != MODEL_MINIMUM_LENGTH[model_mode]:
        raise ValueError(
            "DeepMicroClass2 metadata minimum length disagrees with its model mode"
        )
    if row["class_thresholds"] != THRESHOLD_METADATA:
        raise ValueError(
            "DeepMicroClass2 metadata thresholds disagree with the viSUM decision rule"
        )
    if row["run_status"] not in {"completed", "completed_no_eligible_sequences"}:
        raise ValueError(
            f"Unrecognized DeepMicroClass2 run status: {row['run_status']}"
        )
    if not clean_missing(row["deepmicroclass2_revision"]):
        raise ValueError("DeepMicroClass2 metadata contains an empty revision")
    if not clean_missing(row["score_file"]):
        raise ValueError("DeepMicroClass2 metadata contains an empty score filename")
    return row


def taxonomy_fields(
    deepmicroclass2_class: str,
    classification: str,
) -> dict[str, str]:
    taxonomy = unclassified_taxonomy(classification)
    if deepmicroclass2_class in DOMAIN_MAP:
        taxonomy["d__Domain"] = f"d__{DOMAIN_MAP[deepmicroclass2_class]}"
    return taxonomy


def main() -> None:
    args = parse_args()
    header_records = load_header_map(args.header_map, args.sample_id)
    metadata = load_run_metadata(
        args.run_metadata,
        args.sample_id,
        args.input_type,
    )
    score_rows = read_tsv(args.score_table, SCORE_REQUIRED)

    input_count = parse_integer(
        metadata["input_sequence_count"], "input sequence count", args.sample_id
    )
    eligible_count = parse_integer(
        metadata["eligible_sequence_count"],
        "eligible sequence count",
        args.sample_id,
    )
    prediction_count = parse_integer(
        metadata["prediction_count"], "prediction count", args.sample_id
    )
    if input_count != len(header_records):
        raise ValueError("DeepMicroClass2 metadata and header-map input counts disagree")
    if prediction_count != len(score_rows):
        raise ValueError(
            "DeepMicroClass2 metadata and score-table prediction counts disagree"
        )

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
    if eligible_count != len(expected_ids):
        raise ValueError(
            "DeepMicroClass2 metadata and header-map eligible counts disagree"
        )

    score_records: dict[str, dict[str, str]] = {}
    for row in score_rows:
        sequence_id = row["contig"].strip()
        if not sequence_id:
            raise ValueError("DeepMicroClass2 score table contains an empty contig")
        if sequence_id in score_records:
            raise ValueError(f"Duplicate DeepMicroClass2 contig: {sequence_id}")
        if sequence_id not in header_records:
            raise ValueError(
                f"DeepMicroClass2 contig is absent from the header map: {sequence_id}"
            )
        score_records[sequence_id] = row

    if set(score_records) != expected_ids:
        missing = sorted(expected_ids.difference(score_records))
        unexpected = sorted(set(score_records).difference(expected_ids))
        details = []
        if missing:
            details.append(f"missing eligible sequences: {missing[:5]}")
        if unexpected:
            details.append(f"unexpected ineligible sequences: {unexpected[:5]}")
        raise ValueError(
            "DeepMicroClass2 score-table coverage mismatch; " + "; ".join(details)
        )

    expected_status = (
        "completed" if score_rows else "completed_no_eligible_sequences"
    )
    if metadata["run_status"] != expected_status:
        raise ValueError(
            "DeepMicroClass2 metadata status disagrees with the prediction count"
        )

    evidence_rows: list[dict[str, str]] = []
    for sequence_id, row in score_records.items():
        probabilities: dict[str, float] = {}
        probability_text: dict[str, str] = {}
        for class_name in DEEPMICROCLASS2_CLASSES:
            numeric, cleaned = parse_probability(
                row[class_name], f"{class_name} probability", sequence_id
            )
            probabilities[class_name] = numeric
            probability_text[class_name] = cleaned
        if abs(sum(probabilities.values()) - 1.0) > 5e-4:
            raise ValueError(
                f"DeepMicroClass2 probabilities do not sum to 1 for {sequence_id}"
            )

        assigned_label = clean_missing(row["label"])
        if assigned_label not in DEEPMICROCLASS2_CLASSES:
            raise ValueError(
                f"Unrecognized DeepMicroClass2 assigned label for {sequence_id}: "
                f"{assigned_label}"
            )
        assigned_confidence, _ = parse_probability(
            row["confidence"], "assigned confidence", sequence_id
        )
        if not math.isclose(
            assigned_confidence,
            probabilities[assigned_label],
            rel_tol=0.0,
            abs_tol=5e-5,
        ):
            raise ValueError(
                f"DeepMicroClass2 label and confidence disagree for {sequence_id}"
            )

        top_score = max(probabilities.values())
        top_classes = [
            class_name
            for class_name, probability in probabilities.items()
            if math.isclose(probability, top_score, rel_tol=0.0, abs_tol=1e-12)
        ]
        if len(top_classes) != 1:
            continue

        top_class = top_classes[0]
        if top_score < CLASS_THRESHOLDS[top_class]:
            continue

        classification = CLASSIFICATION_MAP[top_class]
        header_row = header_records[sequence_id]
        output = {
            "sample_id": args.sample_id,
            "sequence_id": sequence_id,
            "parent_sequence_id": "",
            "record_type": header_row["record_type"].strip(),
            "coordinates": "",
            "tool": "deepmicroclass2",
            "classification": classification,
            "score": probability_text[top_class],
            "score_type": "deepmicroclass2_top_class_probability",
            "length": str(
                parse_integer(
                    header_row["length"],
                    "header-map length",
                    sequence_id,
                    minimum=1,
                )
            ),
            "topology": "",
            "n_genes": "",
            "n_hallmarks": "",
            "deepmicroclass2_class": top_class,
            "evidence_strength": "qualified",
            "strength_basis": "deepmicroclass2_top_class_official_threshold",
        }
        output.update(taxonomy_fields(top_class, classification))
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

    virus_count = sum(row["classification"] == "virus" for row in evidence_rows)
    cellular_count = sum(
        row["classification"] == "cellular" for row in evidence_rows
    )
    plasmid_count = sum(
        row["classification"] == "plasmid" for row in evidence_rows
    )
    print(
        f"Wrote {args.output}: {len(score_rows)} predictions, "
        f"{len(evidence_rows)} confident calls "
        f"({virus_count} viral, {cellular_count} cellular, "
        f"{plasmid_count} plasmid)"
    )


if __name__ == "__main__":
    main()
