#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from evidence_schema import (
    CORE_EVIDENCE_COLUMNS,
    RANK_PREFIXES,
    TAXONOMY_COLUMNS,
    unclassified_taxonomy,
)

OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "fdr",
    "genetic_code",
    "marker_enrichment",
    "n_uscg",
    "evidence_strength",
    "strength_basis",
]

VIRUS_REQUIRED = {
    "seq_name",
    "length",
    "topology",
    "coordinates",
    "n_genes",
    "genetic_code",
    "virus_score",
    "fdr",
    "n_hallmarks",
    "marker_enrichment",
    "taxonomy",
}

PLASMID_REQUIRED = {
    "seq_name",
    "length",
    "topology",
    "n_genes",
    "genetic_code",
    "plasmid_score",
    "fdr",
    "n_hallmarks",
    "marker_enrichment",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert geNomad summaries into standardized viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--virus-summary", required=True, type=Path)
    parser.add_argument("--virus-genes", required=True, type=Path)
    parser.add_argument("--plasmid-summary", required=True, type=Path)
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


def load_header_map(path: Path, sample_id: str) -> dict[str, dict[str, str]]:
    rows = read_tsv(path, {"sample_id", "sequence_id", "record_type"})
    records: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["sample_id"] != sample_id:
            raise ValueError(
                f"Header-map sample '{row['sample_id']}' does not match '{sample_id}'"
            )
        sequence_id = row["sequence_id"]
        if not sequence_id:
            raise ValueError("Header map contains an empty sequence_id")
        if sequence_id in records:
            raise ValueError(f"Duplicate sequence_id in header map: {sequence_id}")
        records[sequence_id] = row
    if not records:
        raise ValueError("Header map contains no sequence records")
    return records


def load_run_metadata(
    path: Path, sample_id: str, input_type: str
) -> dict[str, str]:
    rows = read_tsv(
        path,
        {
            "sample_id",
            "input_type",
            "genomad_version",
            "input_sequence_count",
            "score_calibration_requested",
            "score_calibration_applied",
            "run_status",
            "virus_call_count",
            "plasmid_call_count",
        },
    )
    if len(rows) != 1:
        raise ValueError(f"Expected one geNomad metadata row in {path}; found {len(rows)}")
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError("geNomad metadata does not match the requested sample and input type")
    if not clean_missing(row["genomad_version"]):
        raise ValueError("geNomad metadata contains an empty version")
    if row["run_status"] not in {
        "completed_with_virus_calls",
        "completed_no_viruses_detected",
    }:
        raise ValueError(f"Unrecognized geNomad run status: {row['run_status']}")
    parse_nonnegative_int(row["virus_call_count"], "virus call count", sample_id)
    parse_nonnegative_int(row["plasmid_call_count"], "plasmid call count", sample_id)
    parse_nonnegative_int(row["input_sequence_count"], "input sequence count", sample_id)
    parse_bool(row["score_calibration_requested"], "score_calibration_requested")
    calibration_applied = parse_bool(
        row["score_calibration_applied"], "score_calibration_applied"
    )
    expected_applied = (
        parse_bool(
            row["score_calibration_requested"], "score_calibration_requested"
        )
        and int(row["input_sequence_count"]) >= 1000
    )
    if calibration_applied != expected_applied:
        raise ValueError(
            "geNomad calibration metadata disagrees with its >=1,000-sequence rule"
        )
    return row


def clean_missing(value: str | None) -> str:
    value = "" if value is None else value.strip()
    return "" if value.lower() in {"", "na", "nan", "none"} else value


def parse_bool(value: str, label: str) -> bool:
    normalized = clean_missing(value).lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"{label} must be true or false; found: {value}")
    return normalized == "true"


def validate_score(value: str, label: str, sequence_id: str) -> str:
    value = clean_missing(value)
    if not value:
        raise ValueError(f"Missing {label} for {sequence_id}")
    score = float(value)
    if not 0.0 <= score <= 1.0:
        raise ValueError(f"{label} outside 0-1 for {sequence_id}: {value}")
    return value


def parse_nonnegative_int(value: str, label: str, sequence_id: str) -> int:
    value = clean_missing(value)
    if not value:
        raise ValueError(f"Missing {label} for {sequence_id}")
    numeric = int(value)
    if numeric < 0:
        raise ValueError(f"{label} must be nonnegative for {sequence_id}: {value}")
    return numeric


def parse_float(value: str, label: str, sequence_id: str) -> float:
    value = clean_missing(value)
    if not value:
        raise ValueError(f"Missing {label} for {sequence_id}")
    return float(value)


def load_virus_gene_metrics(
    path: Path, virus_rows: list[dict[str, str]]
) -> dict[str, int]:
    gene_rows = read_tsv(path, {"gene", "uscg"})
    expected_gene_counts = {
        row["seq_name"].strip(): parse_nonnegative_int(
            row["n_genes"], "gene count", row["seq_name"].strip()
        )
        for row in virus_rows
    }
    observed_gene_counts = {sequence_id: 0 for sequence_id in expected_gene_counts}
    uscg_counts = {sequence_id: 0 for sequence_id in expected_gene_counts}

    for row in gene_rows:
        gene_id = row["gene"].strip()
        sequence_id, separator, gene_number = gene_id.rpartition("_")
        if not separator or not gene_number.isdigit():
            raise ValueError(f"Unexpected geNomad virus gene identifier: {gene_id}")
        if sequence_id not in expected_gene_counts:
            raise ValueError(
                f"Virus gene does not map to a virus-summary sequence: {gene_id}"
            )
        uscg = clean_missing(row["uscg"])
        if uscg not in {"0", "1"}:
            raise ValueError(f"Invalid USCG indicator for {gene_id}: {row['uscg']}")
        observed_gene_counts[sequence_id] += 1
        uscg_counts[sequence_id] += int(uscg)

    for sequence_id, expected_count in expected_gene_counts.items():
        observed_count = observed_gene_counts[sequence_id]
        if observed_count != expected_count:
            raise ValueError(
                "geNomad virus summary and virus-gene table disagree for "
                f"{sequence_id}: expected {expected_count}, found {observed_count}"
            )
    return uscg_counts


def resolve_identity(
    sequence_id: str,
    coordinates: str,
    header_records: dict[str, dict[str, str]],
) -> tuple[str, str, str]:
    if sequence_id in header_records:
        return "", header_records[sequence_id]["record_type"], clean_missing(coordinates)

    provirus_token = "|provirus_"
    if provirus_token in sequence_id:
        parent_sequence_id = sequence_id.split(provirus_token, 1)[0]
        if parent_sequence_id not in header_records:
            raise ValueError(
                f"Provirus parent is absent from the header map: {parent_sequence_id}"
            )
        clean_coordinates = clean_missing(coordinates)
        if not clean_coordinates:
            clean_coordinates = sequence_id.split(provirus_token, 1)[1].replace("_", "-")
        return parent_sequence_id, "provirus", clean_coordinates

    raise ValueError(f"geNomad sequence is absent from the header map: {sequence_id}")


def taxonomy_fields(raw_taxonomy: str, classification: str) -> dict[str, str]:
    output = unclassified_taxonomy(classification)
    raw_taxonomy = clean_missing(raw_taxonomy)

    if classification == "virus":
        if raw_taxonomy and raw_taxonomy.lower() != "unclassified":
            tokens = [token.strip() for token in raw_taxonomy.split(";")]
            if tokens[0].lower() != "viruses":
                raise ValueError(f"Unexpected geNomad virus taxonomy: {raw_taxonomy}")
            for index, token in enumerate(tokens[:7]):
                rank = clean_missing(token) or "unclassified"
                output[TAXONOMY_COLUMNS[index]] = f"{RANK_PREFIXES[index]}{rank}"

    return output


def standardize_row(
    row: dict[str, str],
    classification: str,
    sample_id: str,
    header_records: dict[str, dict[str, str]],
    calibration_applied: bool,
    n_uscg: int | None = None,
) -> dict[str, str]:
    sequence_id = row["seq_name"].strip()
    if not sequence_id:
        raise ValueError("geNomad summary contains an empty seq_name")

    coordinates = row.get("coordinates", "")
    parent_id, record_type, coordinates = resolve_identity(
        sequence_id, coordinates, header_records
    )
    source_score_type = "virus_score" if classification == "virus" else "plasmid_score"
    score = validate_score(row[source_score_type], source_score_type, sequence_id)
    score_type = (
        f"calibrated_{classification}_probability"
        if calibration_applied
        else source_score_type
    )
    taxonomy = clean_missing(row.get("taxonomy", ""))
    n_genes = parse_nonnegative_int(row["n_genes"], "gene count", sequence_id)
    n_hallmarks = parse_nonnegative_int(
        row["n_hallmarks"], "hallmark count", sequence_id
    )
    marker_enrichment = parse_float(
        row["marker_enrichment"], "marker enrichment", sequence_id
    )
    fdr = clean_missing(row.get("fdr", ""))
    parsed_fdr: float | None = None
    if fdr:
        parsed_fdr = float(validate_score(fdr, "FDR", sequence_id))
    if calibration_applied and parsed_fdr is None:
        raise ValueError(f"Calibrated geNomad call lacks an FDR for {sequence_id}")
    if not calibration_applied and parsed_fdr is not None:
        raise ValueError(f"Uncalibrated geNomad call unexpectedly has an FDR: {sequence_id}")

    evidence_strength = "qualified"
    strength_basis = "genomad_default_summary"
    if classification == "virus":
        if n_uscg is None:
            raise ValueError(f"Missing USCG count for geNomad virus: {sequence_id}")
        passes_conservative = (
            float(score) >= 0.80
            and (parsed_fdr is None or parsed_fdr <= 0.05)
            and n_genes >= 1
            and marker_enrichment >= 1.50
            and n_hallmarks >= 1
            and n_uscg <= 2
        )
        if passes_conservative:
            evidence_strength = "strong"
            strength_basis = "genomad_conservative_preset"

    output = {
        "sample_id": sample_id,
        "sequence_id": sequence_id,
        "parent_sequence_id": parent_id,
        "record_type": record_type,
        "coordinates": coordinates,
        "tool": "genomad",
        "classification": classification,
        "score": score,
        "score_type": score_type,
        "length": row["length"].strip(),
        "topology": clean_missing(row.get("topology", "")),
        "n_genes": str(n_genes),
        "n_hallmarks": str(n_hallmarks),
        "fdr": fdr,
        "genetic_code": row["genetic_code"].strip(),
        "marker_enrichment": row["marker_enrichment"].strip(),
        "n_uscg": "" if n_uscg is None else str(n_uscg),
        "evidence_strength": evidence_strength,
        "strength_basis": strength_basis,
    }
    output.update(taxonomy_fields(taxonomy, classification))
    return output


def main() -> None:
    args = parse_args()
    header_records = load_header_map(args.header_map, args.sample_id)
    metadata = load_run_metadata(args.run_metadata, args.sample_id, args.input_type)

    evidence_rows: list[dict[str, str]] = []
    virus_rows = read_tsv(args.virus_summary, VIRUS_REQUIRED)
    plasmid_rows = read_tsv(args.plasmid_summary, PLASMID_REQUIRED)
    virus_uscg_counts = load_virus_gene_metrics(args.virus_genes, virus_rows)
    calibration_applied = parse_bool(
        metadata["score_calibration_applied"], "score_calibration_applied"
    )

    virus_call_count = parse_nonnegative_int(
        metadata["virus_call_count"], "virus call count", args.sample_id
    )
    plasmid_call_count = parse_nonnegative_int(
        metadata["plasmid_call_count"], "plasmid call count", args.sample_id
    )
    if virus_call_count != len(virus_rows):
        raise ValueError("geNomad metadata and virus-summary call counts disagree")
    if plasmid_call_count != len(plasmid_rows):
        raise ValueError("geNomad metadata and plasmid-summary call counts disagree")
    expected_status = (
        "completed_with_virus_calls"
        if virus_rows
        else "completed_no_viruses_detected"
    )
    if metadata["run_status"] != expected_status:
        raise ValueError("geNomad metadata status disagrees with the virus summary")

    for row in virus_rows:
        evidence_rows.append(
            standardize_row(
                row,
                "virus",
                args.sample_id,
                header_records,
                calibration_applied,
                virus_uscg_counts[row["seq_name"].strip()],
            )
        )
    for row in plasmid_rows:
        evidence_rows.append(
            standardize_row(
                row,
                "plasmid",
                args.sample_id,
                header_records,
                calibration_applied,
            )
        )

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
    plasmid_count = sum(row["classification"] == "plasmid" for row in evidence_rows)
    print(
        f"Wrote {args.output}: {virus_count} virus calls, "
        f"{plasmid_count} plasmid calls"
    )


if __name__ == "__main__":
    main()
