#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


OUTPUT_COLUMNS = [
    "sample_id",
    "input_type",
    "sequence_id",
    "parent_sequence_id",
    "record_type",
    "coordinates",
    "tool",
    "tool_version",
    "classification",
    "score",
    "score_type",
    "score_calibrated",
    "fdr",
    "length",
    "topology",
    "n_genes",
    "genetic_code",
    "n_hallmarks",
    "marker_enrichment",
    "taxonomy",
    "classification_string",
    "d__Domain",
    "r__Realm",
    "k__Kingdom",
    "p__Phylum",
    "c__Class",
    "o__Order",
    "f__Family",
    "g__Genus",
    "s__Species",
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
    parser.add_argument("--plasmid-summary", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--score-calibrated", required=True, choices=("true", "false"))
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


def load_tool_version(path: Path, sample_id: str, input_type: str) -> str:
    rows = read_tsv(path, {"sample_id", "input_type", "genomad_version"})
    if len(rows) != 1:
        raise ValueError(f"Expected one geNomad metadata row in {path}; found {len(rows)}")
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError("geNomad metadata does not match the requested sample and input type")
    return row["genomad_version"].removeprefix("geNomad, version ").strip()


def clean_missing(value: str | None) -> str:
    value = "" if value is None else value.strip()
    return "" if value.lower() in {"", "na", "nan", "none"} else value


def validate_score(value: str, label: str, sequence_id: str) -> str:
    value = clean_missing(value)
    if not value:
        raise ValueError(f"Missing {label} for {sequence_id}")
    score = float(value)
    if not 0.0 <= score <= 1.0:
        raise ValueError(f"{label} outside 0-1 for {sequence_id}: {value}")
    return value


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
    ranks = ["unclassified"] * 9
    raw_taxonomy = clean_missing(raw_taxonomy)

    if classification == "virus":
        ranks[0] = "Viruses"
        if raw_taxonomy and raw_taxonomy.lower() != "unclassified":
            tokens = [token.strip() for token in raw_taxonomy.split(";")]
            if tokens[0].lower() != "viruses":
                raise ValueError(f"Unexpected geNomad virus taxonomy: {raw_taxonomy}")
            for index, token in enumerate(tokens[:7]):
                ranks[index] = clean_missing(token) or "unclassified"

    prefixes = ("d__", "r__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")
    values = [f"{prefix}{rank}" for prefix, rank in zip(prefixes, ranks)]
    return {
        "classification_string": ";".join(values),
        "d__Domain": values[0],
        "r__Realm": values[1],
        "k__Kingdom": values[2],
        "p__Phylum": values[3],
        "c__Class": values[4],
        "o__Order": values[5],
        "f__Family": values[6],
        "g__Genus": values[7],
        "s__Species": values[8],
    }


def standardize_row(
    row: dict[str, str],
    classification: str,
    sample_id: str,
    input_type: str,
    tool_version: str,
    score_calibrated: str,
    header_records: dict[str, dict[str, str]],
) -> dict[str, str]:
    sequence_id = row["seq_name"].strip()
    if not sequence_id:
        raise ValueError("geNomad summary contains an empty seq_name")

    coordinates = row.get("coordinates", "")
    parent_id, record_type, coordinates = resolve_identity(
        sequence_id, coordinates, header_records
    )
    score_type = "virus_score" if classification == "virus" else "plasmid_score"
    score = validate_score(row[score_type], score_type, sequence_id)
    taxonomy = clean_missing(row.get("taxonomy", ""))

    output = {
        "sample_id": sample_id,
        "input_type": input_type,
        "sequence_id": sequence_id,
        "parent_sequence_id": parent_id,
        "record_type": record_type,
        "coordinates": coordinates,
        "tool": "genomad",
        "tool_version": tool_version,
        "classification": classification,
        "score": score,
        "score_type": score_type,
        "score_calibrated": score_calibrated,
        "fdr": clean_missing(row.get("fdr", "")),
        "length": row["length"].strip(),
        "topology": clean_missing(row.get("topology", "")),
        "n_genes": row["n_genes"].strip(),
        "genetic_code": row["genetic_code"].strip(),
        "n_hallmarks": row["n_hallmarks"].strip(),
        "marker_enrichment": row["marker_enrichment"].strip(),
        "taxonomy": taxonomy,
    }
    output.update(taxonomy_fields(taxonomy, classification))
    return output


def main() -> None:
    args = parse_args()
    header_records = load_header_map(args.header_map, args.sample_id)
    tool_version = load_tool_version(args.run_metadata, args.sample_id, args.input_type)

    evidence_rows: list[dict[str, str]] = []
    for row in read_tsv(args.virus_summary, VIRUS_REQUIRED):
        evidence_rows.append(
            standardize_row(
                row,
                "virus",
                args.sample_id,
                args.input_type,
                tool_version,
                args.score_calibrated,
                header_records,
            )
        )
    for row in read_tsv(args.plasmid_summary, PLASMID_REQUIRED):
        evidence_rows.append(
            standardize_row(
                row,
                "plasmid",
                args.sample_id,
                args.input_type,
                tool_version,
                args.score_calibrated,
                header_records,
            )
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
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
