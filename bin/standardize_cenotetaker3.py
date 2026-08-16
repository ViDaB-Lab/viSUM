#!/usr/bin/env python3

"""Convert Cenote-Taker 3 calls into sparse, standardized viSUM evidence."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from evidence_schema import (
    CORE_EVIDENCE_COLUMNS,
    RANK_PREFIXES,
    TAXONOMY_COLUMNS,
    unclassified_taxonomy,
)


OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "parent_length",
    "virion_hallmark_count",
    "dna_rep_hallmark_count",
    "rdrp_hallmark_count",
    "virion_hallmark_genes",
    "dna_rep_hallmark_genes",
    "rdrp_hallmark_genes",
    "genetic_code",
    "evidence_strength",
    "strength_basis",
]

SUMMARY_REQUIRED = {
    "contig",
    "input_name",
    "virus_seq_length",
    "end_feature",
    "gene_count",
    "virion_hallmark_count",
    "rep_hallmark_count",
    "RDRP_hallmark_count",
    "virion_hallmark_genes",
    "rep_hallmark_genes",
    "RDRP_hallmark_genes",
    "taxonomy_hierarchy",
    "gcode",
}

PRUNE_REQUIRED = {
    "contig",
    "contig_length",
    "chunk_length",
    "chunk_name",
    "chunk_start",
    "chunk_stop",
}

GENE_REQUIRED = {
    "contig",
    "gene_start",
    "gene_stop",
    "gene_name",
    "chunk_name",
}

RANK_CODE_TO_INDEX = {
    "d": 0,
    "r": 1,
    "k": 2,
    "p": 3,
    "c": 4,
    "o": 5,
    "f": 6,
    "g": 7,
    "s": 8,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Cenote-Taker 3 outputs into standardized viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--virus-summary", required=True, type=Path)
    parser.add_argument("--virus-fasta", required=True, type=Path)
    parser.add_argument("--prune-summary", required=True, type=Path)
    parser.add_argument("--gene-annotations", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def clean_missing(value: str | None) -> str:
    value = "" if value is None else value.strip()
    return "" if value.lower() in {"", "na", "nan", "none"} else value


def read_tsv(path: Path, required_columns: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        missing = required_columns.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return list(reader)


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
    if not numeric.is_integer():
        raise ValueError(f"{label} is not an integer for {sequence_id}: {cleaned}")
    integer = int(numeric)
    if integer < minimum:
        raise ValueError(
            f"{label} must be at least {minimum} for {sequence_id}: {cleaned}"
        )
    return integer


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


def load_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_length = 0

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    if current_length < 1:
                        raise ValueError(f"Empty FASTA record in {path}: {current_id}")
                    lengths[current_id] = current_length
                current_id = line[1:].split(maxsplit=1)[0]
                if not current_id:
                    raise ValueError(f"Empty FASTA identifier in {path} at line {line_number}")
                if current_id in lengths:
                    raise ValueError(f"Duplicate FASTA identifier in {path}: {current_id}")
                current_length = 0
            else:
                if current_id is None:
                    raise ValueError(f"Sequence data precedes a FASTA header in {path}")
                current_length += len(line)

    if current_id is not None:
        if current_length < 1:
            raise ValueError(f"Empty FASTA record in {path}: {current_id}")
        lengths[current_id] = current_length
    return lengths


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
            "cenotetaker3_version",
            "run_status",
            "virus_call_count",
        },
    )
    if len(rows) != 1:
        raise ValueError(f"Expected one CT3 metadata row in {path}; found {len(rows)}")
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError("CT3 metadata does not match the requested sample and input type")
    if not clean_missing(row["cenotetaker3_version"]):
        raise ValueError("CT3 metadata contains an empty version")
    if row["run_status"] not in {
        "completed_with_virus_calls",
        "completed_no_viruses_detected",
    }:
        raise ValueError(f"Unrecognized CT3 run status: {row['run_status']}")
    parse_integer(row["virus_call_count"], "virus call count", sample_id)
    return row


def unique_rows(
    rows: list[dict[str, str]],
    key: str,
    label: str,
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


def parse_taxonomy(value: str | None) -> dict[str, str]:
    taxonomy = unclassified_taxonomy("virus")
    ranks = [taxonomy[column].split("__", 1)[1] for column in TAXONOMY_COLUMNS]

    for token in clean_missing(value).split(";"):
        token = token.strip()
        if "_" not in token:
            continue
        code, label = token.split("_", 1)
        label = clean_missing(label.lstrip("_"))
        if not label or label.lower() == "unclassified virus":
            continue

        if code == "-":
            if label.lower() == "viruses":
                ranks[0] = "Viruses"
            elif label.lower().endswith("viria") and ranks[1] == "unclassified":
                ranks[1] = label
            continue

        rank_index = RANK_CODE_TO_INDEX.get(code.lower())
        if rank_index is not None:
            ranks[rank_index] = label

    ranks[0] = "Viruses"
    return {
        column: f"{prefix}{rank}"
        for column, prefix, rank in zip(TAXONOMY_COLUMNS, RANK_PREFIXES, ranks)
    }


def gene_call_key(row: dict[str, str]) -> str:
    contig = row["contig"].strip()
    if not contig:
        raise ValueError("CT3 gene annotations contain an empty contig")
    chunk_name = clean_missing(row["chunk_name"])
    return f"{contig}@{chunk_name}" if chunk_name else contig


def load_gene_counts(
    rows: list[dict[str, str]],
    fasta_lengths: dict[str, int],
) -> dict[str, int]:
    genes_by_call: dict[str, dict[str, tuple[int, int]]] = defaultdict(dict)
    for row in rows:
        call_key = gene_call_key(row)
        if call_key not in fasta_lengths:
            raise ValueError(f"Gene annotation has no matching CT3 virus sequence: {call_key}")
        gene_name = row["gene_name"].strip()
        if not gene_name:
            raise ValueError(f"CT3 gene annotation has an empty gene_name for {call_key}")
        start = parse_integer(row["gene_start"], "gene start", call_key, minimum=1)
        stop = parse_integer(row["gene_stop"], "gene stop", call_key, minimum=1)
        if stop < start or stop > fasta_lengths[call_key]:
            raise ValueError(f"Invalid gene coordinates for {call_key}: {start}-{stop}")
        coordinates = (start, stop)
        previous_coordinates = genes_by_call[call_key].get(gene_name)
        if previous_coordinates is not None and previous_coordinates != coordinates:
            raise ValueError(
                f"CT3 gene {gene_name} has conflicting coordinates for {call_key}"
            )
        genes_by_call[call_key][gene_name] = coordinates
    return {call_key: len(genes) for call_key, genes in genes_by_call.items()}


def load_pruned_regions(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    regions: dict[str, dict[str, str]] = {}
    for row in rows:
        chunk_name = clean_missing(row["chunk_name"])
        if not chunk_name:
            continue
        contig = row["contig"].strip()
        if not contig:
            raise ValueError("CT3 prune summary contains an empty contig")
        key = f"{contig}@{chunk_name}"
        if key in regions:
            raise ValueError(f"Duplicate CT3 pruned region: {key}")
        regions[key] = row
    return regions


def resolve_identity(
    row: dict[str, str],
    actual_length: int,
    header_records: dict[str, dict[str, str]],
    pruned_regions: dict[str, dict[str, str]],
) -> tuple[str, str, str, str, int]:
    call_id = row["contig"].strip()
    input_name = row["input_name"].strip()
    if input_name not in header_records:
        raise ValueError(f"CT3 input_name is absent from the header map: {input_name}")
    parent_length = parse_integer(
        header_records[input_name]["length"], "parent length", input_name, minimum=1
    )

    if "@" not in call_id:
        return (
            input_name,
            "",
            header_records[input_name]["record_type"],
            "",
            parent_length,
        )

    region = pruned_regions.get(call_id)
    if region is None:
        raise ValueError(f"Missing CT3 prune-summary region for {call_id}")
    region_parent_length = parse_integer(
        region["contig_length"], "prune parent length", call_id, minimum=1
    )
    if region_parent_length != parent_length:
        raise ValueError(f"Header-map and CT3 parent lengths disagree for {call_id}")

    raw_start = parse_integer(region["chunk_start"], "chunk start", call_id)
    raw_stop = parse_integer(region["chunk_stop"], "chunk stop", call_id, minimum=1)
    reported_chunk_length = parse_integer(
        region["chunk_length"], "chunk length", call_id, minimum=1
    )
    start = raw_start + 1
    end = raw_start + actual_length
    if end > parent_length:
        raise ValueError(f"CT3 pruned region exceeds its parent for {call_id}: {start}-{end}")

    # CT3 3.4.4 can report the stop and length one base too large for chunks
    # ending at the parent boundary. The emitted FASTA length is authoritative.
    if abs(raw_stop - end) > 1 or abs(reported_chunk_length - actual_length) > 1:
        raise ValueError(f"CT3 prune coordinates and FASTA length disagree for {call_id}")

    sequence_id = f"{input_name}|provirus_{start}_{end}"
    return sequence_id, input_name, "provirus", f"{start}-{end}", parent_length


def main() -> None:
    args = parse_args()
    header_records = load_header_map(args.header_map, args.sample_id)
    metadata = load_run_metadata(args.run_metadata, args.sample_id, args.input_type)
    summary_rows = read_tsv(args.virus_summary, SUMMARY_REQUIRED)
    prune_rows = read_tsv(args.prune_summary, PRUNE_REQUIRED)
    gene_rows = read_tsv(args.gene_annotations, GENE_REQUIRED)
    fasta_lengths = load_fasta_lengths(args.virus_fasta)

    metadata_call_count = parse_integer(
        metadata["virus_call_count"], "virus call count", args.sample_id
    )
    if metadata_call_count != len(summary_rows):
        raise ValueError("CT3 metadata and summary virus-call counts disagree")
    expected_status = (
        "completed_with_virus_calls"
        if summary_rows
        else "completed_no_viruses_detected"
    )
    if metadata["run_status"] != expected_status:
        raise ValueError("CT3 metadata status disagrees with the summary")

    summary_records = unique_rows(summary_rows, "contig", "CT3 virus summary")
    if set(summary_records) != set(fasta_lengths):
        raise ValueError("CT3 summary and virus FASTA identifiers disagree")

    gene_counts = load_gene_counts(gene_rows, fasta_lengths)
    if set(gene_counts).difference(summary_records):
        raise ValueError("CT3 gene annotations contain calls absent from the summary")
    pruned_regions = load_pruned_regions(prune_rows)

    evidence_rows: list[dict[str, str]] = []
    evidence_ids: set[str] = set()
    for call_id, row in summary_records.items():
        actual_length = fasta_lengths[call_id]
        reported_length = parse_integer(
            row["virus_seq_length"], "virus sequence length", call_id, minimum=1
        )
        allowed_length_error = 1 if "@" in call_id else 0
        if abs(reported_length - actual_length) > allowed_length_error:
            raise ValueError(f"CT3 summary and FASTA lengths disagree for {call_id}")

        reported_gene_count = parse_integer(row["gene_count"], "gene count", call_id)
        if gene_counts.get(call_id, 0) != reported_gene_count:
            raise ValueError(f"CT3 summary and annotation gene counts disagree for {call_id}")

        virion_count = parse_integer(
            row["virion_hallmark_count"], "virion hallmark count", call_id
        )
        dna_rep_count = parse_integer(
            row["rep_hallmark_count"], "DNA-replication hallmark count", call_id
        )
        rdrp_count = parse_integer(
            row["RDRP_hallmark_count"], "RdRP hallmark count", call_id
        )
        total_hallmarks = virion_count + dna_rep_count + rdrp_count
        if total_hallmarks < 1:
            # CT3 permits a user-selected zero-hallmark discovery threshold.
            # Preserve such calls in CT3's native outputs, but do not treat
            # them as qualified viral evidence in viSUM.
            continue

        sequence_id, parent_id, record_type, coordinates, parent_length = (
            resolve_identity(
                row,
                actual_length,
                header_records,
                pruned_regions,
            )
        )
        if sequence_id in evidence_ids:
            raise ValueError(f"Duplicate standardized CT3 sequence_id: {sequence_id}")
        evidence_ids.add(sequence_id)

        topology = clean_missing(row["end_feature"])
        is_linear = record_type == "provirus" or not topology
        if is_linear and total_hallmarks >= 2:
            evidence_strength = "strong"
            strength_basis = "cenotetaker3_two_hallmark_linear_recommendation"
        else:
            evidence_strength = "qualified"
            strength_basis = "cenotetaker3_default_hallmark_requirement"

        output = {
            "sample_id": args.sample_id,
            "sequence_id": sequence_id,
            "parent_sequence_id": parent_id,
            "record_type": record_type,
            "coordinates": coordinates,
            "tool": "cenotetaker3",
            "classification": "virus",
            "score": "",
            "score_type": "",
            "length": str(actual_length),
            "topology": topology,
            "n_genes": str(reported_gene_count),
            "n_hallmarks": str(total_hallmarks),
            "parent_length": str(parent_length),
            "virion_hallmark_count": str(virion_count),
            "dna_rep_hallmark_count": str(dna_rep_count),
            "rdrp_hallmark_count": str(rdrp_count),
            "virion_hallmark_genes": clean_missing(row["virion_hallmark_genes"]),
            "dna_rep_hallmark_genes": clean_missing(row["rep_hallmark_genes"]),
            "rdrp_hallmark_genes": clean_missing(row["RDRP_hallmark_genes"]),
            "genetic_code": clean_missing(row["gcode"]),
            "evidence_strength": evidence_strength,
            "strength_basis": strength_basis,
        }
        output.update(parse_taxonomy(row["taxonomy_hierarchy"]))
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
