#!/usr/bin/env python3

"""Convert CheckV quality and host-boundary calls into viSUM evidence."""

import argparse
import csv
from pathlib import Path
from typing import Iterator, TextIO

from evidence_schema import CORE_EVIDENCE_COLUMNS, unclassified_taxonomy


OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "evidence_strength",
    "strength_basis",
    "checkv_quality",
    "miuvig_quality",
    "completeness",
    "completeness_method",
    "contamination",
    "provirus",
    "proviral_length",
    "viral_genes",
    "host_genes",
    "complete_genome_prediction",
    "complete_genome_confidence",
    "aai_confidence",
    "aai_identity",
    "aai_alignment_fraction",
    "aai_hit_count",
    "hmm_completeness_lower",
    "hmm_completeness_upper",
    "hmm_hit_count",
    "warnings",
]

QUALITY_REQUIRED = {
    "contig_id",
    "contig_length",
    "provirus",
    "proviral_length",
    "gene_count",
    "viral_genes",
    "host_genes",
    "checkv_quality",
    "miuvig_quality",
    "completeness",
    "completeness_method",
    "contamination",
    "warnings",
}

COMPLETENESS_REQUIRED = {
    "contig_id",
    "contig_length",
    "viral_length",
    "aai_confidence",
    "aai_id",
    "aai_af",
    "aai_num_hits",
    "hmm_completeness_lower",
    "hmm_completeness_upper",
    "hmm_num_hits",
}

CONTAMINATION_REQUIRED = {
    "contig_id",
    "contig_length",
    "provirus",
    "proviral_length",
    "region_types",
    "region_coords_bp",
}

COMPLETE_GENOMES_REQUIRED = {
    "contig_id",
    "contig_length",
    "prediction_type",
    "confidence_level",
}

METADATA_REQUIRED = {
    "sample_id",
    "input_type",
    "candidate_sequence_count",
    "quality_summary_row_count",
    "determined_quality_count",
    "provirus_count",
    "checkv_version",
    "run_status",
}

DETERMINED_QUALITIES = {"complete", "high-quality", "medium-quality", "low-quality"}
MISSING_VALUES = {"", "na", "nan", "none"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert CheckV outputs into standardized viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--candidate-fasta", required=True, type=Path)
    parser.add_argument("--quality-summary", required=True, type=Path)
    parser.add_argument("--completeness", required=True, type=Path)
    parser.add_argument("--contamination", required=True, type=Path)
    parser.add_argument("--complete-genomes", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def clean_missing(value: str | None) -> str:
    cleaned = "" if value is None else value.strip()
    return "" if cleaned.lower() in MISSING_VALUES else cleaned


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return list(reader)


def unique_rows(
    rows: list[dict[str, str]], key: str, label: str
) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = row[key].strip()
        if not identifier:
            raise ValueError(f"{label} contains an empty {key}")
        if identifier in records:
            raise ValueError(f"Duplicate {label} row: {identifier}")
        records[identifier] = row
    return records


def parse_nonnegative_integer(value: str, label: str) -> int:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing {label}")
    try:
        parsed = int(cleaned)
    except ValueError as error:
        raise ValueError(f"{label} is not an integer: {value}") from error
    if parsed < 0:
        raise ValueError(f"{label} must be nonnegative: {value}")
    return parsed


def parse_positive_integer(value: str, label: str) -> int:
    parsed = parse_nonnegative_integer(value, label)
    if parsed < 1:
        raise ValueError(f"{label} must be positive: {value}")
    return parsed


def read_fasta(handle: TextIO) -> Iterator[tuple[str, str]]:
    identifier: str | None = None
    sequence_parts: list[str] = []
    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                yield identifier, "".join(sequence_parts)
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Empty FASTA header at line {line_number}")
            identifier = header.split(maxsplit=1)[0]
            sequence_parts = []
        else:
            if identifier is None:
                raise ValueError(
                    f"Sequence data precedes the first FASTA header at line {line_number}"
                )
            sequence_parts.append("".join(line.split()))
    if identifier is not None:
        yield identifier, "".join(sequence_parts)


def load_fasta_lengths(path: Path) -> dict[str, int]:
    records: dict[str, int] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for identifier, sequence in read_fasta(handle):
            if identifier in records:
                raise ValueError(f"Duplicate FASTA identifier: {identifier}")
            if not sequence:
                raise ValueError(f"Empty FASTA sequence: {identifier}")
            records[identifier] = len(sequence)
    return records


def parse_metadata(path: Path, sample_id: str, input_type: str) -> dict[str, str]:
    rows = read_tsv(path, METADATA_REQUIRED)
    if len(rows) != 1:
        raise ValueError(f"Expected one CheckV metadata row in {path}; found {len(rows)}")
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError("CheckV metadata does not match the requested sample and input type")
    if row["run_status"] not in {"completed", "skipped_no_discovery_candidates"}:
        raise ValueError(f"Unrecognized CheckV run status: {row['run_status']}")
    if not clean_missing(row["checkv_version"]):
        raise ValueError("CheckV metadata contains an empty version")
    return row


def parse_regions(row: dict[str, str], parent_length: int) -> list[tuple[int, int]]:
    region_types = clean_missing(row["region_types"])
    region_coordinates = clean_missing(row["region_coords_bp"])
    if not region_types or not region_coordinates:
        raise ValueError(f"CheckV provirus lacks region coordinates: {row['contig_id']}")

    types = [value.strip().lower() for value in region_types.split(",")]
    coordinates = [value.strip() for value in region_coordinates.split(",")]
    if len(types) != len(coordinates):
        raise ValueError(
            f"CheckV region types and coordinates disagree: {row['contig_id']}"
        )

    viral_regions: list[tuple[int, int]] = []
    for region_type, coordinate in zip(types, coordinates):
        if "-" not in coordinate:
            raise ValueError(
                f"Invalid CheckV region coordinate for {row['contig_id']}: {coordinate}"
            )
        start_text, end_text = coordinate.split("-", 1)
        start = parse_positive_integer(start_text, "CheckV region start")
        end = parse_positive_integer(end_text, "CheckV region end")
        if start > end or end > parent_length:
            raise ValueError(
                f"CheckV region exceeds {row['contig_id']}: {start}-{end}"
            )
        if region_type == "viral":
            viral_regions.append((start, end))

    if not viral_regions:
        raise ValueError(f"CheckV provirus contains no viral region: {row['contig_id']}")
    return viral_regions


def is_confident_complete_call(row: dict[str, str] | None) -> bool:
    if row is None:
        return False
    confidence = clean_missing(row["confidence_level"]).lower()
    return bool(confidence and confidence != "not-determined")


def build_evidence_row(
    sample_id: str,
    quality: dict[str, str],
    completeness_detail: dict[str, str],
    complete_call: dict[str, str] | None,
    sequence_id: str,
    parent_sequence_id: str,
    record_type: str,
    coordinates: str,
    length: int,
) -> dict[str, str]:
    completeness = clean_missing(quality["completeness"])
    checkv_quality = clean_missing(quality["checkv_quality"])
    aai_confidence = clean_missing(completeness_detail["aai_confidence"])
    confident_aai = aai_confidence.lower() in {"medium", "high"}
    high_completeness_tier = checkv_quality.lower() in {"complete", "high-quality"}
    viral_genes = parse_nonnegative_integer(
        clean_missing(quality["viral_genes"]) or "0", "CheckV viral gene count"
    )
    host_genes = parse_nonnegative_integer(
        clean_missing(quality["host_genes"]) or "0", "CheckV host gene count"
    )

    # CheckV quality describes estimated completeness, not probability of viral
    # origin.  Origin strength therefore requires gene-content support.
    if viral_genes == 0 and host_genes > 0:
        classification = "cellular"
        evidence_strength = "qualified"
        strength_basis = "checkv_host_genes_without_viral_genes"
    elif viral_genes == 0:
        classification = "unclassified"
        evidence_strength = "weak"
        strength_basis = "checkv_annotation_only_no_informative_genes"
    elif high_completeness_tier and confident_aai:
        classification = "virus"
        evidence_strength = "strong"
        strength_basis = "checkv_high_quality_confident_aai"
    elif record_type == "provirus":
        classification = "virus"
        evidence_strength = "qualified"
        strength_basis = "checkv_provirus_boundary_with_viral_genes"
    elif high_completeness_tier and host_genes == 0:
        classification = "virus"
        evidence_strength = "qualified"
        strength_basis = "checkv_high_quality_viral_genes_without_confident_aai"
    elif checkv_quality.lower() == "medium-quality" and host_genes == 0:
        classification = "virus"
        evidence_strength = "qualified"
        strength_basis = "checkv_medium_quality_viral_genes_without_host_genes"
    elif host_genes == 0:
        classification = "virus"
        evidence_strength = "weak"
        strength_basis = "checkv_low_quality_viral_genes_without_host_genes"
    else:
        classification = "virus"
        evidence_strength = "weak"
        strength_basis = "checkv_mixed_host_and_viral_genes_without_boundary"
    output = {
        "sample_id": sample_id,
        "sequence_id": sequence_id,
        "parent_sequence_id": parent_sequence_id,
        "record_type": record_type,
        "coordinates": coordinates,
        "tool": "checkv",
        "classification": classification,
        "score": completeness,
        "score_type": "checkv_completeness_percent" if completeness else "",
        "length": str(length),
        "topology": (
            "Provirus"
            if record_type == "provirus"
            else clean_missing(complete_call["prediction_type"])
            if is_confident_complete_call(complete_call)
            else ""
        ),
        "n_genes": clean_missing(quality["gene_count"]),
        "n_hallmarks": "",
        "evidence_strength": evidence_strength,
        "strength_basis": strength_basis,
        "checkv_quality": checkv_quality,
        "miuvig_quality": clean_missing(quality["miuvig_quality"]),
        "completeness": completeness,
        "completeness_method": clean_missing(quality["completeness_method"]),
        "contamination": clean_missing(quality["contamination"]),
        "provirus": clean_missing(quality["provirus"]),
        "proviral_length": clean_missing(quality["proviral_length"]),
        "viral_genes": clean_missing(quality["viral_genes"]),
        "host_genes": clean_missing(quality["host_genes"]),
        "complete_genome_prediction": (
            clean_missing(complete_call["prediction_type"]) if complete_call else ""
        ),
        "complete_genome_confidence": (
            clean_missing(complete_call["confidence_level"]) if complete_call else ""
        ),
        "aai_confidence": aai_confidence,
        "aai_identity": clean_missing(completeness_detail["aai_id"]),
        "aai_alignment_fraction": clean_missing(completeness_detail["aai_af"]),
        "aai_hit_count": clean_missing(completeness_detail["aai_num_hits"]),
        "hmm_completeness_lower": clean_missing(
            completeness_detail["hmm_completeness_lower"]
        ),
        "hmm_completeness_upper": clean_missing(
            completeness_detail["hmm_completeness_upper"]
        ),
        "hmm_hit_count": clean_missing(completeness_detail["hmm_num_hits"]),
        "warnings": clean_missing(quality["warnings"]),
    }
    output.update(unclassified_taxonomy(classification))
    return output


def run(args: argparse.Namespace) -> None:
    metadata = parse_metadata(args.run_metadata, args.sample_id, args.input_type)
    fasta_lengths = load_fasta_lengths(args.candidate_fasta)
    quality_rows = unique_rows(
        read_tsv(args.quality_summary, QUALITY_REQUIRED), "contig_id", "quality summary"
    )
    completeness_rows = unique_rows(
        read_tsv(args.completeness, COMPLETENESS_REQUIRED),
        "contig_id",
        "completeness summary",
    )
    contamination_rows = unique_rows(
        read_tsv(args.contamination, CONTAMINATION_REQUIRED),
        "contig_id",
        "contamination summary",
    )
    complete_calls = unique_rows(
        read_tsv(args.complete_genomes, COMPLETE_GENOMES_REQUIRED),
        "contig_id",
        "complete-genomes summary",
    )

    candidate_count = parse_nonnegative_integer(
        metadata["candidate_sequence_count"], "candidate sequence count"
    )
    metadata_quality_count = parse_nonnegative_integer(
        metadata["quality_summary_row_count"], "quality-summary row count"
    )
    if candidate_count != len(fasta_lengths) or metadata_quality_count != len(quality_rows):
        raise ValueError("CheckV metadata counts disagree with staged inputs")
    if set(quality_rows) != set(fasta_lengths):
        raise ValueError("CheckV quality summary and candidate FASTA identifiers disagree")
    if set(completeness_rows) != set(fasta_lengths):
        raise ValueError("CheckV completeness and candidate FASTA identifiers disagree")
    if set(contamination_rows) != set(fasta_lengths):
        raise ValueError("CheckV contamination and candidate FASTA identifiers disagree")
    if set(complete_calls).difference(fasta_lengths):
        raise ValueError("CheckV complete-genome calls contain an unknown candidate")

    evidence_rows: list[dict[str, str]] = []
    observed_proviruses = 0
    for contig_id, quality in quality_rows.items():
        parent_length = fasta_lengths[contig_id]
        for source_row, label in (
            (quality, "quality"),
            (completeness_rows[contig_id], "completeness"),
            (contamination_rows[contig_id], "contamination"),
        ):
            reported_length = parse_positive_integer(
                source_row["contig_length"], f"{label} length for {contig_id}"
            )
            if reported_length != parent_length:
                raise ValueError(f"CheckV {label} length disagrees for {contig_id}")

        complete_call = complete_calls.get(contig_id)
        if complete_call is not None:
            reported_length = parse_positive_integer(
                complete_call["contig_length"], f"complete-genome length for {contig_id}"
            )
            if reported_length != parent_length:
                raise ValueError(f"CheckV complete-genome length disagrees for {contig_id}")

        is_provirus = clean_missing(quality["provirus"]).lower() == "yes"
        contamination_is_provirus = (
            clean_missing(contamination_rows[contig_id]["provirus"]).lower() == "yes"
        )
        if is_provirus != contamination_is_provirus:
            raise ValueError(f"CheckV provirus flags disagree for {contig_id}")

        if is_provirus:
            observed_proviruses += 1
            regions = parse_regions(contamination_rows[contig_id], parent_length)
            reported_proviral_length = parse_positive_integer(
                quality["proviral_length"], f"proviral length for {contig_id}"
            )
            total_region_length = sum(end - start + 1 for start, end in regions)
            if total_region_length != reported_proviral_length:
                raise ValueError(f"CheckV proviral coordinates disagree for {contig_id}")
            for start, end in regions:
                sequence_id = f"{contig_id}|provirus_{start}_{end}"
                evidence_rows.append(
                    build_evidence_row(
                        args.sample_id,
                        quality,
                        completeness_rows[contig_id],
                        complete_call,
                        sequence_id,
                        contig_id,
                        "provirus",
                        f"{start}-{end}",
                        end - start + 1,
                    )
                )
        else:
            # Completeness availability must not discard measured gene content.
            # build_evidence_row keeps zero informative genes unclassified/weak;
            # an undetermined quality category is not itself negative evidence.
            evidence_rows.append(
                build_evidence_row(
                    args.sample_id,
                    quality,
                    completeness_rows[contig_id],
                    complete_call,
                    contig_id,
                    "",
                    "input_contig",
                    "",
                    parent_length,
                )
            )

    metadata_provirus_count = parse_nonnegative_integer(
        metadata["provirus_count"], "metadata provirus count"
    )
    if metadata_provirus_count != observed_proviruses:
        raise ValueError("CheckV metadata and standardized provirus counts disagree")
    observed_determined_count = sum(
        clean_missing(row["checkv_quality"]).lower() in DETERMINED_QUALITIES
        for row in quality_rows.values()
    )
    metadata_determined_count = parse_nonnegative_integer(
        metadata["determined_quality_count"], "metadata determined-quality count"
    )
    if metadata_determined_count != observed_determined_count:
        raise ValueError("CheckV metadata and determined-quality counts disagree")

    expected_status = (
        "skipped_no_discovery_candidates" if candidate_count == 0 else "completed"
    )
    if metadata["run_status"] != expected_status:
        raise ValueError("CheckV run status disagrees with the candidate count")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(evidence_rows)

    print(
        f"Standardized CheckV sample={args.sample_id} evidence={len(evidence_rows)} "
        f"proviruses={observed_proviruses}"
    )


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
