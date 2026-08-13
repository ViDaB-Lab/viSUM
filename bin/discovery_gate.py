#!/usr/bin/env python3

"""Route contigs with qualified viral evidence into viSUM refinement."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Iterator, TextIO


AUDIT_COLUMNS = [
    "sample_id",
    "sequence_id",
    "input_type",
    "length",
    "discovery_status",
    "advance_to_refinement",
    "viral_tool_count",
    "viral_tools",
    "viral_evidence_row_count",
    "cellular_tool_count",
    "cellular_tools",
    "cellular_evidence_row_count",
    "plasmid_tool_count",
    "plasmid_tools",
    "plasmid_evidence_row_count",
    "conflicting_origin_evidence",
    "decision_reason",
]

SUMMARY_COLUMNS = [
    "sample_id",
    "input_type",
    "input_sequence_count",
    "candidate_sequence_count",
    "noncandidate_sequence_count",
    "viral_count",
    "likely_viral_count",
    "ambiguous_count",
    "likely_nonviral_count",
    "unresolved_count",
    "evidence_file_count",
]

REQUIRED_EVIDENCE_COLUMNS = {
    "sample_id",
    "sequence_id",
    "parent_sequence_id",
    "tool",
    "classification",
}
ALLOWED_CLASSIFICATIONS = {"virus", "cellular", "plasmid"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate standardized discovery evidence at input-contig level and "
            "write the candidate FASTA used for refinement."
        )
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument(
        "--evidence",
        nargs="*",
        default=[],
        type=Path,
        help="Standardized sparse evidence TSV files for this sample.",
    )
    parser.add_argument("--output-candidates", required=True, type=Path)
    parser.add_argument("--output-noncandidates", required=True, type=Path)
    parser.add_argument("--output-audit", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return list(reader)


def parse_positive_integer(value: str, label: str) -> int:
    try:
        integer = int(value)
    except ValueError as error:
        raise ValueError(f"{label} is not an integer: {value}") from error
    if integer < 1:
        raise ValueError(f"{label} must be positive: {value}")
    return integer


def load_header_map(path: Path, sample_id: str) -> dict[str, dict[str, str]]:
    rows = read_tsv(
        path,
        {"sample_id", "sequence_id", "parent_sequence_id", "record_type", "length"},
    )
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
        if row["parent_sequence_id"].strip():
            raise ValueError(
                f"Discovery gate expects input-contig records, not child record {sequence_id}"
            )
        if row["record_type"].strip() != "input_contig":
            raise ValueError(f"Unexpected header-map record type for {sequence_id}")
        parse_positive_integer(row["length"], f"length for {sequence_id}")
        records[sequence_id] = row
    if not records:
        raise ValueError("Header map contains no input contigs")
    return records


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


def load_fasta(path: Path, header_records: dict[str, dict[str, str]]) -> dict[str, str]:
    records: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for sequence_id, sequence in read_fasta(handle):
            if sequence_id in records:
                raise ValueError(f"Duplicate FASTA sequence ID: {sequence_id}")
            if sequence_id not in header_records:
                raise ValueError(f"FASTA sequence is absent from header map: {sequence_id}")
            if not sequence:
                raise ValueError(f"FASTA sequence is empty: {sequence_id}")
            expected_length = parse_positive_integer(
                header_records[sequence_id]["length"], f"length for {sequence_id}"
            )
            if len(sequence) != expected_length:
                raise ValueError(
                    f"FASTA and header-map lengths disagree for {sequence_id}"
                )
            records[sequence_id] = sequence
    if set(records) != set(header_records):
        missing = sorted(set(header_records).difference(records))
        raise ValueError(f"Header-map sequences absent from FASTA: {missing[:5]}")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for sequence_id, sequence in records:
            handle.write(f">{sequence_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def root_sequence_id(
    row: dict[str, str], header_records: dict[str, dict[str, str]], path: Path
) -> str:
    sequence_id = row["sequence_id"].strip()
    parent_id = row["parent_sequence_id"].strip()
    root_id = parent_id or sequence_id
    if not root_id:
        raise ValueError(f"Evidence row has no sequence identity in {path}")
    if root_id not in header_records:
        raise ValueError(
            f"Evidence root sequence is absent from the header map in {path}: {root_id}"
        )
    if parent_id and sequence_id == parent_id:
        raise ValueError(f"Evidence child equals its parent in {path}: {sequence_id}")
    return root_id


def load_evidence(
    paths: list[Path],
    sample_id: str,
    header_records: dict[str, dict[str, str]],
) -> dict[str, dict[str, list[str]]]:
    evidence: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    seen_paths: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen_paths:
            raise ValueError(f"Evidence file supplied more than once: {path}")
        seen_paths.add(resolved)
        for row in read_tsv(path, REQUIRED_EVIDENCE_COLUMNS):
            if row["sample_id"] != sample_id:
                raise ValueError(
                    f"Evidence sample '{row['sample_id']}' in {path} does not match "
                    f"'{sample_id}'"
                )
            tool = row["tool"].strip()
            classification = row["classification"].strip()
            if not tool:
                raise ValueError(f"Evidence row has an empty tool in {path}")
            if classification not in ALLOWED_CLASSIFICATIONS:
                raise ValueError(
                    f"Unsupported evidence classification in {path}: {classification}"
                )
            root_id = root_sequence_id(row, header_records, path)
            evidence[root_id][classification].append(tool)
    return evidence


def classify(tools_by_classification: dict[str, list[str]]) -> tuple[str, bool, str]:
    viral_tools = set(tools_by_classification.get("virus", []))
    cellular_tools = set(tools_by_classification.get("cellular", []))
    plasmid_tools = set(tools_by_classification.get("plasmid", []))
    opposing_tools = cellular_tools | plasmid_tools

    if viral_tools and opposing_tools:
        return (
            "ambiguous",
            True,
            "qualified viral evidence conflicts with cellular or plasmid evidence",
        )
    if len(viral_tools) >= 2:
        return "viral", True, "qualified viral evidence from multiple tools"
    if len(viral_tools) == 1:
        return "likely_viral", True, "qualified viral evidence from one tool"
    if opposing_tools:
        return (
            "likely_nonviral",
            False,
            "qualified cellular or plasmid evidence with no qualified viral evidence",
        )
    return "unresolved", False, "no qualified discovery evidence"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> None:
    header_records = load_header_map(args.header_map, args.sample_id)
    fasta_records = load_fasta(args.fasta, header_records)
    evidence = load_evidence(args.evidence, args.sample_id, header_records)

    audit_rows: list[dict[str, object]] = []
    candidate_records: list[tuple[str, str]] = []
    noncandidate_records: list[tuple[str, str]] = []
    status_counts: dict[str, int] = defaultdict(int)

    for sequence_id, sequence in fasta_records.items():
        tools_by_classification = evidence.get(sequence_id, {})
        status, advance, reason = classify(tools_by_classification)
        status_counts[status] += 1

        viral_calls = tools_by_classification.get("virus", [])
        cellular_calls = tools_by_classification.get("cellular", [])
        plasmid_calls = tools_by_classification.get("plasmid", [])
        viral_tools = sorted(set(viral_calls))
        cellular_tools = sorted(set(cellular_calls))
        plasmid_tools = sorted(set(plasmid_calls))
        conflicting = bool(viral_tools and (cellular_tools or plasmid_tools))
        audit_rows.append(
            {
                "sample_id": args.sample_id,
                "sequence_id": sequence_id,
                "input_type": args.input_type,
                "length": len(sequence),
                "discovery_status": status,
                "advance_to_refinement": "true" if advance else "false",
                "viral_tool_count": len(viral_tools),
                "viral_tools": ",".join(viral_tools),
                "viral_evidence_row_count": len(viral_calls),
                "cellular_tool_count": len(cellular_tools),
                "cellular_tools": ",".join(cellular_tools),
                "cellular_evidence_row_count": len(cellular_calls),
                "plasmid_tool_count": len(plasmid_tools),
                "plasmid_tools": ",".join(plasmid_tools),
                "plasmid_evidence_row_count": len(plasmid_calls),
                "conflicting_origin_evidence": "true" if conflicting else "false",
                "decision_reason": reason,
            }
        )
        target = candidate_records if advance else noncandidate_records
        target.append((sequence_id, sequence))

    write_fasta(args.output_candidates, candidate_records)
    write_fasta(args.output_noncandidates, noncandidate_records)
    write_tsv(args.output_audit, AUDIT_COLUMNS, audit_rows)
    write_tsv(
        args.output_summary,
        SUMMARY_COLUMNS,
        [
            {
                "sample_id": args.sample_id,
                "input_type": args.input_type,
                "input_sequence_count": len(fasta_records),
                "candidate_sequence_count": len(candidate_records),
                "noncandidate_sequence_count": len(noncandidate_records),
                "viral_count": status_counts["viral"],
                "likely_viral_count": status_counts["likely_viral"],
                "ambiguous_count": status_counts["ambiguous"],
                "likely_nonviral_count": status_counts["likely_nonviral"],
                "unresolved_count": status_counts["unresolved"],
                "evidence_file_count": len(args.evidence),
            }
        ],
    )
    print(
        f"Discovery gate sample={args.sample_id} input={len(fasta_records)} "
        f"candidates={len(candidate_records)} noncandidates={len(noncandidate_records)}"
    )


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
