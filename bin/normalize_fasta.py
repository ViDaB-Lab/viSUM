#!/usr/bin/env python3

"""Normalize FASTA identifiers and write the viSUM record-provenance table."""

import argparse
import csv
from pathlib import Path
from typing import Iterator, TextIO


MAP_COLUMNS = [
    "sample_id",
    "sequence_id",
    "parent_sequence_id",
    "record_type",
    "original_id",
    "original_header",
    "length",
    "coordinates",
    "extraction_tool",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create stable viSUM sequence IDs while preserving FASTA provenance."
    )
    parser.add_argument("--input", required=True, type=Path, help="Input FASTA")
    parser.add_argument("--sample-id", required=True, help="Validated sample prefix")
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-map", required=True, type=Path)
    return parser.parse_args()


def read_fasta(handle: TextIO) -> Iterator[tuple[int, str, str]]:
    """Yield line number, full header, and whitespace-free sequence."""
    header: str | None = None
    header_line = 0
    sequence_parts: list[str] = []

    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.rstrip("\r\n")
        if not line.strip():
            continue

        if line.startswith(">"):
            if header is not None:
                yield header_line, header, "".join(sequence_parts)

            header = line[1:].strip()
            header_line = line_number
            sequence_parts = []

            if not header:
                raise ValueError(f"Empty FASTA header at line {line_number}")
        else:
            if header is None:
                raise ValueError(
                    f"Sequence data found before the first FASTA header at line {line_number}"
                )
            sequence_parts.append("".join(line.split()))

    if header is not None:
        yield header_line, header, "".join(sequence_parts)


def write_wrapped_fasta(handle: TextIO, sequence_id: str, sequence: str) -> None:
    handle.write(f">{sequence_id}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def normalize(args: argparse.Namespace) -> int:
    record_count = 0

    with (
        args.input.open("r", encoding="utf-8", errors="replace") as source,
        args.output_fasta.open("w", encoding="utf-8", newline="\n") as fasta_out,
        args.output_map.open("w", encoding="utf-8", newline="") as map_out,
    ):
        writer = csv.DictWriter(
            map_out,
            fieldnames=MAP_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for record_count, (header_line, header, sequence) in enumerate(
            read_fasta(source), start=1
        ):
            if not sequence:
                raise ValueError(
                    f"FASTA record '{header}' at line {header_line} has an empty sequence"
                )

            sequence_id = f"{args.sample_id}__c{record_count:06d}"
            original_id = header.split(maxsplit=1)[0]

            write_wrapped_fasta(fasta_out, sequence_id, sequence)
            writer.writerow(
                {
                    "sample_id": args.sample_id,
                    "sequence_id": sequence_id,
                    "parent_sequence_id": "",
                    "record_type": "input_contig",
                    "original_id": original_id,
                    "original_header": header,
                    "length": len(sequence),
                    "coordinates": "",
                    "extraction_tool": "",
                }
            )

    if record_count == 0:
        raise ValueError(f"No FASTA records found in '{args.input}'")

    return record_count


def main() -> None:
    args = parse_args()
    try:
        count = normalize(args)
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error

    print(f"Normalized {count} FASTA record(s) for sample '{args.sample_id}'.")


if __name__ == "__main__":
    main()
