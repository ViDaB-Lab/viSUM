#!/usr/bin/env python3

"""Window oversized nucleotide sequences before TEsorter element-mode analysis."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DEFAULT_MAX_LENGTH = 270_000
DEFAULT_OVERLAP = 30_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-map", required=True, type=Path)
    parser.add_argument("--max-length", type=int, default=DEFAULT_MAX_LENGTH)
    parser.add_argument("--overlap", type=int, default=DEFAULT_OVERLAP)
    return parser.parse_args()


def read_fasta(path: Path):
    identifier = None
    description = ""
    sequence: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, description, "".join(sequence)
                description = line[1:].strip()
                identifier = description.split()[0] if description else ""
                if not identifier or identifier in seen:
                    raise ValueError(f"Missing or duplicate FASTA identifier: {identifier}")
                seen.add(identifier)
                sequence = []
            else:
                if identifier is None:
                    raise ValueError("FASTA sequence encountered before its header")
                sequence.append(line)
    if identifier is not None:
        yield identifier, description, "".join(sequence)


def write_record(handle, identifier: str, sequence: str) -> None:
    handle.write(f">{identifier}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + "\n")


def main() -> None:
    args = parse_args()
    if args.max_length <= 0:
        raise ValueError("--max-length must be positive")
    if args.overlap < 0 or args.overlap >= args.max_length:
        raise ValueError("--overlap must be nonnegative and smaller than --max-length")

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_map.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    input_count = 0
    window_count = 0
    split_count = 0
    step = args.max_length - args.overlap

    with args.output_fasta.open("w", encoding="utf-8", newline="\n") as output:
        for sequence_id, _description, sequence in read_fasta(args.input):
            input_count += 1
            if not sequence:
                raise ValueError(f"Empty FASTA sequence: {sequence_id}")
            length = len(sequence)
            if length <= args.max_length:
                windows = [(sequence_id, 0, length)]
            else:
                split_count += 1
                windows = []
                start = 0
                number = 1
                while start < length:
                    end = min(start + args.max_length, length)
                    window_id = (
                        f"{sequence_id}__tesorter_window_{number:04d}_"
                        f"{start + 1}_{end}"
                    )
                    windows.append((window_id, start, end))
                    if end == length:
                        break
                    start += step
                    number += 1

            for window_id, start, end in windows:
                fragment = sequence[start:end]
                write_record(output, window_id, fragment)
                window_count += 1
                rows.append({
                    "tesorter_sequence_id": window_id,
                    "sequence_id": sequence_id,
                    "window_start": str(start + 1),
                    "window_end": str(end),
                    "window_length": str(end - start),
                    "original_length": str(length),
                    "was_split": "true" if length > args.max_length else "false",
                })

    columns = [
        "tesorter_sequence_id", "sequence_id", "window_start", "window_end",
        "window_length", "original_length", "was_split",
    ]
    with args.output_map.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"TESORTER_INPUT sequences={input_count} windows={window_count} "
        f"split_sequences={split_count} max_length={args.max_length} "
        f"overlap={args.overlap}"
    )


if __name__ == "__main__":
    main()
