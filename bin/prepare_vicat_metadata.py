#!/usr/bin/env python3
"""Create the compact MetaVR metadata table used by the viCAT database build."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


REQUIRED_COLUMNS = [
    "uvig",
    "votu",
    "ictv_taxonomy",
    "ictv_taxonomy_method",
]
OPTIONAL_COLUMNS = [
    "genome_type",
    "viral_confidence",
    "quality",
    "genomad_score",
    "source",
    "length",
    "cds_count",
    "n_virus_hallmarks",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-rows", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"MetaVR metadata file not found: {args.input}")
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {args.output}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with gzip.open(args.input, "rt", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("MetaVR metadata has no header")
        missing = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise SystemExit("MetaVR metadata is missing columns: " + ", ".join(missing))

        retained = REQUIRED_COLUMNS + [
            name for name in OPTIONAL_COLUMNS if name in reader.fieldnames
        ]
        with gzip.open(args.output, "wt", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(
                target,
                fieldnames=retained,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            for row in reader:
                writer.writerow(row)
                row_count += 1

    if args.expected_rows is not None and row_count != args.expected_rows:
        args.output.unlink(missing_ok=True)
        raise SystemExit(
            f"MetaVR metadata row count mismatch: observed {row_count}, "
            f"expected {args.expected_rows}"
        )
    print(f"MetaVR metadata rows: {row_count}")
    print(f"Compact metadata: {args.output}")


if __name__ == "__main__":
    main()
