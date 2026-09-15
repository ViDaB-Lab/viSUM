#!/usr/bin/env python3

"""Reconcile a benchmark truth TSV with a more complete curated truth TSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Table has no header: {path}")
        return list(reader.fieldnames), list(reader)


def unique_by_key(
    rows: list[dict[str, str]], key: str, path: Path
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = row.get(key, "").strip()
        if not identifier:
            raise ValueError(f"Blank {key} in {path}")
        if identifier in indexed:
            raise ValueError(f"Duplicate {key}={identifier} in {path}")
        indexed[identifier] = row
    return indexed


def equivalent_value(left: str, right: str) -> bool:
    if left == right:
        return True
    missing = {"", "NA", "N/A"}
    return left.strip().upper() in missing and right.strip().upper() in missing


def reconcile(
    target: Path,
    source: Path,
    output: Path,
    key: str = "benchmark_sequence_id",
    required_ids: set[str] | None = None,
) -> list[str]:
    columns, target_rows = read_tsv(target)
    source_columns, source_rows = read_tsv(source)
    if key not in columns:
        raise ValueError(f"Target lacks key column {key}: {target}")
    missing_columns = [column for column in columns if column not in source_columns]
    if missing_columns:
        raise ValueError(
            f"Source lacks target columns {missing_columns}: {source}"
        )

    target_by_id = unique_by_key(target_rows, key, target)
    source_by_id = unique_by_key(source_rows, key, source)
    required_ids = required_ids or set()
    comparison_ids = (
        sorted(required_ids.intersection(target_by_id))
        if required_ids else sorted(target_by_id)
    )
    for identifier in comparison_ids:
        target_row = target_by_id[identifier]
        source_row = source_by_id.get(identifier)
        if source_row is None:
            continue
        differences = [
            column for column in columns
            if not equivalent_value(
                target_row.get(column, ""), source_row.get(column, "")
            )
        ]
        if differences:
            raise ValueError(
                f"Conflicting existing truth row {identifier}; columns={differences}"
            )

    missing_ids = [
        row[key] for row in source_rows if row[key] not in target_by_id
    ]
    unavailable = sorted(required_ids.difference(source_by_id))
    if unavailable:
        raise ValueError(f"Required IDs absent from source: {unavailable}")
    already_present = sorted(required_ids.intersection(target_by_id))
    ids_to_add = [
        identifier for identifier in missing_ids
        if not required_ids or identifier in required_ids
    ]
    if required_ids and set(ids_to_add) != required_ids.difference(already_present):
        raise ValueError("Required truth-row reconciliation was incomplete")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(target_rows)
        writer.writerows(source_by_id[identifier] for identifier in ids_to_add)
    return ids_to_add


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Append missing truth rows from a curated TSV while rejecting "
            "duplicate IDs, schema mismatch, and conflicting existing rows."
        )
    )
    parser.add_argument("--target", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--key", default="benchmark_sequence_id")
    parser.add_argument("--require-id", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        added = reconcile(
            args.target, args.source, args.output, args.key, set(args.require_id)
        )
    except (OSError, ValueError, csv.Error) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(f"Reconciled truth rows: added={len(added)} ids={','.join(added) or 'none'}")


if __name__ == "__main__":
    main()
