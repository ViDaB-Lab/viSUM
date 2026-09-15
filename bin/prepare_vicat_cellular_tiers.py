#!/usr/bin/env python3
"""Create deterministic, nested cellular-protein tiers for viCAT benchmarking."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import defaultdict
from pathlib import Path


REQUIRED_COLUMNS = {"protein_id", "cellular_group", "source_accession"}


def open_text(path: Path, mode: str):
    return gzip.open(path, mode + "t", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode, encoding="utf-8", newline="")


def parse_sizes(value: str) -> list[tuple[str, int]]:
    tiers = []
    for item in value.split(","):
        name, raw_size = item.split(":", 1)
        name = name.strip()
        size = int(raw_size)
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or size < 1:
            raise argparse.ArgumentTypeError(f"Invalid tier: {item}")
        tiers.append((name, size))
    if len({name for name, _ in tiers}) != len(tiers):
        raise argparse.ArgumentTypeError("Tier names must be unique")
    if [size for _, size in tiers] != sorted(size for _, size in tiers):
        raise argparse.ArgumentTypeError("Tier sizes must increase")
    return tiers


def stable_key(seed: int, row: dict[str, str]) -> str:
    token = "\t".join((str(seed), row["cellular_group"], row["source_accession"], row["protein_id"]))
    return hashlib.sha256(token.encode()).hexdigest()


def read_metadata(path: Path, seed: int) -> tuple[list[str], list[dict[str, str]]]:
    with open_text(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames or []))
            raise ValueError(f"Cellular metadata is missing columns: {missing}")
        rows = []
        seen = set()
        for row in reader:
            protein_id = row["protein_id"].strip()
            group = row["cellular_group"].strip()
            source = row["source_accession"].strip()
            if not protein_id or not group or not source:
                raise ValueError("protein_id, cellular_group, and source_accession cannot be blank")
            if protein_id in seen:
                raise ValueError(f"Duplicate protein_id in cellular metadata: {protein_id}")
            seen.add(protein_id)
            row = dict(row)
            row.update(protein_id=protein_id, cellular_group=group, source_accession=source)
            row["_key"] = stable_key(seed, row)
            rows.append(row)
    return list(reader.fieldnames), rows


def balanced_order(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Round-robin groups after stable within-group hashing; every tier is a prefix."""
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["cellular_group"]].append(row)
    for values in groups.values():
        values.sort(key=lambda row: (row["_key"], row["protein_id"]))
    names = sorted(groups)
    positions = {name: 0 for name in names}
    ordered = []
    active = names
    while active:
        next_active = []
        for name in active:
            position = positions[name]
            values = groups[name]
            if position < len(values):
                ordered.append(values[position])
                positions[name] += 1
            if positions[name] < len(values):
                next_active.append(name)
        active = next_active
    return ordered


def fasta_records(path: Path):
    with open_text(path, "r") as handle:
        header = None
        sequence = []
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            elif header is not None:
                sequence.append(line.strip())
        if header is not None:
            yield header, "".join(sequence)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cellular-proteins", required=True, type=Path)
    parser.add_argument("--cellular-metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--tiers", default="small:250000,medium:750000,large:1500000")
    parser.add_argument("--seed", type=int, default=41)
    args = parser.parse_args()

    tiers = parse_sizes(args.tiers)
    fields, rows = read_metadata(args.cellular_metadata, args.seed)
    largest = tiers[-1][1]
    if len(rows) < largest:
        raise ValueError(f"Requested {largest:,} proteins but metadata contains only {len(rows):,}")
    selected = balanced_order(rows)[:largest]
    rank = {row["protein_id"]: index + 1 for index, row in enumerate(selected)}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_handles = {}
    metadata_handles = {}
    writers = {}
    try:
        for name, _ in tiers:
            fasta_handles[name] = gzip.open(args.output_dir / f"cellular_{name}.faa.gz", "wt", encoding="utf-8")
            metadata_handles[name] = gzip.open(args.output_dir / f"cellular_{name}.metadata.tsv.gz", "wt", encoding="utf-8", newline="")
            writers[name] = csv.DictWriter(metadata_handles[name], fieldnames=fields + ["tier_rank"], delimiter="\t", lineterminator="\n")
            writers[name].writeheader()

        found = set()
        for header, sequence in fasta_records(args.cellular_proteins):
            protein_id = header.split()[0]
            tier_rank = rank.get(protein_id)
            if tier_rank is None:
                continue
            if protein_id in found:
                raise ValueError(f"Duplicate selected protein in FASTA: {protein_id}")
            found.add(protein_id)
            for name, size in tiers:
                if tier_rank <= size:
                    fasta_handles[name].write(f">{header}\n{sequence}\n")

        missing = set(rank) - found
        if missing:
            examples = ", ".join(sorted(missing)[:5])
            raise ValueError(f"FASTA is missing {len(missing):,} selected proteins: {examples}")

        selected_by_id = {row["protein_id"]: row for row in selected}
        for protein_id, tier_rank in sorted(rank.items(), key=lambda item: item[1]):
            row = {key: value for key, value in selected_by_id[protein_id].items() if not key.startswith("_")}
            row["tier_rank"] = tier_rank
            for name, size in tiers:
                if tier_rank <= size:
                    writers[name].writerow(row)
    finally:
        for handle in fasta_handles.values():
            handle.close()
        for handle in metadata_handles.values():
            handle.close()

    with (args.output_dir / "cellular_tier_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tier", "requested_proteins", "selected_proteins", "seed", "selection_method"])
        for name, size in tiers:
            writer.writerow([name, size, size, args.seed, "stable_hash_group_round_robin_prefix"])
    print(f"Created {len(tiers)} nested tiers in {args.output_dir}")


if __name__ == "__main__":
    main()
