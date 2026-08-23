#!/usr/bin/env python3
"""Build viCAT cellular FASTA and metadata from an NCBI Datasets package."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-metadata", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = arguments()
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"assembly_accession", "cellular_group"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest is missing columns: {sorted(missing)}")
        rows = [dict(row) for row in reader]

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    proteins = assemblies = 0
    seen: set[str] = set()
    absent: list[str] = []

    with gzip.open(args.output_fasta, "wt", encoding="utf-8", newline="\n") as fasta_out, gzip.open(
        args.output_metadata, "wt", encoding="utf-8", newline=""
    ) as metadata_out:
        writer = csv.writer(metadata_out, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "cellular_group", "source_accession"])
        for row in rows:
            accession = row["assembly_accession"].strip()
            group = row["cellular_group"].strip()
            protein_path = args.package_root / "ncbi_dataset" / "data" / accession / "protein.faa"
            if not protein_path.is_file():
                absent.append(accession)
                continue
            assemblies += 1
            header: str | None = None
            sequence: list[str] = []

            def emit() -> None:
                nonlocal proteins, header, sequence
                if header is None:
                    return
                original = header.split()[0]
                protein_id = f"CELLULAR|{accession}|{original}"
                if protein_id in seen:
                    raise ValueError(f"Duplicate generated protein ID: {protein_id}")
                seen.add(protein_id)
                fasta_out.write(f">{protein_id}\n{''.join(sequence)}\n")
                writer.writerow([protein_id, group, accession])
                proteins += 1

            with protein_path.open(encoding="utf-8") as protein_handle:
                for raw in protein_handle:
                    line = raw.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        emit()
                        header, sequence = line[1:], []
                    else:
                        sequence.append(line)
                emit()

    if absent:
        preview = ", ".join(absent[:10])
        raise ValueError(f"Missing protein.faa for {len(absent)} assemblies: {preview}")
    if proteins == 0:
        raise ValueError("No proteins were written")
    print(f"Assemblies: {assemblies}")
    print(f"Proteins: {proteins}")
    print(f"FASTA: {args.output_fasta}")
    print(f"Metadata: {args.output_metadata}")


if __name__ == "__main__":
    main()
