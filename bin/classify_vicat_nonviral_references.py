#!/usr/bin/env python3
"""Classify an NCBI protein catalogue by replicon for a viCAT decoy database."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
from collections import Counter, defaultdict
from pathlib import Path


OUTPUT_CLASSES = (
    "CELLULAR_CHROMOSOME",
    "CELLULAR_UNPLACED",
    "PLASMID",
    "PLASTID",
    "MITOCHONDRIAL",
    "SHARED_NONVIRAL",
)
REQUIRED_MANIFEST_COLUMNS = {"assembly_accession", "cellular_group"}
REQUIRED_FEATURE_COLUMNS = {
    "feature",
    "seq_type",
    "genomic_accession",
    "product_accession",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--metadata-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest is missing columns: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("Manifest contains no assemblies")
    accessions = [row["assembly_accession"].strip() for row in rows]
    if any(not accession for accession in accessions):
        raise ValueError("Manifest contains a blank assembly accession")
    if len(accessions) != len(set(accessions)):
        raise ValueError("Manifest contains duplicate assembly accessions")
    return rows


def reference_class(sequence_type: str) -> str | None:
    normalized = " ".join(sequence_type.strip().lower().replace("_", " ").split())
    if normalized in {"chromosome", "linkage group", "unlocalized scaffold on chromosome"}:
        return "CELLULAR_CHROMOSOME"
    if normalized == "unplaced scaffold":
        return "CELLULAR_UNPLACED"
    if normalized == "plasmid":
        return "PLASMID"
    if normalized in {"chloroplast", "plastid"}:
        return "PLASTID"
    if normalized in {"mitochondrion", "mitochondrial"}:
        return "MITOCHONDRIAL"
    return None


def load_feature_occurrences(path: Path) -> dict[str, list[tuple[str | None, str, str]]]:
    occurrences: dict[str, list[tuple[str | None, str, str]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = [field.lstrip("# ").strip() for field in next(reader)]
        columns = {name: index for index, name in enumerate(header)}
        missing = REQUIRED_FEATURE_COLUMNS - set(columns)
        if missing:
            raise ValueError(f"Feature table is missing columns: {sorted(missing)}")
        for fields in reader:
            if len(fields) < len(header) or fields[columns["feature"]] != "CDS":
                continue
            protein_id = fields[columns["product_accession"]].strip()
            if not protein_id or protein_id == "-":
                continue
            sequence_type = fields[columns["seq_type"]].strip()
            replicon = fields[columns["genomic_accession"]].strip()
            occurrences[protein_id].append((reference_class(sequence_type), sequence_type, replicon))
    return occurrences


def fasta_records(path: Path):
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:], []
            elif header is not None:
                sequence.append(line)
        if header is not None:
            yield header, "".join(sequence)


def classify_occurrences(occurrences: list[tuple[str | None, str, str]]) -> tuple[str | None, str]:
    if not occurrences:
        return None, "protein_absent_from_feature_table"
    unsupported = sorted({raw_type or "blank" for label, raw_type, _ in occurrences if label is None})
    if unsupported:
        return None, "unsupported_sequence_type:" + ",".join(unsupported)
    labels = {label for label, _, _ in occurrences if label is not None}
    if len(labels) == 1:
        return next(iter(labels)), ""
    return "SHARED_NONVIRAL", ""


def open_outputs(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    handles = {}
    for label in OUTPUT_CLASSES:
        final = output_dir / f"{label.lower()}.faa.gz"
        temporary = final.with_name(final.name + ".tmp")
        handles[label] = (final, temporary, gzip.open(temporary, "wt", encoding="utf-8", newline="\n"))
    metadata_final = output_dir / "vicat_nonviral_reference_metadata.tsv.gz"
    metadata_tmp = metadata_final.with_name(metadata_final.name + ".tmp")
    excluded_final = output_dir / "vicat_nonviral_exclusions.tsv.gz"
    excluded_tmp = excluded_final.with_name(excluded_final.name + ".tmp")
    return handles, (metadata_final, metadata_tmp), (excluded_final, excluded_tmp)


def main() -> None:
    args = parse_args()
    try:
        manifest = read_manifest(args.manifest)
        fasta_outputs, metadata_paths, excluded_paths = open_outputs(args.output_dir)
        metadata_handle = gzip.open(metadata_paths[1], "wt", encoding="utf-8", newline="")
        excluded_handle = gzip.open(excluded_paths[1], "wt", encoding="utf-8", newline="")
        metadata_columns = [
            "reference_id", "source_protein_id", "reference_class", "source_classes",
            "cellular_group", "source_accession", "replicon_accessions", "replicon_types",
            "provirus_flank_eligible", "organism_name", "taxid",
        ]
        exclusion_columns = [
            "source_accession", "cellular_group", "source_protein_id", "reason", "replicon_types",
        ]
        metadata_writer = csv.DictWriter(metadata_handle, fieldnames=metadata_columns, delimiter="\t", lineterminator="\n")
        exclusion_writer = csv.DictWriter(excluded_handle, fieldnames=exclusion_columns, delimiter="\t", lineterminator="\n")
        metadata_writer.writeheader()
        exclusion_writer.writeheader()

        class_counts = Counter()
        exclusion_counts = Counter()
        group_counts = Counter()
        input_proteins = 0
        assemblies = 0

        try:
            for row in manifest:
                accession = row["assembly_accession"].strip()
                group = row["cellular_group"].strip()
                feature = args.metadata_root / group / accession / f"{accession}.feature_table.txt.gz"
                proteins = args.package_root / "ncbi_dataset" / "data" / accession / "protein.faa"
                if not feature.is_file():
                    raise ValueError(f"Feature table not found: {feature}")
                if not proteins.is_file():
                    raise ValueError(f"Protein FASTA not found: {proteins}")
                occurrences = load_feature_occurrences(feature)
                seen = set()
                assemblies += 1
                for header, sequence in fasta_records(proteins):
                    protein_id = header.split()[0]
                    if protein_id in seen:
                        raise ValueError(f"Duplicate protein in {proteins}: {protein_id}")
                    seen.add(protein_id)
                    input_proteins += 1
                    source_occurrences = occurrences.get(protein_id, [])
                    label, reason = classify_occurrences(source_occurrences)
                    raw_types = sorted({raw_type or "blank" for _, raw_type, _ in source_occurrences})
                    if label is None:
                        exclusion_writer.writerow({
                            "source_accession": accession,
                            "cellular_group": group,
                            "source_protein_id": protein_id,
                            "reason": reason,
                            "replicon_types": ",".join(raw_types),
                        })
                        exclusion_counts[reason] += 1
                        continue
                    source_classes = sorted({item for item, _, _ in source_occurrences if item})
                    replicons = sorted({replicon for _, _, replicon in source_occurrences if replicon})
                    reference_id = f"NONVIRAL|{label}|{accession}|{protein_id}"
                    fasta_outputs[label][2].write(f">{reference_id}\n{sequence}\n")
                    metadata_writer.writerow({
                        "reference_id": reference_id,
                        "source_protein_id": protein_id,
                        "reference_class": label,
                        "source_classes": ",".join(source_classes),
                        "cellular_group": group,
                        "source_accession": accession,
                        "replicon_accessions": ",".join(replicons),
                        "replicon_types": ",".join(raw_types),
                        "provirus_flank_eligible": "true" if label == "CELLULAR_CHROMOSOME" else "false",
                        "organism_name": row.get("organism_name", "").strip(),
                        "taxid": (row.get("taxid") or row.get("species_taxid") or "").strip(),
                    })
                    class_counts[label] += 1
                    group_counts[group] += 1
        finally:
            metadata_handle.close()
            excluded_handle.close()
            for _, _, handle in fasta_outputs.values():
                handle.close()

        for final, temporary, _ in fasta_outputs.values():
            os.replace(temporary, final)
        os.replace(metadata_paths[1], metadata_paths[0])
        os.replace(excluded_paths[1], excluded_paths[0])

        summary = args.output_dir / "vicat_nonviral_classification_summary.tsv"
        with summary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["metric", "category", "value"])
            writer.writerow(["assembly_count", "all", assemblies])
            writer.writerow(["input_protein_count", "all", input_proteins])
            writer.writerow(["included_protein_count", "all", sum(class_counts.values())])
            writer.writerow(["excluded_protein_count", "all", sum(exclusion_counts.values())])
            for label in OUTPUT_CLASSES:
                writer.writerow(["reference_class_count", label, class_counts[label]])
            for group, count in sorted(group_counts.items()):
                writer.writerow(["cellular_group_included_count", group, count])
            for reason, count in sorted(exclusion_counts.items()):
                writer.writerow(["exclusion_reason_count", reason, count])

        print(f"Assemblies: {assemblies}")
        print(f"Input proteins: {input_proteins}")
        print(f"Included proteins: {sum(class_counts.values())}")
        print(f"Excluded proteins: {sum(exclusion_counts.values())}")
        print(f"Summary: {summary}")
    except (OSError, EOFError, StopIteration, UnicodeError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
