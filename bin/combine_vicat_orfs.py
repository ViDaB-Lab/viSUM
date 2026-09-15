#!/usr/bin/env python3
"""Combine pyrodigal-gv/rv calls while retaining caller and locus metadata."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Prediction:
    parent: str
    start: int
    end: int
    strand: str
    translation: str
    translation_table: str
    caller: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--gv-proteins", required=True, type=Path)
    parser.add_argument("--gv-gff", required=True, type=Path)
    parser.add_argument("--rv-proteins", type=Path)
    parser.add_argument("--rv-gff", type=Path)
    parser.add_argument("--output-proteins", required=True, type=Path)
    parser.add_argument("--output-map", required=True, type=Path)
    return parser.parse_args()


def parse_attributes(raw: str) -> dict[str, str]:
    output = {}
    for token in raw.strip().split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            output[key] = value
    return output


def load_gff(path: Path) -> dict[str, tuple[str, int, int, str, str]]:
    features = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attributes = parse_attributes(fields[8])
            feature_id = attributes.get("ID")
            if not feature_id:
                raise ValueError(f"GFF CDS is missing ID in {path}: {line.rstrip()}")
            if feature_id in features:
                raise ValueError(f"Duplicate GFF feature ID in {path}: {feature_id}")
            features[feature_id] = (
                fields[0],
                int(fields[3]),
                int(fields[4]),
                fields[6],
                attributes.get("transl_table", ""),
            )
    return features


def fasta_records(path: Path):
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            elif header is not None:
                sequence.append(line.strip())
    if header is not None:
        yield header, "".join(sequence)


def load_predictions(proteins: Path, gff: Path, caller: str) -> list[Prediction]:
    features = load_gff(gff)
    predictions = []
    for header, translation in fasta_records(proteins):
        parts = [part.strip() for part in header.split(" # ")]
        attributes = parse_attributes(parts[-1]) if parts else {}
        # pyrodigal-gv 0.3.2 can emit a short attribute ID (for example
        # ID=1_1) in the protein description while the GFF uses the full
        # FASTA record identifier (for example sample_c000001_1). Prefer the
        # record identifier and retain the attribute ID for compatible CLIs.
        candidate_ids = [parts[0].split()[0] if parts else "", attributes.get("ID", "")]
        feature_id = next((value for value in candidate_ids if value in features), None)
        if feature_id is None:
            # Final compatibility fallback: match the parent sequence and
            # coordinates encoded in the Prodigal-style protein header.
            try:
                parent_id = candidate_ids[0].rsplit("_", 1)[0]
                header_start, header_end = int(parts[1]), int(parts[2])
                header_strand = "+" if int(parts[3]) == 1 else "-"
            except (IndexError, ValueError):
                parent_id, header_start, header_end, header_strand = "", -1, -1, ""
            matches = [
                key for key, value in features.items()
                if value[:4] == (parent_id, header_start, header_end, header_strand)
            ]
            if len(matches) == 1:
                feature_id = matches[0]
        if feature_id is None:
            raise ValueError(
                f"Could not connect protein header to GFF feature in {proteins}: {header}"
            )
        parent, start, end, strand, table = features[feature_id]
        if len(parts) >= 4:
            expected_strand = "+" if parts[3] == "1" else "-"
            if (int(parts[1]), int(parts[2]), expected_strand) != (start, end, strand):
                raise ValueError(
                    f"Protein/GFF coordinates disagree for {feature_id} in {proteins}"
                )
        predictions.append(
            Prediction(parent, start, end, strand, translation.rstrip("*"), table, caller)
        )
    if len(predictions) != len(features):
        raise ValueError(
            f"Protein/GFF count mismatch for {caller}: "
            f"{len(predictions)} proteins and {len(features)} CDS features"
        )
    return predictions


def main() -> None:
    args = parse_args()
    predictions = load_predictions(args.gv_proteins, args.gv_gff, "pyrodigal-gv")
    if args.input_type == "rna":
        if args.rv_proteins is None or args.rv_gff is None:
            raise SystemExit("RNA inputs require both --rv-proteins and --rv-gff")
        predictions.extend(load_predictions(args.rv_proteins, args.rv_gff, "pyrodigal-rv"))

    # Exact coordinate/strand/translation matches are one prediction supported
    # by both callers. Different starts, frames, or translations remain visible.
    collapsed: dict[tuple, set[str]] = {}
    representative: dict[tuple, Prediction] = {}
    for prediction in predictions:
        key = (
            prediction.parent,
            prediction.start,
            prediction.end,
            prediction.strand,
            prediction.translation,
        )
        collapsed.setdefault(key, set()).add(prediction.caller)
        representative[key] = prediction

    ordered = sorted(
        representative.items(),
        key=lambda item: (
            item[1].parent,
            item[1].start,
            item[1].end,
            item[1].strand,
            item[1].translation,
        ),
    )
    args.output_proteins.parent.mkdir(parents=True, exist_ok=True)
    args.output_map.parent.mkdir(parents=True, exist_ok=True)
    with args.output_proteins.open("w", encoding="utf-8") as proteins, args.output_map.open(
        "w", encoding="utf-8", newline=""
    ) as mapping:
        columns = [
            "sample_id", "orf_id", "sequence_id", "start", "end", "strand",
            "protein_length", "translation_table", "callers",
        ]
        writer = csv.DictWriter(mapping, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, (key, prediction) in enumerate(ordered, start=1):
            orf_id = f"{args.sample_id}_vicat_orf{index:09d}"
            callers = ",".join(sorted(collapsed[key]))
            proteins.write(f">{orf_id}\n{prediction.translation}\n")
            writer.writerow(
                {
                    "sample_id": args.sample_id,
                    "orf_id": orf_id,
                    "sequence_id": prediction.parent,
                    "start": prediction.start,
                    "end": prediction.end,
                    "strand": prediction.strand,
                    "protein_length": len(prediction.translation),
                    "translation_table": prediction.translation_table,
                    "callers": callers,
                }
            )

    print(
        f"Combined viCAT ORFs: raw={len(predictions)} unique={len(ordered)} "
        f"input_type={args.input_type}"
    )


if __name__ == "__main__":
    main()
