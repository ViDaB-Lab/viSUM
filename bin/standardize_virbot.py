#!/usr/bin/env python3

"""Convert VirBot calls into sparse, ICTV-validated viSUM evidence."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from evidence_schema import CORE_EVIDENCE_COLUMNS, TAXONOMY_COLUMNS, unclassified_taxonomy


ICTV_RANKS = (
    "Realm",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
)
SCORE_REQUIRED = {
    "Contig_acc",
    "RNA-viral_gene_content",
    "Encoded_proteins_num",
    "Likely_taxa",
}
METADATA_REQUIRED = {
    "sample_id",
    "input_type",
    "input_sequence_count",
    "positive_sequence_count",
    "sensitive_mode",
    "taxa_mode",
    "virbot_version",
    "virbot_revision",
    "run_status",
    "score_file",
    "virus_fasta",
}
OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "evidence_strength",
    "strength_basis",
]
VIRBOT_MINIMUM_GENE_FRACTION = 0.0625


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert VirBot calls into ICTV-validated viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("rna",))
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--score-table", required=True, type=Path)
    parser.add_argument("--virus-fasta", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--ictv-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def clean_missing(value: str | None) -> str:
    value = "" if value is None else value.strip()
    return "" if value.lower() in {"", "na", "nan", "none"} else value


def clean_taxon(value: str | None) -> str:
    return clean_missing(value).rstrip(".").strip()


def read_table(
    path: Path,
    required_columns: set[str],
    delimiter: str,
) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Table has no header: {path}")
        missing = required_columns.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return list(reader)


def parse_integer(
    value: str | None,
    label: str,
    sequence_id: str,
    minimum: int = 0,
) -> int:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing {label} for {sequence_id}")
    numeric = float(cleaned)
    if not math.isfinite(numeric) or not numeric.is_integer():
        raise ValueError(f"{label} is not an integer for {sequence_id}: {cleaned}")
    integer = int(numeric)
    if integer < minimum:
        raise ValueError(
            f"{label} must be at least {minimum} for {sequence_id}: {cleaned}"
        )
    return integer


def parse_score(value: str | None, sequence_id: str) -> tuple[float, str]:
    cleaned = clean_missing(value)
    if not cleaned:
        raise ValueError(f"Missing RNA-viral gene fraction for {sequence_id}")
    numeric = float(cleaned)
    if not math.isfinite(numeric) or not 0.0 <= numeric <= 1.0:
        raise ValueError(
            f"RNA-viral gene fraction outside 0-1 for {sequence_id}: {cleaned}"
        )
    if numeric < VIRBOT_MINIMUM_GENE_FRACTION:
        raise ValueError(
            f"VirBot-positive sequence is below VirBot's internal cutoff for "
            f"{sequence_id}: {cleaned}"
        )
    return numeric, cleaned


def load_header_map(path: Path, sample_id: str) -> dict[str, dict[str, str]]:
    rows = read_table(
        path,
        {"sample_id", "sequence_id", "record_type", "length"},
        "\t",
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
        parse_integer(row["length"], "header-map length", sequence_id, minimum=1)
        records[sequence_id] = row
    if not records:
        raise ValueError("Header map contains no sequence records")
    return records


def load_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    sequence_id: str | None = None
    length = 0

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_id is not None:
                    if length < 1:
                        raise ValueError(f"Empty FASTA record: {sequence_id}")
                    lengths[sequence_id] = length
                sequence_id = line[1:].split(maxsplit=1)[0]
                if not sequence_id:
                    raise ValueError(f"Empty FASTA identifier at line {line_number}")
                if sequence_id in lengths:
                    raise ValueError(f"Duplicate FASTA identifier: {sequence_id}")
                length = 0
            else:
                if sequence_id is None:
                    raise ValueError("Sequence data precedes the first FASTA header")
                length += len(line)

    if sequence_id is not None:
        if length < 1:
            raise ValueError(f"Empty FASTA record: {sequence_id}")
        lengths[sequence_id] = length
    return lengths


def load_metadata(
    path: Path,
    sample_id: str,
    input_type: str,
) -> dict[str, str]:
    rows = read_table(path, METADATA_REQUIRED, "\t")
    if len(rows) != 1:
        raise ValueError(f"Expected one VirBot metadata row in {path}; found {len(rows)}")
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError("VirBot metadata does not match the requested sample and type")
    parse_integer(row["input_sequence_count"], "input sequence count", sample_id)
    parse_integer(row["positive_sequence_count"], "positive sequence count", sample_id)
    if row["sensitive_mode"].strip().lower() not in {"true", "false"}:
        raise ValueError("VirBot metadata has an invalid sensitive_mode")
    if row["taxa_mode"].strip().upper() not in {"TOP", "LCA"}:
        raise ValueError("VirBot metadata has an invalid taxa_mode")
    if row["run_status"] not in {"completed", "completed_no_virus_calls"}:
        raise ValueError(f"Unrecognized VirBot run status: {row['run_status']}")
    for column in ("virbot_version", "virbot_revision", "score_file", "virus_fasta"):
        if not clean_missing(row[column]):
            raise ValueError(f"VirBot metadata contains an empty {column}")
    return row


def load_ictv_taxonomy(
    path: Path,
) -> dict[str, list[tuple[int, tuple[str, ...]]]]:
    rows = read_table(path, set(ICTV_RANKS), ",")
    occurrences: dict[str, dict[int, set[tuple[str, ...]]]] = defaultdict(
        lambda: defaultdict(set)
    )

    for row in rows:
        lineage = tuple(clean_taxon(row[rank]) for rank in ICTV_RANKS)
        for rank_index, taxon in enumerate(lineage):
            if taxon:
                occurrences[taxon.casefold()][rank_index].add(lineage)

    if not occurrences:
        raise ValueError(f"ICTV taxonomy table contains no taxon names: {path}")

    # Collapse each taxon to the lineage shared by all of its current ICTV
    # species. This makes per-sequence resolution constant-time and prevents
    # ambiguous lower ranks from being inferred.
    index: dict[str, list[tuple[int, tuple[str, ...]]]] = {}
    for taxon, ranks in occurrences.items():
        index[taxon] = []
        for rank_index, lineages in ranks.items():
            consensus = tuple(
                next(iter(values)) if len(values) == 1 else ""
                for values in (
                    {lineage[index] for lineage in lineages if lineage[index]}
                    for index in range(len(ICTV_RANKS))
                )
            )
            index[taxon].append((rank_index, consensus))
    return index


def validated_taxonomy(
    raw_lineage: str | None,
    ictv_index: dict[str, list[tuple[int, tuple[str, ...]]]],
) -> dict[str, str]:
    """Resolve the deepest recognized VirBot taxon onto one current ICTV path."""
    taxonomy = unclassified_taxonomy("virus")
    tokens = [
        clean_taxon(token)
        for token in clean_missing(raw_lineage).split(";")
        if clean_taxon(token)
    ]
    if not tokens:
        return taxonomy
    if tokens[0].casefold() != "viruses":
        raise ValueError(f"VirBot lineage does not begin with Viruses: {raw_lineage}")

    matches: list[tuple[str, int, tuple[str, ...]]] = []
    for token in tokens[1:]:
        for rank_index, current_lineage in ictv_index.get(token.casefold(), []):
            matches.append((token, rank_index, current_lineage))

    if not matches:
        return taxonomy

    deepest_rank = max(rank_index for _, rank_index, _ in matches)
    deepest_tokens = {
        token.casefold() for token, rank_index, _ in matches if rank_index == deepest_rank
    }
    if len(deepest_tokens) != 1:
        raise ValueError(
            "VirBot lineage maps to multiple ICTV taxa at the same deepest rank: "
            f"{raw_lineage}"
        )

    candidate_lineages = {
        current_lineage
        for token, rank_index, current_lineage in matches
        if rank_index == deepest_rank and token.casefold() in deepest_tokens
    }
    for rank_index in range(deepest_rank + 1):
        values = {
            lineage[rank_index]
            for lineage in candidate_lineages
            if lineage[rank_index]
        }
        if len(values) == 1:
            column = TAXONOMY_COLUMNS[rank_index + 1]
            prefix = column.split("__", 1)[0]
            taxonomy[column] = f"{prefix}__{values.pop()}"

    return taxonomy


def main() -> None:
    args = parse_args()
    header_records = load_header_map(args.header_map, args.sample_id)
    fasta_lengths = load_fasta_lengths(args.virus_fasta)
    metadata = load_metadata(args.run_metadata, args.sample_id, args.input_type)
    score_rows = read_table(args.score_table, SCORE_REQUIRED, ",")
    ictv_index = load_ictv_taxonomy(args.ictv_csv)

    input_count = parse_integer(
        metadata["input_sequence_count"], "input sequence count", args.sample_id
    )
    positive_count = parse_integer(
        metadata["positive_sequence_count"], "positive sequence count", args.sample_id
    )
    if input_count != len(header_records):
        raise ValueError("VirBot metadata and header-map input counts disagree")
    if positive_count != len(score_rows):
        raise ValueError("VirBot metadata and score-table positive counts disagree")

    expected_status = "completed" if score_rows else "completed_no_virus_calls"
    if metadata["run_status"] != expected_status:
        raise ValueError(
            f"VirBot run status '{metadata['run_status']}' disagrees with "
            f"the number of calls"
        )

    evidence_rows: list[dict[str, str]] = []
    score_ids: set[str] = set()
    for row in score_rows:
        sequence_id = row["Contig_acc"].strip()
        if not sequence_id:
            raise ValueError("VirBot score table contains an empty Contig_acc")
        if sequence_id in score_ids:
            raise ValueError(f"Duplicate VirBot sequence ID: {sequence_id}")
        if sequence_id not in header_records:
            raise ValueError(f"VirBot sequence is absent from header map: {sequence_id}")
        score_ids.add(sequence_id)

        _, score_text = parse_score(row["RNA-viral_gene_content"], sequence_id)
        n_genes = parse_integer(
            row["Encoded_proteins_num"], "encoded protein count", sequence_id, minimum=1
        )
        length = parse_integer(
            header_records[sequence_id]["length"],
            "header-map length",
            sequence_id,
            minimum=1,
        )
        if fasta_lengths.get(sequence_id) != length:
            raise ValueError(
                f"VirBot FASTA and header-map lengths disagree for {sequence_id}"
            )

        output = {
            "sample_id": args.sample_id,
            "sequence_id": sequence_id,
            "parent_sequence_id": "",
            "record_type": header_records[sequence_id]["record_type"].strip(),
            "coordinates": "",
            "tool": "virbot",
            "classification": "virus",
            "score": score_text,
            "score_type": "virbot_rna_viral_gene_fraction",
            "length": str(length),
            "topology": "",
            "n_genes": str(n_genes),
            "n_hallmarks": "",
            "evidence_strength": "qualified",
            "strength_basis": "virbot_adaptive_protein_cutoffs_and_gene_fraction",
        }
        output.update(validated_taxonomy(row["Likely_taxa"], ictv_index))
        evidence_rows.append(output)

    if set(fasta_lengths) != score_ids:
        missing = sorted(score_ids.difference(fasta_lengths))
        unexpected = sorted(set(fasta_lengths).difference(score_ids))
        raise ValueError(
            "VirBot score/FASTA identifiers disagree; "
            f"missing={missing[:5]}, unexpected={unexpected[:5]}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(evidence_rows)

    print(
        f"Wrote {args.output}: {len(evidence_rows)} ICTV-validated VirBot calls"
    )


if __name__ == "__main__":
    main()
