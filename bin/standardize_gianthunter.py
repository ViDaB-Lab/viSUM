#!/usr/bin/env python3

"""Convert GiantHunter NCLDV calls into sparse, ICTV-validated evidence."""

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from evidence_schema import CORE_EVIDENCE_COLUMNS, TAXONOMY_COLUMNS, unclassified_taxonomy


PREDICTION_COLUMNS = {
    "Accession",
    "Length",
    "GiantVirus",
    "PotentialLineage",
    "Score",
}
NO_REFERENCE_HIT_COLUMNS = {
    "Accession",
    "Length",
    "PotentialLineage",
    "Score",
    "Genus",
    "GenusCluster",
}
METADATA_COLUMNS = {
    "sample_id",
    "input_type",
    "input_sequence_count",
    "eligible_sequence_count",
    "giant_virus_call_count",
    "gianthunter_version",
    "minimum_length",
    "run_status",
}
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
GIANTHUNTER_TO_ICTV_RANK = {
    "clade": 0,
    "kingdom": 1,
    "phylum": 2,
    "class": 3,
    "order": 4,
    "family": 5,
    "genus": 6,
    "species": 7,
}
RUN_STATUSES = {
    "skipped_no_sequences_meeting_minimum_length",
    "completed_no_reference_protein_hits",
    "completed_no_giant_virus_calls",
    "completed_with_giant_virus_calls",
}
OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert GiantHunter NCLDV calls into standardized viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna",))
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--prediction-table", required=True, type=Path)
    parser.add_argument("--virus-fasta", required=True, type=Path)
    parser.add_argument("--gene-annotations", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--ictv-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def clean_missing(value: str | None) -> str:
    value = "" if value is None else value.strip()
    return "" if value.lower() in {"", "na", "nan", "none"} else value


def clean_taxon(value: str | None) -> str:
    return clean_missing(value).rstrip(".").strip()


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


def parse_score(value: str | None, sequence_id: str) -> str:
    cleaned = clean_missing(value)
    if not cleaned or cleaned == "-":
        raise ValueError(f"Missing GiantHunter score for {sequence_id}")
    numeric = float(cleaned)
    if not math.isfinite(numeric) or not 0.0 <= numeric <= 1.0:
        raise ValueError(
            f"GiantHunter score outside 0-1 for {sequence_id}: {cleaned}"
        )
    return cleaned


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


def load_metadata(
    path: Path,
    sample_id: str,
    input_type: str,
) -> dict[str, str]:
    rows = read_table(path, METADATA_COLUMNS, "\t")
    if len(rows) != 1:
        raise ValueError(
            f"Expected one GiantHunter metadata row in {path}; found {len(rows)}"
        )
    row = rows[0]
    if row["sample_id"] != sample_id or row["input_type"] != input_type:
        raise ValueError(
            "GiantHunter metadata does not match the requested sample and input type"
        )
    if input_type != "dna":
        raise ValueError("GiantHunter evidence is supported only for DNA inputs")
    if not clean_missing(row["gianthunter_version"]):
        raise ValueError("GiantHunter metadata contains an empty version")
    if row["run_status"] not in RUN_STATUSES:
        raise ValueError(f"Unrecognized GiantHunter run status: {row['run_status']}")

    input_count = parse_integer(
        row["input_sequence_count"], "input sequence count", sample_id
    )
    eligible_count = parse_integer(
        row["eligible_sequence_count"], "eligible sequence count", sample_id
    )
    call_count = parse_integer(
        row["giant_virus_call_count"], "giant-virus call count", sample_id
    )
    parse_integer(row["minimum_length"], "minimum length", sample_id, minimum=1)
    if eligible_count > input_count:
        raise ValueError("GiantHunter eligible count exceeds its input count")
    if call_count > eligible_count:
        raise ValueError("GiantHunter call count exceeds its eligible count")

    status = row["run_status"]
    if status == "skipped_no_sequences_meeting_minimum_length":
        if eligible_count != 0 or call_count != 0:
            raise ValueError("GiantHunter skipped status disagrees with its counts")
    elif status == "completed_with_giant_virus_calls":
        if eligible_count == 0 or call_count == 0:
            raise ValueError("GiantHunter positive status disagrees with its counts")
    elif eligible_count == 0 or call_count != 0:
        raise ValueError("GiantHunter no-call status disagrees with its counts")
    return row


def load_prediction_table(
    path: Path,
) -> tuple[str, list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"GiantHunter prediction table has no header: {path}")
        columns = set(reader.fieldnames)
        rows = list(reader)

    if PREDICTION_COLUMNS.issubset(columns):
        return "canonical", rows
    if NO_REFERENCE_HIT_COLUMNS.issubset(columns):
        return "no_reference_hits", rows
    raise ValueError(
        f"Unexpected GiantHunter prediction columns in {path}: {reader.fieldnames}"
    )


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


def load_gene_counts(path: Path) -> Counter[str]:
    if path.stat().st_size == 0:
        return Counter()
    rows = read_table(path, {"Genome"}, "\t")
    counts: Counter[str] = Counter()
    for row in rows:
        sequence_id = row["Genome"].strip()
        if not sequence_id:
            raise ValueError("GiantHunter gene annotations contain an empty Genome")
        counts[sequence_id] += 1
    return counts


def load_ictv_taxonomy(
    path: Path,
) -> dict[tuple[int, str], tuple[str, ...]]:
    rows = read_table(path, set(ICTV_RANKS), ",")
    occurrences: dict[tuple[int, str], set[tuple[str, ...]]] = defaultdict(set)
    for row in rows:
        lineage = tuple(clean_taxon(row[rank]) for rank in ICTV_RANKS)
        for rank_index, taxon in enumerate(lineage):
            if taxon:
                occurrences[(rank_index, taxon.casefold())].add(lineage)
    if not occurrences:
        raise ValueError(f"ICTV taxonomy table contains no taxon names: {path}")

    index: dict[tuple[int, str], tuple[str, ...]] = {}
    for key, lineages in occurrences.items():
        index[key] = tuple(
            next(iter(values)) if len(values) == 1 else ""
            for values in (
                {lineage[index] for lineage in lineages if lineage[index]}
                for index in range(len(ICTV_RANKS))
            )
        )
    return index


def ncldv_taxonomy() -> dict[str, str]:
    taxonomy = unclassified_taxonomy("virus")
    taxonomy["r__Realm"] = "r__Varidnaviria"
    taxonomy["k__Kingdom"] = "k__Bamfordvirae"
    taxonomy["p__Phylum"] = "p__Nucleocytoviricota"
    return taxonomy


def validated_taxonomy(
    raw_lineage: str,
    sequence_id: str,
    ictv_index: dict[tuple[int, str], tuple[str, ...]],
) -> dict[str, str]:
    taxonomy = ncldv_taxonomy()
    if raw_lineage.casefold() == "unclassified":
        return taxonomy

    tokens: dict[str, str] = {}
    for raw_token in raw_lineage.split(";"):
        token = raw_token.strip()
        if not token:
            continue
        if ":" not in token:
            raise ValueError(
                f"Malformed GiantHunter lineage token for {sequence_id}: {token}"
            )
        raw_rank, raw_taxon = token.split(":", 1)
        rank = raw_rank.strip().casefold()
        taxon = clean_taxon(raw_taxon)
        if not taxon:
            continue
        if rank in tokens and tokens[rank].casefold() != taxon.casefold():
            raise ValueError(
                f"Duplicate GiantHunter lineage rank for {sequence_id}: {raw_rank}"
            )
        tokens[rank] = taxon

    if tokens.get("superkingdom", "").casefold() != "viruses":
        raise ValueError(
            f"GiantHunter lineage does not begin with Viruses for {sequence_id}: "
            f"{raw_lineage}"
        )

    matches: list[tuple[int, tuple[str, ...]]] = []
    for raw_rank, rank_index in GIANTHUNTER_TO_ICTV_RANK.items():
        taxon = tokens.get(raw_rank)
        if not taxon:
            continue
        current_lineage = ictv_index.get((rank_index, taxon.casefold()))
        if current_lineage is not None:
            matches.append((rank_index, current_lineage))
    if not matches:
        return taxonomy

    deepest_rank = max(rank_index for rank_index, _ in matches)
    candidate_lineages = {
        lineage for rank_index, lineage in matches if rank_index == deepest_rank
    }
    if len(candidate_lineages) != 1:
        raise ValueError(
            f"Ambiguous current ICTV lineage for GiantHunter call {sequence_id}"
        )
    current_lineage = candidate_lineages.pop()

    expected_prefix = ("Varidnaviria", "Bamfordvirae", "Nucleocytoviricota")
    if tuple(current_lineage[:3]) != expected_prefix:
        raise ValueError(
            f"GiantHunter call {sequence_id} resolved outside Nucleocytoviricota"
        )

    for rank_index in range(deepest_rank + 1):
        taxon = current_lineage[rank_index]
        if taxon:
            column = TAXONOMY_COLUMNS[rank_index + 1]
            prefix = column.split("__", 1)[0]
            taxonomy[column] = f"{prefix}__{taxon}"
    return taxonomy


def main() -> None:
    args = parse_args()
    header_records = load_header_map(args.header_map, args.sample_id)
    metadata = load_metadata(args.run_metadata, args.sample_id, args.input_type)
    prediction_format, prediction_rows = load_prediction_table(args.prediction_table)
    fasta_lengths = load_fasta_lengths(args.virus_fasta)
    gene_counts = load_gene_counts(args.gene_annotations)
    ictv_index = load_ictv_taxonomy(args.ictv_csv)

    input_count = parse_integer(
        metadata["input_sequence_count"], "input sequence count", args.sample_id
    )
    eligible_count = parse_integer(
        metadata["eligible_sequence_count"], "eligible sequence count", args.sample_id
    )
    expected_call_count = parse_integer(
        metadata["giant_virus_call_count"], "giant-virus call count", args.sample_id
    )
    minimum_length = parse_integer(
        metadata["minimum_length"], "minimum length", args.sample_id, minimum=1
    )
    if input_count != len(header_records):
        raise ValueError("GiantHunter metadata and header-map input counts disagree")

    status = metadata["run_status"]
    if prediction_format == "no_reference_hits":
        if status != "completed_no_reference_protein_hits":
            raise ValueError(
                "GiantHunter alternate no-hit table disagrees with run metadata"
            )
        positive_rows: list[dict[str, str]] = []
    else:
        if status == "completed_no_reference_protein_hits":
            raise ValueError("GiantHunter no-hit metadata lacks its alternate table")
        if eligible_count > 0 and len(prediction_rows) != input_count:
            raise ValueError(
                "Completed GiantHunter prediction and input row counts disagree"
            )

        seen_ids: set[str] = set()
        positive_rows = []
        for row in prediction_rows:
            sequence_id = row["Accession"].strip()
            if not sequence_id:
                raise ValueError("GiantHunter prediction contains an empty Accession")
            if sequence_id in seen_ids:
                raise ValueError(f"Duplicate GiantHunter sequence ID: {sequence_id}")
            if sequence_id not in header_records:
                raise ValueError(
                    f"GiantHunter sequence is absent from header map: {sequence_id}"
                )
            seen_ids.add(sequence_id)

            prediction_length = parse_integer(
                row["Length"], "prediction length", sequence_id, minimum=1
            )
            header_length = parse_integer(
                header_records[sequence_id]["length"],
                "header-map length",
                sequence_id,
                minimum=1,
            )
            if prediction_length != header_length:
                raise ValueError(
                    f"GiantHunter and header-map lengths disagree for {sequence_id}"
                )

            label = row["GiantVirus"].strip()
            if label == "GiantVirus":
                if prediction_length < minimum_length:
                    raise ValueError(
                        f"GiantHunter positive is shorter than its minimum for {sequence_id}"
                    )
                parse_score(row["Score"], sequence_id)
                lineage = row["PotentialLineage"].strip()
                if not lineage or lineage == "-":
                    raise ValueError(
                        f"GiantHunter positive lacks a lineage state for {sequence_id}"
                    )
                positive_rows.append(row)
            elif label == "Non-GiantVirus":
                if row["PotentialLineage"].strip() != "-" or row["Score"].strip() != "-":
                    raise ValueError(
                        f"GiantHunter negative retains positive evidence for {sequence_id}"
                    )
            else:
                raise ValueError(
                    f"Unrecognized GiantHunter label for {sequence_id}: {label}"
                )

    if len(positive_rows) != expected_call_count:
        raise ValueError("GiantHunter metadata and prediction call counts disagree")

    positive_ids = {row["Accession"].strip() for row in positive_rows}
    if set(fasta_lengths) != positive_ids:
        missing = sorted(positive_ids.difference(fasta_lengths))
        unexpected = sorted(set(fasta_lengths).difference(positive_ids))
        raise ValueError(
            "GiantHunter prediction/FASTA identifiers disagree; "
            f"missing={missing[:5]}, unexpected={unexpected[:5]}"
        )

    evidence_rows: list[dict[str, str]] = []
    for row in positive_rows:
        sequence_id = row["Accession"].strip()
        length = parse_integer(row["Length"], "prediction length", sequence_id, minimum=1)
        if fasta_lengths[sequence_id] != length:
            raise ValueError(
                f"GiantHunter prediction and FASTA lengths disagree for {sequence_id}"
            )
        n_genes = gene_counts.get(sequence_id, 0)
        if n_genes < 1:
            raise ValueError(
                f"GiantHunter positive lacks gene annotations: {sequence_id}"
            )

        raw_lineage = row["PotentialLineage"].strip()
        score_type = (
            "gianthunter_model_score"
            if raw_lineage.casefold() == "unclassified"
            else "gianthunter_weighted_lca_support"
        )
        output = {
            "sample_id": args.sample_id,
            "sequence_id": sequence_id,
            "parent_sequence_id": "",
            "record_type": header_records[sequence_id]["record_type"].strip(),
            "coordinates": "",
            "tool": "gianthunter",
            "classification": "virus",
            "score": parse_score(row["Score"], sequence_id),
            "score_type": score_type,
            "length": str(length),
            "topology": "",
            "n_genes": str(n_genes),
            "n_hallmarks": "",
        }
        output.update(
            validated_taxonomy(raw_lineage, sequence_id, ictv_index)
        )
        evidence_rows.append(output)

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
        f"Wrote {args.output}: {len(evidence_rows)} GiantHunter NCLDV calls"
    )


if __name__ == "__main__":
    main()
