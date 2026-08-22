#!/usr/bin/env python3

"""Convert TEsorter classifications into auxiliary viSUM evidence."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from evidence_schema import CORE_EVIDENCE_COLUMNS, unclassified_taxonomy


CLASSIFICATION_REQUIRED = {
    "#TE", "Order", "Superfamily", "Clade", "Complete", "Strand", "Domains"
}
DOMAIN_REQUIRED = {"#id", "length", "evalue", "coverge", "probability", "score"}
REGION_REQUIRED = {
    "sample_id", "input_type", "sequence_id", "parent_sequence_id",
    "record_type", "coordinates", "refined_length",
}
METADATA_REQUIRED = {
    "sample_id", "input_type", "te_classification_count", "domain_row_count",
    "run_status",
}

VIRAL_LIKE_LABELS = {"retrovirus", "pararetrovirus"}
AMBIGUOUS_LABELS = {"mixture", "maverick", "polinton"}

OUTPUT_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "evidence_strength",
    "strength_basis",
    "tesorter_evidence_category",
    "tesorter_order",
    "tesorter_superfamily",
    "tesorter_clade",
    "tesorter_complete",
    "tesorter_strand",
    "assignment_method",
    "domain_count",
    "domains",
    "maximum_domain_probability",
    "maximum_domain_score",
    "maximum_domain_coverage",
    "minimum_domain_evalue",
]

AUDIT_COLUMNS = OUTPUT_COLUMNS + ["classification_row_number"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert TEsorter output into auxiliary viSUM evidence."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--region-map", required=True, type=Path)
    parser.add_argument("--classifications", required=True, type=Path)
    parser.add_argument("--domains", required=True, type=Path)
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--output-evidence", required=True, type=Path)
    parser.add_argument("--output-audit", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path, required: set[str], allow_empty: bool = False) -> list[dict[str, str]]:
    if path.stat().st_size == 0:
        if allow_empty:
            return []
        raise ValueError(f"Required TSV is empty: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return list(reader)


def parse_float(value: str, label: str, sequence_id: str) -> float:
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(f"Invalid {label} for {sequence_id}: {value}") from error


def format_number(value: float) -> str:
    return f"{value:.12g}"


def domain_sequence_id(domain_id: str) -> str:
    sequence_id = domain_id.split("|", 1)[0].strip()
    if not sequence_id:
        raise ValueError(f"Cannot extract sequence ID from TEsorter domain ID: {domain_id}")
    return sequence_id


def category_and_strength(row: dict[str, str], domain_count: int) -> tuple[str, str, str, str]:
    labels = {
        row["Order"].strip().lower(),
        row["Superfamily"].strip().lower(),
        row["Clade"].strip().lower(),
    }
    domains = row["Domains"].strip()
    direct_hmm = bool(domains and domains.lower() != "none" and domain_count > 0)

    if labels.intersection(VIRAL_LIKE_LABELS):
        category = "viral_like_mobile_element"
    elif labels.intersection(AMBIGUOUS_LABELS):
        category = "ambiguous_mobile_element"
    else:
        category = "retroelement"

    if not direct_hmm:
        return category, "weak", "tesorter_second_pass_similarity_transfer", "second_pass"
    if row["Complete"].strip().lower() == "yes":
        return category, "strong", "tesorter_complete_domain_architecture", "direct_hmm"
    return category, "qualified", "tesorter_direct_domain_support", "direct_hmm"


def main() -> None:
    args = parse_args()
    regions = read_tsv(args.region_map, REGION_REQUIRED)
    metadata_rows = read_tsv(args.run_metadata, METADATA_REQUIRED)
    if len(metadata_rows) != 1:
        raise ValueError("TEsorter run metadata must contain exactly one data row")
    metadata = metadata_rows[0]
    if metadata["sample_id"].strip() != args.sample_id:
        raise ValueError("TEsorter run metadata sample_id does not match --sample-id")
    if metadata["input_type"].strip() != args.input_type:
        raise ValueError("TEsorter run metadata input_type does not match --input-type")

    region_by_id: dict[str, dict[str, str]] = {}
    for row in regions:
        sequence_id = row["sequence_id"].strip()
        if not sequence_id or sequence_id in region_by_id:
            raise ValueError(f"Missing or duplicate region-map sequence ID: {sequence_id}")
        region_by_id[sequence_id] = row

    classifications = read_tsv(
        args.classifications, CLASSIFICATION_REQUIRED, allow_empty=True
    )
    domain_rows = read_tsv(args.domains, DOMAIN_REQUIRED, allow_empty=True)
    domains_by_sequence: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in domain_rows:
        sequence_id = domain_sequence_id(row["#id"])
        if sequence_id not in region_by_id:
            raise ValueError(f"TEsorter domain references an unknown refined sequence: {sequence_id}")
        domains_by_sequence[sequence_id].append(row)

    evidence_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for row_number, row in enumerate(classifications, start=2):
        sequence_id = row["#TE"].strip()
        if not sequence_id or sequence_id in seen:
            raise ValueError(f"Missing or duplicate TEsorter classification ID: {sequence_id}")
        seen.add(sequence_id)
        if sequence_id not in region_by_id:
            raise ValueError(f"TEsorter classification references an unknown refined sequence: {sequence_id}")

        region = region_by_id[sequence_id]
        sequence_domains = domains_by_sequence.get(sequence_id, [])
        reported_domains = row["Domains"].strip().lower()
        if reported_domains == "none" and sequence_domains:
            raise ValueError(
                f"TEsorter second-pass row unexpectedly has domain evidence: {sequence_id}"
            )
        if reported_domains not in {"", "none"} and not sequence_domains:
            raise ValueError(
                f"TEsorter direct-domain row has no matching domain records: {sequence_id}"
            )
        category, strength, basis, method = category_and_strength(
            row, len(sequence_domains)
        )
        probabilities = [
            parse_float(item["probability"], "domain probability", sequence_id)
            for item in sequence_domains
        ]
        scores = [
            parse_float(item["score"], "domain score", sequence_id)
            for item in sequence_domains
        ]
        coverages = [
            parse_float(item["coverge"], "domain coverage", sequence_id)
            for item in sequence_domains
        ]
        evalues = [
            parse_float(item["evalue"], "domain evalue", sequence_id)
            for item in sequence_domains
        ]

        output = {
            "sample_id": args.sample_id,
            "sequence_id": sequence_id,
            "parent_sequence_id": region["parent_sequence_id"].strip(),
            "record_type": region["record_type"].strip(),
            "coordinates": region["coordinates"].strip(),
            "tool": "tesorter",
            "classification": category,
            "score": format_number(max(probabilities)) if probabilities else "",
            "score_type": "maximum_tesorter_domain_probability" if probabilities else "",
            "length": region["refined_length"].strip(),
            "topology": "",
            "n_genes": "",
            "n_hallmarks": "",
            "evidence_strength": strength,
            "strength_basis": basis,
            "tesorter_evidence_category": category,
            "tesorter_order": row["Order"].strip(),
            "tesorter_superfamily": row["Superfamily"].strip(),
            "tesorter_clade": row["Clade"].strip(),
            "tesorter_complete": row["Complete"].strip(),
            "tesorter_strand": row["Strand"].strip(),
            "assignment_method": method,
            "domain_count": str(len(sequence_domains)),
            "domains": row["Domains"].strip(),
            "maximum_domain_probability": format_number(max(probabilities)) if probabilities else "",
            "maximum_domain_score": format_number(max(scores)) if scores else "",
            "maximum_domain_coverage": format_number(max(coverages)) if coverages else "",
            "minimum_domain_evalue": format_number(min(evalues)) if evalues else "",
            "classification_row_number": str(row_number),
        }
        output.update(unclassified_taxonomy(category))
        evidence_rows.append(output)

    orphan_domain_ids = sorted(set(domains_by_sequence).difference(seen))
    if orphan_domain_ids:
        raise ValueError(
            "TEsorter domain table contains sequences without classifications: "
            + ", ".join(orphan_domain_ids[:5])
        )

    expected_classifications = int(metadata["te_classification_count"])
    expected_domains = int(metadata["domain_row_count"])
    if len(evidence_rows) != expected_classifications:
        raise ValueError(
            "TEsorter classification count disagrees with run metadata: "
            f"observed={len(evidence_rows)} metadata={expected_classifications}"
        )
    if len(domain_rows) != expected_domains:
        raise ValueError(
            "TEsorter domain count disagrees with run metadata: "
            f"observed={len(domain_rows)} metadata={expected_domains}"
        )

    for path, columns in (
        (args.output_evidence, OUTPUT_COLUMNS),
        (args.output_audit, AUDIT_COLUMNS),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(evidence_rows)

    category_counts: dict[str, int] = defaultdict(int)
    strength_counts: dict[str, int] = defaultdict(int)
    for row in evidence_rows:
        category_counts[row["classification"]] += 1
        strength_counts[row["evidence_strength"]] += 1
    print(
        f"STANDARDIZED sample={args.sample_id} tool=tesorter rows={len(evidence_rows)} "
        f"categories={dict(sorted(category_counts.items()))} "
        f"strengths={dict(sorted(strength_counts.items()))}"
    )


if __name__ == "__main__":
    main()
