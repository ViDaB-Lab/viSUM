#!/usr/bin/env python3

"""Standardize vConTACT3 taxonomy calls and preserve project group membership."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

from evidence_schema import CORE_EVIDENCE_COLUMNS, unclassified_taxonomy


RANKS = ("realm", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus")
FORMAL_RANKS = ("realm", "kingdom", "phylum", "class", "order", "family", "genus")
RANK_COLUMNS = {
    "realm": "r__Realm",
    "kingdom": "k__Kingdom",
    "phylum": "p__Phylum",
    "class": "c__Class",
    "order": "o__Order",
    "family": "f__Family",
    "genus": "g__Genus",
}
SYNTHETIC_PREFIXES = ("novel_", "unplaced_")
MISSING_VALUES = {"", "na", "nan", "none", "null", "np.nan", "-"}
DOMAIN_PATTERN = re.compile(
    r"\.vcontact3_(prokaryotes|eukaryotes)_final_assignments\.csv$"
)

REGION_REQUIRED = {
    "sample_id",
    "input_type",
    "sequence_id",
    "parent_sequence_id",
    "record_type",
    "coordinates",
    "refined_length",
}
METADATA_REQUIRED = {
    "sample_id",
    "input_type",
    "database_domain",
    "refined_sequence_count",
    "database_version",
    "run_status",
}
ASSIGNMENT_REQUIRED = {
    "Genome",
    "GenomeName",
    "Proteins",
    "Reference",
    "Size_Kb",
    "realm_prediction",
    "host_domain",
    "vog_host_domain",
}
for _rank in RANKS[1:]:
    ASSIGNMENT_REQUIRED.update(
        {
            f"{_rank}_prediction",
            f"{_rank}_network_support",
            f"{_rank}_evidence",
        }
    )

EVIDENCE_COLUMNS = CORE_EVIDENCE_COLUMNS + [
    "classification_rank",
    "vcontact3_assignment_method",
    "vcontact3_database_domain",
    "vcontact3_database_version",
    "vcontact3_host_domain",
    "vcontact3_vog_host_domain",
]

GROUP_COLUMNS = [
    "sample_id",
    "input_type",
    "sequence_id",
    "parent_sequence_id",
    "record_type",
    "coordinates",
    "length",
    "database_domain",
    "database_version",
    "candidate_assignment_status",
    "raw_reference_flag",
    "proteins",
    "size_kb",
    "has_named_realm",
    "has_reference_taxonomy",
    "deepest_reference_rank",
    "novel_group_ranks",
    "unplaced_group_ranks",
    "standardization_decision",
]
for _rank in RANKS:
    GROUP_COLUMNS.append(f"{_rank}_prediction")
    if _rank != "realm":
        GROUP_COLUMNS.extend(
            [
                f"{_rank}_evidence",
                f"{_rank}_network_support",
                f"{_rank}_vog_support",
            ]
        )
GROUP_COLUMNS.extend(["host_domain", "vog_host_domain"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert vConTACT3 assignments into sparse formal evidence and a "
            "complete harmonizer-facing project-group table."
        )
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--region-map", required=True, type=Path)
    parser.add_argument(
        "--assignments", required=True, action="append", type=Path
    )
    parser.add_argument("--run-metadata", required=True, type=Path)
    parser.add_argument("--output-evidence", required=True, type=Path)
    parser.add_argument("--output-groups", required=True, type=Path)
    return parser.parse_args()


def clean(value: str | None) -> str:
    cleaned = "" if value is None else str(value).strip()
    return "" if cleaned.lower() in MISSING_VALUES else cleaned


def read_table(
    path: Path, required: set[str], delimiter: str
) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Table has no header: {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return list(reader)


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_nonnegative_integer(value: str, label: str) -> int:
    cleaned = clean(value)
    if not cleaned:
        raise ValueError(f"Missing {label}")
    try:
        parsed = int(cleaned)
    except ValueError as error:
        raise ValueError(f"{label} is not an integer: {value}") from error
    if parsed < 0:
        raise ValueError(f"{label} must be nonnegative: {value}")
    return parsed


def validate_optional_score(value: str, label: str) -> str:
    cleaned = clean(value)
    if not cleaned:
        return ""
    try:
        score = float(cleaned)
    except ValueError as error:
        raise ValueError(f"{label} is not numeric: {value}") from error
    if not math.isfinite(score) or score < 0:
        raise ValueError(f"{label} must be finite and nonnegative: {value}")
    return cleaned


def load_regions(
    path: Path, sample_id: str, input_type: str
) -> dict[str, dict[str, str]]:
    rows = read_table(path, REGION_REQUIRED, "\t")
    regions: dict[str, dict[str, str]] = {}
    for row in rows:
        sequence_id = clean(row["sequence_id"])
        if not sequence_id:
            raise ValueError("Region map contains an empty sequence_id")
        if sequence_id in regions:
            raise ValueError(f"Duplicate region-map sequence: {sequence_id}")
        if clean(row["sample_id"]) != sample_id:
            raise ValueError(f"Region-map sample mismatch for {sequence_id}")
        if clean(row["input_type"]) != input_type:
            raise ValueError(f"Region-map input-type mismatch for {sequence_id}")
        length = parse_nonnegative_integer(
            row["refined_length"], f"refined length for {sequence_id}"
        )
        if length == 0:
            raise ValueError(f"Refined sequence has zero length: {sequence_id}")
        regions[sequence_id] = row
    return regions


def load_metadata(
    path: Path, sample_id: str, input_type: str, expected_count: int
) -> dict[str, dict[str, str]]:
    rows = read_table(path, METADATA_REQUIRED, "\t")
    metadata: dict[str, dict[str, str]] = {}
    for row in rows:
        domain = clean(row["database_domain"])
        if domain not in {"prokaryotes", "eukaryotes"}:
            raise ValueError(f"Unsupported metadata database domain: {domain}")
        if domain in metadata:
            raise ValueError(f"Duplicate metadata row for database domain: {domain}")
        if clean(row["sample_id"]) != sample_id:
            raise ValueError(f"Run-metadata sample mismatch for {domain}")
        if clean(row["input_type"]) != input_type:
            raise ValueError(f"Run-metadata input-type mismatch for {domain}")
        observed_count = parse_nonnegative_integer(
            row["refined_sequence_count"], f"refined sequence count for {domain}"
        )
        if observed_count != expected_count:
            raise ValueError(
                f"Run-metadata refined count mismatch for {domain}: "
                f"expected {expected_count}, observed {observed_count}"
            )
        run_status = clean(row["run_status"])
        valid_status = (
            run_status == "completed"
            if expected_count
            else run_status == "skipped_no_refined_candidates"
        )
        if not valid_status:
            raise ValueError(f"Unexpected vConTACT3 run status for {domain}: {run_status}")
        metadata[domain] = row
    return metadata


def assignment_domain(path: Path) -> str:
    match = DOMAIN_PATTERN.search(path.name)
    if match is None:
        raise ValueError(f"Cannot determine vConTACT3 database domain from: {path.name}")
    return match.group(1)


def is_synthetic_taxon(value: str) -> bool:
    taxon = clean(value)
    lowered = taxon.lower()
    return (
        not taxon
        or lowered in {"default", "singleton", "unclassified"}
        or lowered.startswith(SYNTHETIC_PREFIXES)
        or "|" in taxon
    )


def is_named_taxon(value: str) -> bool:
    return not is_synthetic_taxon(value)


def candidate_identifier(
    row: dict[str, str], candidate_ids: set[str], path: Path
) -> str:
    matches = {
        clean(row.get(column))
        for column in ("Genome", "GenomeName")
        if clean(row.get(column)) in candidate_ids
    }
    if len(matches) > 1:
        raise ValueError(
            f"Assignment row matches multiple candidate IDs in {path}: {sorted(matches)}"
        )
    return next(iter(matches), "")


def load_candidate_assignments(
    path: Path, regions: dict[str, dict[str, str]]
) -> dict[str, dict[str, str]]:
    rows = read_table(path, ASSIGNMENT_REQUIRED, ",")
    candidates: dict[str, dict[str, str]] = {}
    candidate_ids = set(regions)
    for row in rows:
        sequence_id = candidate_identifier(row, candidate_ids, path)
        if not sequence_id:
            continue
        if sequence_id in candidates:
            raise ValueError(f"Duplicate candidate assignment in {path}: {sequence_id}")
        reference_flag = clean(row.get("Reference")).lower()
        if reference_flag in {"true", "1"}:
            raise ValueError(
                f"Candidate ID collides with a vConTACT3 reference row in {path}: "
                f"{sequence_id}"
            )
        candidates[sequence_id] = row
    missing = candidate_ids.difference(candidates)
    if missing:
        preview = ", ".join(sorted(missing)[:5])
        raise ValueError(
            f"vConTACT3 output {path} is missing {len(missing)} refined candidates: "
            f"{preview}"
        )
    return candidates


def formal_taxonomy(
    row: dict[str, str]
) -> tuple[dict[str, str], str, str, str]:
    accepted: dict[str, str] = {}
    realm = clean(row.get("realm_prediction"))
    if is_named_taxon(realm):
        accepted["realm"] = realm

    deepest_formal_reference = ""
    deepest_reference = ""
    for rank in RANKS[1:]:
        prediction = clean(row.get(f"{rank}_prediction"))
        evidence = clean(row.get(f"{rank}_evidence")).lower()
        if evidence == "reference" and is_named_taxon(prediction):
            deepest_reference = rank
            if rank != "subfamily":
                accepted[rank] = prediction
                deepest_formal_reference = rank

    taxonomy = unclassified_taxonomy("virus")
    for rank, taxon in accepted.items():
        if rank in RANK_COLUMNS:
            column = RANK_COLUMNS[rank]
            taxonomy[column] = f"{column[:3]}{taxon}"

    deepest = ""
    for rank in FORMAL_RANKS:
        if rank in accepted:
            deepest = rank
    return taxonomy, deepest, deepest_formal_reference, deepest_reference


def assignment_status(row: dict[str, str]) -> str:
    reference_flag = clean(row.get("Reference")).lower()
    realm = clean(row.get("realm_prediction")).lower()
    if realm == "singleton":
        return "singleton"
    if reference_flag in {"false", "0"}:
        return "clustered"
    if any(clean(row.get(f"{rank}_prediction")) for rank in RANKS):
        return "assigned_without_reference_flag"
    return "unassigned"


def group_decision(
    row: dict[str, str], named_realm: bool, deepest_reference: str
) -> str:
    if deepest_reference:
        return "reference_supported_taxonomy"
    if named_realm:
        return "realm_only_taxonomy"
    if any(
        clean(row.get(f"{rank}_prediction")).lower().startswith("novel_")
        for rank in RANKS
    ):
        return "novel_group_only"
    if clean(row.get("realm_prediction")).lower() == "singleton":
        return "singleton"
    return "no_qualified_taxonomy"


def build_group_row(
    sample_id: str,
    input_type: str,
    sequence_id: str,
    region: dict[str, str],
    row: dict[str, str],
    domain: str,
    database_version: str,
    deepest_reference: str,
) -> dict[str, str]:
    named_realm = is_named_taxon(clean(row.get("realm_prediction")))
    novel_ranks = [
        rank
        for rank in RANKS
        if clean(row.get(f"{rank}_prediction")).lower().startswith("novel_")
    ]
    unplaced_ranks = [
        rank
        for rank in RANKS
        if clean(row.get(f"{rank}_prediction")).lower().startswith("unplaced_")
    ]
    output = {
        "sample_id": sample_id,
        "input_type": input_type,
        "sequence_id": sequence_id,
        "parent_sequence_id": clean(region["parent_sequence_id"]),
        "record_type": clean(region["record_type"]),
        "coordinates": clean(region["coordinates"]),
        "length": clean(region["refined_length"]),
        "database_domain": domain,
        "database_version": database_version,
        "candidate_assignment_status": assignment_status(row),
        "raw_reference_flag": clean(row.get("Reference")),
        "proteins": clean(row.get("Proteins")),
        "size_kb": clean(row.get("Size_Kb")),
        "has_named_realm": str(named_realm).lower(),
        "has_reference_taxonomy": str(bool(deepest_reference)).lower(),
        "deepest_reference_rank": deepest_reference,
        "novel_group_ranks": ",".join(novel_ranks),
        "unplaced_group_ranks": ",".join(unplaced_ranks),
        "standardization_decision": group_decision(
            row, named_realm, deepest_reference
        ),
        "host_domain": clean(row.get("host_domain")),
        "vog_host_domain": clean(row.get("vog_host_domain")),
    }
    for rank in RANKS:
        output[f"{rank}_prediction"] = clean(row.get(f"{rank}_prediction"))
        if rank != "realm":
            output[f"{rank}_evidence"] = clean(row.get(f"{rank}_evidence"))
            output[f"{rank}_network_support"] = validate_optional_score(
                row.get(f"{rank}_network_support", ""),
                f"{rank} network support for {sequence_id} ({domain})",
            )
            output[f"{rank}_vog_support"] = validate_optional_score(
                row.get(f"{rank}_vog_support", ""),
                f"{rank} VOG support for {sequence_id} ({domain})",
            )
    return output


def build_evidence_row(
    sample_id: str,
    sequence_id: str,
    region: dict[str, str],
    row: dict[str, str],
    domain: str,
    database_version: str,
    taxonomy: dict[str, str],
    deepest: str,
    deepest_reference: str,
) -> dict[str, str]:
    score = ""
    score_type = ""
    if deepest_reference:
        score = validate_optional_score(
            row.get(f"{deepest_reference}_network_support", ""),
            f"{deepest_reference} network support for {sequence_id} ({domain})",
        )
        if score:
            score_type = "vcontact3_network_support"

    output = {
        "sample_id": sample_id,
        "sequence_id": sequence_id,
        "parent_sequence_id": clean(region["parent_sequence_id"]),
        "record_type": clean(region["record_type"]),
        "coordinates": clean(region["coordinates"]),
        "tool": "vcontact3",
        "classification": "virus",
        "score": score,
        "score_type": score_type,
        "length": clean(region["refined_length"]),
        "topology": "",
        "n_genes": clean(row.get("Proteins")),
        "n_hallmarks": "",
        "classification_rank": deepest,
        "vcontact3_assignment_method": (
            "reference_supported" if deepest_reference else "realm_only"
        ),
        "vcontact3_database_domain": domain,
        "vcontact3_database_version": database_version,
        "vcontact3_host_domain": clean(row.get("host_domain")),
        "vcontact3_vog_host_domain": clean(row.get("vog_host_domain")),
    }
    output.update(taxonomy)
    return output


def main() -> None:
    args = parse_args()
    regions = load_regions(args.region_map, args.sample_id, args.input_type)
    metadata = load_metadata(
        args.run_metadata, args.sample_id, args.input_type, len(regions)
    )

    assignment_paths: dict[str, Path] = {}
    for path in args.assignments:
        domain = assignment_domain(path)
        if domain in assignment_paths:
            raise ValueError(f"Duplicate assignment input for database domain: {domain}")
        assignment_paths[domain] = path
    if set(assignment_paths) != set(metadata):
        raise ValueError(
            "Assignment/metadata domain mismatch: "
            f"assignments={sorted(assignment_paths)}, metadata={sorted(metadata)}"
        )

    evidence_rows: list[dict[str, str]] = []
    group_rows: list[dict[str, str]] = []
    if regions:
        for domain in sorted(assignment_paths):
            candidates = load_candidate_assignments(assignment_paths[domain], regions)
            database_version = clean(metadata[domain]["database_version"])
            for sequence_id, region in regions.items():
                row = candidates[sequence_id]
                (
                    taxonomy,
                    deepest,
                    deepest_formal_reference,
                    deepest_reference,
                ) = formal_taxonomy(row)
                group_rows.append(
                    build_group_row(
                        args.sample_id,
                        args.input_type,
                        sequence_id,
                        region,
                        row,
                        domain,
                        database_version,
                        deepest_reference,
                    )
                )
                if deepest:
                    evidence_rows.append(
                        build_evidence_row(
                            args.sample_id,
                            sequence_id,
                            region,
                            row,
                            domain,
                            database_version,
                            taxonomy,
                            deepest,
                            deepest_formal_reference,
                        )
                    )

    write_tsv(args.output_evidence, EVIDENCE_COLUMNS, evidence_rows)
    write_tsv(args.output_groups, GROUP_COLUMNS, group_rows)
    print(
        f"Wrote {args.output_evidence}: {len(evidence_rows)} qualified domain calls; "
        f"{args.output_groups}: {len(group_rows)} candidate-domain memberships"
    )


if __name__ == "__main__":
    main()
