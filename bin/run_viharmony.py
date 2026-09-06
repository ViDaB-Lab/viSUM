#!/usr/bin/env python3

"""Combine viSUM evidence into final sequences, decisions, and taxonomy calls."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO


RANKS = ("domain", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species")
RANK_COLUMNS = {
    "domain": "d__Domain", "realm": "r__Realm", "kingdom": "k__Kingdom",
    "phylum": "p__Phylum", "class": "c__Class", "order": "o__Order",
    "family": "f__Family", "genus": "g__Genus", "species": "s__Species",
}
PREFIXES = {rank: column[:3] for rank, column in RANK_COLUMNS.items()}
UNCLASSIFIED = {rank: f"{PREFIXES[rank]}unclassified" for rank in RANKS}

ORIGIN_TOOLS = {
    "genomad", "virsorter2", "cenotetaker3", "deep6", "deepmicroclass2",
    "virbot", "gianthunter", "vicat", "checkv",
}
TAXONOMY_TOOLS = {
    "genomad", "cenotetaker3", "virbot", "gianthunter", "vicat", "vitap", "vcontact3",
}
METHOD_FAMILY = {
    "genomad": "marker_gene", "cenotetaker3": "hallmark_gene",
    "virbot": "protein_model", "gianthunter": "protein_model",
    "vicat": "protein_homology", "vitap": "network_homology",
    "vcontact3": "gene_sharing_network",
}
TIER_SCORE = {"provisional": 1, "moderate": 2, "high": 3}

METADATA_COLUMNS = [
    "original_contig_name", "normalized_name", "final_sequence_id", "record_type",
    "provirus_coordinates", "original_length", "refined_length", "discovery_status",
    "viral_decision",
    "viral_confidence", "strong_tools", "qualified_tools", "cellular_conflict_tools",
    "checkv_complete_viral_contig_support",
    "plasmid_conflict_tools", "strict_taxonomy", "strict_taxonomy_rank",
    "analysis_taxonomy", "analysis_taxonomy_rank", "taxonomy_confidence",
    "taxonomy_supporting_tools", "taxonomy_conflict", "vcontact3_groups",
    "vcontact3_group_status", "sequence_interpretation", "tesorter_status",
    "tesorter_evidence_strength", "tesorter_categories", "tesorter_orders",
    "tesorter_superfamilies", "tesorter_assignment_methods",
    "vicat_origin_pattern", "vicat_provirus_status", "vicat_evidence_scope",
    "vicat_viral_supported_loci", "vicat_cellular_supported_loci",
    "vicat_nonviral_supported_classes", "vicat_dominant_nonviral_class",
]

STRENGTH_ORDER = {"": 0, "weak": 1, "qualified": 2, "strong": 3}

PRIMARY_VIRAL_DECISIONS = {"retained_viral", "retained_provirus"}
PROVISIONAL_VIRAL_DECISIONS = {"provisional_viral", "provisional_provirus"}


def parse_optional_interval(value: str) -> tuple[int, int] | None:
    coordinate = value.strip()
    if not coordinate or coordinate.upper() == "NA" or coordinate.count("-") != 1:
        return None
    left, right = coordinate.split("-", 1)
    try:
        start, end = int(left), int(right)
    except ValueError:
        return None
    return (start, end) if start >= 1 and end >= start else None


def intervals_overlap(
    left: tuple[int, int] | None, right: tuple[int, int] | None
) -> bool:
    return bool(
        left and right and left[0] <= right[1] and right[0] <= left[1]
    )


def evidence_applies_to_final(
    exact_id: str,
    parent_id: str,
    row_record_type: str,
    row_interval: tuple[int, int] | None,
    final_id: str,
    final_parent: str,
    final_record_type: str,
    final_interval: tuple[int, int] | None,
    final_length: int,
    known_final_ids: set[str],
) -> bool:
    """Determine whether one evidence row applies to one final sequence."""
    if exact_id in known_final_ids:
        return exact_id == final_id
    region_specific = (
        row_record_type in {"provirus", "viral_region"}
        and row_interval is not None
    )
    if region_specific:
        if parent_id != final_parent:
            return False
        if final_record_type == "input_contig":
            # A selected exact full-span call is represented downstream by its
            # parent ID. Carry that evidence across the identity collapse, but
            # never attach a partial/unselected region call to the full contig.
            return row_interval == (1, final_length)
        return intervals_overlap(row_interval, final_interval)
    return parent_id == final_parent or exact_id == final_parent


def adjudicate_viral_decision(
    discovery_status: str,
    record_type: str,
    strong_tools: set[str],
    qualified_tools: set[str],
    cellular_tools: list[str],
    plasmid_tools: list[str],
    sequence_interpretation: str,
) -> str:
    """Convert gathered evidence into a conservative final disposition.

    Refinement establishes the sequence interval to evaluate; it does not prove
    viral origin.  Refined regions therefore pass through the same origin and
    mobile-element conflict checks as intact contigs.  Only a refined region
    that survives those checks receives the ``retained_provirus`` disposition.
    """
    refined_region = record_type in {"provirus", "viral_region"}
    if sequence_interpretation == "likely_retroelement":
        return "likely_retroelement"
    if sequence_interpretation in {
        "viral_retroelement_conflict", "ambiguous_mobile_element"
    }:
        return "ambiguous_review"

    conflicts = bool(cellular_tools or plasmid_tools)
    if conflicts:
        if discovery_status == "ambiguous" and not strong_tools and len(qualified_tools) <= 1:
            return "likely_nonviral"
        return "ambiguous_review"
    if strong_tools or len(qualified_tools) >= 2:
        return "retained_provirus" if refined_region else "retained_viral"
    if qualified_tools:
        return "provisional_provirus" if refined_region else "provisional_viral"
    return "ambiguous_review"


def checkv_supports_complete_viral_contig(
    rows: list[dict[str, str]], record_type: str
) -> bool:
    """Recognize CheckV's strongest intact-contig viral-origin pattern."""
    if record_type != "input_contig":
        return False
    for row in rows:
        if row.get("tool", "").strip().lower() != "checkv":
            continue
        if row.get("classification", "").strip().lower() != "virus":
            continue
        if row.get("checkv_quality", "").strip().lower() not in {
            "complete", "high-quality"
        }:
            continue
        if row.get("provirus", "").strip().lower() != "no":
            continue
        try:
            viral_genes = int(float(row.get("viral_genes", "0") or 0))
            host_genes = int(float(row.get("host_genes", "0") or 0))
        except ValueError:
            continue
        if viral_genes > 0 and host_genes == 0:
            return True
    return False


def summarize_tesorter(
    rows: list[dict[str, str]], has_qualified_viral_evidence: bool,
    has_high_viral_consensus: bool = False,
) -> dict[str, str]:
    """Summarize TEsorter as auxiliary mobile-element evidence.

    TEsorter never contributes an ICTV taxonomy vote. Its result changes only
    the reported biological interpretation and review status.
    """
    tesorter_rows = [
        row for row in rows if row.get("tool", "").strip().lower() == "tesorter"
    ]
    if not tesorter_rows:
        return {
            "sequence_interpretation": "viral_candidate",
            "tesorter_status": "no_tesorter_evidence",
            "tesorter_evidence_strength": "NA",
            "tesorter_categories": "NA",
            "tesorter_orders": "NA",
            "tesorter_superfamilies": "NA",
            "tesorter_assignment_methods": "NA",
        }

    categories = sorted({
        row.get("tesorter_evidence_category", "").strip()
        for row in tesorter_rows
        if row.get("tesorter_evidence_category", "").strip()
    })
    strengths = [
        row.get("evidence_strength", "").strip().lower() for row in tesorter_rows
    ]
    maximum_strength = max(strengths, key=lambda value: STRENGTH_ORDER.get(value, 0))
    retroelement_strengths = [
        row.get("evidence_strength", "").strip().lower()
        for row in tesorter_rows
        if row.get("tesorter_evidence_category", "").strip() == "retroelement"
    ]
    retroelement_strength = (
        max(retroelement_strengths, key=lambda value: STRENGTH_ORDER.get(value, 0))
        if retroelement_strengths else ""
    )
    has_supported_retroelement = STRENGTH_ORDER.get(retroelement_strength, 0) >= 2
    has_weak_retroelement = retroelement_strength == "weak"
    has_viral_like = "viral_like_mobile_element" in categories
    has_ambiguous = "ambiguous_mobile_element" in categories

    # An explicit Retrovirus/Pararetrovirus classification is biologically
    # compatible with viral origin.  When an independent viSUM caller also
    # supports viral origin, do not let a generic retroelement annotation from
    # another TEsorter window take precedence over the explicit viral label.
    if has_viral_like and has_qualified_viral_evidence:
        interpretation = "retrovirus_compatible"
        status = (
            "viral_like_mobile_element_with_mixed_retroelement_evidence"
            if has_supported_retroelement
            else "viral_like_mobile_element_with_viral_support"
        )
    elif (
        (has_supported_retroelement or has_ambiguous)
        and has_high_viral_consensus
    ):
        interpretation = "viral_with_mobile_element_features"
        status = "mobile_element_annotation_overridden_by_high_viral_consensus"
    elif has_supported_retroelement and has_qualified_viral_evidence:
        interpretation = "viral_retroelement_conflict"
        status = f"{retroelement_strength}_retroelement_conflict"
    elif has_supported_retroelement:
        interpretation = "likely_retroelement"
        status = f"{retroelement_strength}_retroelement_evidence"
    elif has_viral_like:
        interpretation = "viral_like_mobile_element"
        status = "viral_like_mobile_element"
    elif has_ambiguous:
        interpretation = "ambiguous_mobile_element"
        status = "ambiguous_mobile_element_review"
    elif has_weak_retroelement and has_qualified_viral_evidence:
        interpretation = "viral_candidate_with_weak_retroelement_signal"
        status = "weak_retroelement_annotation"
    elif has_weak_retroelement:
        interpretation = "weak_retroelement_signal"
        status = "weak_retroelement_annotation"
    else:
        interpretation = "viral_candidate"
        status = "tesorter_unresolved"

    def joined(column: str) -> str:
        values = sorted({
            row.get(column, "").strip() for row in tesorter_rows
            if row.get(column, "").strip() and row.get(column, "").strip().lower() != "unknown"
        })
        return ",".join(values) or "NA"

    return {
        "sequence_interpretation": interpretation,
        "tesorter_status": status,
        "tesorter_evidence_strength": maximum_strength or "NA",
        "tesorter_categories": ",".join(categories) or "NA",
        "tesorter_orders": joined("tesorter_order"),
        "tesorter_superfamilies": joined("tesorter_superfamily"),
        "tesorter_assignment_methods": joined("assignment_method"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the viHARMONY decision engine.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", choices=("dna", "rna"), required=True)
    parser.add_argument("--normalized-fasta", required=True, type=Path)
    parser.add_argument("--header-map", required=True, type=Path)
    parser.add_argument("--discovery-gate", required=True, type=Path)
    parser.add_argument("--refined-fasta", required=True, type=Path)
    parser.add_argument("--region-map", required=True, type=Path)
    parser.add_argument("--ictv-msl", required=True, type=Path)
    parser.add_argument("--evidence", nargs="*", default=[], type=Path)
    parser.add_argument("--vcontact3-groups", nargs="*", default=[], type=Path)
    parser.add_argument("--vcontact3-min-taxonomy-length", type=int, default=1000)
    parser.add_argument("--audit-mode", choices=("none", "compact", "full"), default="compact")
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Table has no header: {path}")
        return list(reader.fieldnames), list(reader)


def require(columns: Iterable[str], required: set[str], path: Path) -> None:
    missing = required.difference(columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(parts)
                identifier = line[1:].split(maxsplit=1)[0]
                if not identifier or identifier in records:
                    raise ValueError(f"Invalid or duplicate FASTA identifier at {path}:{line_number}")
                parts = []
            elif identifier is None:
                raise ValueError(f"Sequence precedes FASTA header at {path}:{line_number}")
            else:
                parts.append("".join(line.split()))
    if identifier is not None:
        records[identifier] = "".join(parts)
    if any(not sequence for sequence in records.values()):
        raise ValueError(f"Empty FASTA sequence in {path}")
    return records


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_taxon(value: str, rank: str) -> str:
    value = (value or "").strip()
    prefix = PREFIXES[rank]
    if value.startswith(prefix):
        value = value[len(prefix):]
    if not value or value.lower() in {"unclassified", "na", "nan", "none", "-"}:
        return ""
    return value


def load_msl(path: Path) -> tuple[set[tuple[str, ...]], str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"ICTV MSL has no header: {path}")
        normalized = {name.strip().lower(): name for name in reader.fieldnames}
        required = [rank for rank in RANKS if rank != "domain"]
        missing = [rank for rank in required if rank not in normalized]
        if missing:
            raise ValueError(f"ICTV MSL lacks rank columns: {missing}")
        lineages = {
            tuple((row.get(normalized[rank]) or "").strip() for rank in required)
            for row in reader
        }
    if not lineages:
        raise ValueError(f"ICTV MSL contains no taxonomy rows: {path}")
    release = path.stem
    return lineages, release


def lineage_is_valid(lineage: dict[str, str], msl: set[tuple[str, ...]]) -> bool:
    named = [(rank, lineage.get(rank, "")) for rank in RANKS[1:] if lineage.get(rank)]
    if not named:
        return False
    rank_positions = {rank: index for index, rank in enumerate(RANKS[1:])}
    return any(all(msl_row[rank_positions[rank]] == taxon for rank, taxon in named) for msl_row in msl)


def vote_tier(row: dict[str, str]) -> str:
    tool = row.get("tool", "").strip().lower()
    if tool == "vitap":
        confidence = row.get("vitap_confidence_level", "").lower()
        if "high" in confidence:
            return "high"
        if "medium" in confidence or "uniref" in confidence:
            return "moderate"
        return "provisional"
    if tool == "vcontact3":
        return "moderate" if row.get("vcontact3_assignment_method") == "reference_supported" else "provisional"
    if tool == "vicat":
        try:
            loci = int(float(row.get("taxonomy_supporting_loci", "0") or 0))
        except ValueError:
            loci = 0
        return "moderate" if loci >= 2 and row.get("taxonomy_conflict", "").lower() not in {"true", "1", "yes"} else "provisional"
    if tool in {"cenotetaker3", "virbot"}:
        try:
            support = int(float(row.get("n_hallmarks", "0") or row.get("n_genes", "0") or 0))
        except ValueError:
            support = 0
        return "moderate" if support >= 2 else "provisional"
    if tool == "gianthunter":
        return "moderate" if "lca" in row.get("score_type", "").lower() else "provisional"
    if tool == "genomad":
        return "moderate"
    return "provisional"


def root_for_evidence(row: dict[str, str], final_to_parent: dict[str, str]) -> str:
    sequence_id = row.get("sequence_id", "").strip()
    parent = row.get("parent_sequence_id", "").strip()
    if sequence_id in final_to_parent:
        return sequence_id
    target_parent = parent or sequence_id
    matches = [final_id for final_id, final_parent in final_to_parent.items() if final_parent == target_parent]
    return matches[0] if len(matches) == 1 else target_parent


def select_vicat_scope(
    rows: list[dict[str, str]], final_id: str, parent_id: str
) -> tuple[list[dict[str, str]], list[dict[str, str]], bool]:
    """Prefer region-scoped viCAT evidence without erasing parent provenance."""
    region_rows = [
        row for row in rows
        if row.get("tool", "").strip().lower() == "vicat"
        and row.get("evidence_scope", "").strip() == "refined_region"
        and row.get("sequence_id", "").strip() == final_id
    ]
    suppress_parent = final_id != parent_id or bool(region_rows)
    selected = [
        row for row in rows
        if not (
            suppress_parent
            and row.get("tool", "").strip().lower() == "vicat"
            and row.get("evidence_scope", "parent_discovery").strip()
                == "parent_discovery"
        )
    ]
    return selected, region_rows, suppress_parent


def taxonomy_string(lineage: dict[str, str]) -> str:
    return ";".join(f"{PREFIXES[rank]}{lineage.get(rank) or 'unclassified'}" for rank in RANKS)


def synthetic_group_prediction(value: str) -> bool:
    return (value or "").strip().lower().startswith(("novel_", "unplaced_"))


def group_parent_context_agrees(
    group: dict[str, str], target_rank: str, strict: dict[str, str]
) -> bool:
    """Require at least one named vConTACT3 parent and agreement with strict taxonomy."""
    target_index = RANKS.index(target_rank)
    named_context: list[tuple[str, str]] = []
    for parent_rank in RANKS[1:target_index]:
        prediction = (group.get(f"{parent_rank}_prediction") or "").strip()
        if (
            prediction
            and prediction.lower() not in {"default", "singleton", "unclassified"}
            and not synthetic_group_prediction(prediction)
        ):
            named_context.append((parent_rank, prediction))
    return bool(named_context) and all(
        strict.get(rank, "") == taxon for rank, taxon in named_context
    )


def apply_vcontact3_groups(
    groups: list[dict[str, str]],
    strict: dict[str, str],
    sequence_length: int,
    minimum_length: int,
    sample_id: str,
) -> tuple[dict[str, str], list[str], list[str], list[dict[str, object]]]:
    """Add only unambiguous, ancestry-compatible project groups."""
    analysis = dict(strict)
    labels: list[str] = []
    statuses: set[str] = set()
    audit_rows: list[dict[str, object]] = []

    for rank in RANKS[1:-1]:
        candidates: list[tuple[str, dict[str, str]]] = []
        for group in groups:
            prediction = (group.get(f"{rank}_prediction") or "").strip()
            if synthetic_group_prediction(prediction):
                labels.append(f"{rank}:{prediction}")
                candidates.append((prediction, group))
        if not candidates:
            continue

        unique_predictions = sorted({prediction for prediction, _ in candidates})
        if sequence_length < minimum_length:
            reason = "below_vcontact3_minimum_length"
            statuses.add(reason)
        elif len(unique_predictions) > 1:
            reason = "ambiguous_vcontact3_groups"
            statuses.add(reason)
        elif analysis.get(rank):
            reason = "strict_rank_already_classified"
            statuses.add(reason)
        elif any(
            group_parent_context_agrees(group, rank, strict)
            for _, group in candidates
        ):
            reason = "selected_compatible_project_group"
            statuses.add(reason)
            analysis[rank] = f"viharmony_{sample_id}_{unique_predictions[0]}"
        else:
            reason = "incompatible_or_unknown_parent_context"
            statuses.add(reason)

        accepted_prediction = (
            unique_predictions[0]
            if reason == "selected_compatible_project_group"
            else ""
        )
        for prediction in unique_predictions:
            audit_rows.append({
                "tool": "vcontact3_group",
                "method_family": "gene_sharing_network",
                "rank": rank,
                "taxon": prediction,
                "tier": "auxiliary",
                "valid": reason == "selected_compatible_project_group",
                "accepted": prediction == accepted_prediction,
                "reason": reason,
            })

    if not labels:
        statuses.add("no_project_group")
    return analysis, sorted(set(labels)), sorted(statuses), audit_rows


def decide_taxonomy(
    rows: list[dict[str, str]],
    msl: set[tuple[str, ...]],
    vcontact3_minimum_length: int = 1000,
) -> tuple[dict[str, str], str, str, list[str], list[dict[str, object]]]:
    votes: list[dict[str, object]] = []
    for row in rows:
        tool = row.get("tool", "").strip().lower()
        if tool not in TAXONOMY_TOOLS or row.get("classification") != "virus":
            continue
        lineage = {rank: clean_taxon(row.get(column, ""), rank) for rank, column in RANK_COLUMNS.items()}
        lineage["domain"] = "Viruses"
        valid = lineage_is_valid(lineage, msl)
        reason = "pending" if valid else "invalid_configured_ictv_msl"
        if tool == "vcontact3":
            try:
                sequence_length = int(float(row.get("length", "0") or 0))
            except ValueError:
                sequence_length = 0
            if sequence_length < vcontact3_minimum_length:
                valid = False
                reason = "below_vcontact3_minimum_length"
        tier = vote_tier(row)
        for rank in RANKS[1:]:
            taxon = lineage.get(rank, "")
            if taxon:
                votes.append({"tool": tool, "method_family": METHOD_FAMILY[tool], "rank": rank, "taxon": taxon, "tier": tier, "valid": valid, "accepted": False, "reason": reason, "_lineage": lineage})

    accepted = {rank: "" for rank in RANKS}
    accepted["domain"] = "Viruses"
    supporting_tools: set[str] = set()
    conflict = "false"
    deepest = "domain"
    for rank in RANKS[1:]:
        earlier_ranks = RANKS[1:RANKS.index(rank)]
        rank_votes = [
            vote for vote in votes
            if vote["rank"] == rank
            and vote["valid"]
            and all(
                not accepted[ancestor]
                or vote["_lineage"].get(ancestor) == accepted[ancestor]
                for ancestor in earlier_ranks
            )
        ]
        if not rank_votes:
            break
        grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
        for vote in rank_votes:
            grouped[str(vote["taxon"])].append(vote)
        scores = {
            taxon: (
                max(TIER_SCORE[str(vote["tier"])] for vote in taxon_votes),
                len({str(vote["method_family"]) for vote in taxon_votes}),
                len({str(vote["tool"]) for vote in taxon_votes}),
            )
            for taxon, taxon_votes in grouped.items()
        }
        best_score = max(scores.values())
        winners = [taxon for taxon, score in scores.items() if score == best_score]
        if len(winners) != 1:
            conflict = "true"
            for vote in rank_votes:
                vote["reason"] = "rank_tie"
            break
        winner = winners[0]
        accepted[rank] = winner
        deepest = rank
        for vote in rank_votes:
            if vote["taxon"] == winner:
                vote["accepted"] = True
                vote["reason"] = "selected"
                supporting_tools.add(str(vote["tool"]))
            else:
                vote["reason"] = "lower_method_aware_support"

    deepest_votes = [vote for vote in votes if vote["rank"] == deepest and vote["accepted"]]
    families = {str(vote["method_family"]) for vote in deepest_votes}
    max_tier = max((TIER_SCORE[str(vote["tier"])] for vote in deepest_votes), default=0)
    if len(families) >= 2:
        confidence = "multi_method_consensus"
    elif max_tier == 3:
        confidence = "single_method_high"
    elif deepest == "domain":
        confidence = "unclassified"
    else:
        confidence = "provisional"
    return accepted, deepest, confidence, sorted(supporting_tools), votes


def original_output_id(original_id: str, record_type: str, coordinates: str) -> str:
    if record_type in {"provirus", "viral_region"}:
        return f"{original_id}|provirus_{coordinates.replace('-', '_')}"
    return original_id


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(args: argparse.Namespace) -> None:
    minimum_vcontact3_length = getattr(args, "vcontact3_min_taxonomy_length", 1000)
    if minimum_vcontact3_length < 0:
        raise ValueError("vConTACT3 minimum taxonomy length must be nonnegative")
    normalized = read_fasta(args.normalized_fasta)
    refined = read_fasta(args.refined_fasta)
    header_columns, header_rows = read_tsv(args.header_map)
    require(header_columns, {"sample_id", "sequence_id", "original_id", "original_header", "length"}, args.header_map)
    headers = {row["sequence_id"]: row for row in header_rows}
    gate_columns, gate_rows = read_tsv(args.discovery_gate)
    require(gate_columns, {"sequence_id", "advance_to_refinement", "discovery_status"}, args.discovery_gate)
    gate = {row["sequence_id"]: row for row in gate_rows}
    region_columns, region_rows = read_tsv(args.region_map)
    require(region_columns, {"sequence_id", "parent_sequence_id", "record_type", "coordinates", "original_length", "refined_length"}, args.region_map)
    regions = {row["sequence_id"]: row for row in region_rows}
    if set(refined) != set(regions):
        raise ValueError("Refined FASTA and provirus-region map contain different sequence IDs")

    final_to_parent = {sequence_id: row["parent_sequence_id"].strip() or sequence_id for sequence_id, row in regions.items()}
    known_final_ids = set(final_to_parent)
    evidence_rows: list[dict[str, str]] = []
    evidence_columns: list[str] = []
    for path in args.evidence:
        columns, rows = read_tsv(path)
        require(columns, {"sample_id", "sequence_id", "parent_sequence_id", "tool", "classification"}, path)
        evidence_columns.extend(column for column in columns if column not in evidence_columns)
        for row in rows:
            if row["sample_id"] != args.sample_id:
                raise ValueError(f"Evidence sample mismatch in {path}: {row['sample_id']}")
            row["_source_file"] = path.name
            evidence_rows.append(row)

    groups_by_sequence: dict[str, list[dict[str, str]]] = defaultdict(list)
    for path in args.vcontact3_groups:
        columns, rows = read_tsv(path)
        if "sequence_id" not in columns:
            continue
        for row in rows:
            groups_by_sequence[row["sequence_id"]].append(row)

    msl, msl_release = load_msl(args.ictv_msl)
    rows_by_final: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in evidence_rows:
        exact = row.get("sequence_id", "").strip()
        parent = row.get("parent_sequence_id", "").strip() or exact
        row_record_type = row.get("record_type", "").strip().lower()
        row_interval = parse_optional_interval(row.get("coordinates", ""))
        for final_id, final_parent in final_to_parent.items():
            # Region-specific refinement calls must never leak into a sibling
            # region from the same parent. Parent-contig discovery evidence is
            # inherited by every refined child from that parent.
            applies = evidence_applies_to_final(
                exact,
                parent,
                row_record_type,
                row_interval,
                final_id,
                final_parent,
                regions[final_id].get("record_type", "").strip(),
                parse_optional_interval(regions[final_id].get("coordinates", "")),
                len(refined[final_id]),
                known_final_ids,
            )
            if applies:
                rows_by_final[final_id].append(row)

    metadata_rows: list[dict[str, object]] = []
    review_rows: list[dict[str, object]] = []
    evidence_audit_rows: list[dict[str, object]] = []
    taxonomy_audit_rows: list[dict[str, object]] = []
    normalized_fasta_rows: list[tuple[str, str]] = []
    original_fasta_rows: list[tuple[str, str]] = []
    provisional_normalized_rows: list[tuple[str, str]] = []
    provisional_original_rows: list[tuple[str, str]] = []
    used_original_output_ids: set[str] = set()
    map_rows: list[dict[str, object]] = []

    for final_id, sequence in refined.items():
        region = regions[final_id]
        record_type = region["record_type"].strip()
        parent_id = final_to_parent[final_id]
        header = headers[parent_id]
        all_rows = rows_by_final.get(final_id, [])
        rows, region_vicat_rows, suppress_parent_vicat = select_vicat_scope(
            all_rows, final_id, parent_id
        )
        for source_row in all_rows:
            audit = dict(source_row)
            is_suppressed_parent_vicat = (
                source_row.get("tool", "").strip().lower() == "vicat"
                and source_row.get("evidence_scope", "parent_discovery").strip()
                    == "parent_discovery"
                and suppress_parent_vicat
            )
            if (
                source_row.get("tool", "").strip().lower() == "deep6"
                and source_row.get("classification", "").strip().lower() == "virus"
                and not source_row.get("evidence_strength", "").strip()
            ):
                audit["evidence_strength"] = "qualified"
                audit["strength_basis"] = "legacy_deep6_confident_prediction"
            audit["final_sequence_id"] = final_id
            audit["applied_scope"] = (
                "region" if source_row.get("sequence_id") == final_id else "parent"
            )
            audit["evidence_application"] = (
                "superseded_parent_vicat"
                if is_suppressed_parent_vicat else "used_for_final_decision"
            )
            evidence_audit_rows.append(audit)
        tools_by_class: dict[str, set[str]] = defaultdict(set)
        strong_tools: set[str] = set()
        qualified_tools: set[str] = set()
        checkv_complete_support = checkv_supports_complete_viral_contig(
            rows, record_type
        )
        for row in rows:
            tool = row.get("tool", "").strip().lower()
            classification = row.get("classification", "").strip().lower()
            effective_strength = row.get("evidence_strength", "").strip().lower()
            if tool == "deep6" and classification == "virus" and not effective_strength:
                effective_strength = "qualified"
            if tool in ORIGIN_TOOLS and effective_strength in {"qualified", "strong"}:
                tools_by_class[classification].add(tool)
                if classification == "virus":
                    if effective_strength == "strong":
                        strong_tools.add(tool)
                    elif effective_strength == "qualified":
                        qualified_tools.add(tool)

        qualified_tools.difference_update(strong_tools)
        if checkv_complete_support:
            strong_tools.add("checkv")
            qualified_tools.discard("checkv")
        if len(strong_tools) >= 2:
            viral_confidence = "high"
        elif strong_tools or len(qualified_tools) >= 2:
            viral_confidence = "supported"
        else:
            viral_confidence = "provisional"
        tesorter_summary = summarize_tesorter(
            rows,
            bool(strong_tools or qualified_tools),
            viral_confidence == "high",
        )

        strict, strict_rank, tax_confidence, tax_tools, vote_rows = decide_taxonomy(
            rows, msl, minimum_vcontact3_length
        )
        for vote in vote_rows:
            vote["sample_id"] = args.sample_id
            vote["final_sequence_id"] = final_id
            taxonomy_audit_rows.append(vote)

        analysis, group_labels, group_statuses, group_audit = apply_vcontact3_groups(
            groups_by_sequence.get(final_id, []),
            strict,
            len(sequence),
            minimum_vcontact3_length,
            args.sample_id,
        )
        for vote in group_audit:
            vote["sample_id"] = args.sample_id
            vote["final_sequence_id"] = final_id
            taxonomy_audit_rows.append(vote)

        coords = region["coordinates"].strip() or "NA"
        discovery_status = gate[parent_id]["discovery_status"].strip()
        strict_taxonomy = taxonomy_string(strict)
        analysis_taxonomy = taxonomy_string(analysis)
        original_name = header["original_id"]
        cellular = sorted(tools_by_class.get("cellular", set()))
        plasmid = sorted(tools_by_class.get("plasmid", set()))
        selected_vicat_rows = [
            row for row in rows if row.get("tool", "").strip().lower() == "vicat"
        ]
        vicat_patterns = sorted({
            row.get("origin_pattern", "").strip() for row in selected_vicat_rows
            if row.get("origin_pattern", "").strip()
        })
        vicat_scopes = sorted({
            row.get("evidence_scope", "parent_discovery").strip()
            for row in selected_vicat_rows
        })
        vicat_nonviral_classes = sorted({
            item.strip()
            for row in selected_vicat_rows
            for item in row.get("nonviral_supported_classes", "").split(",")
            if item.strip()
        })
        vicat_dominant_nonviral_classes = sorted({
            row.get("dominant_nonviral_class", "").strip()
            for row in selected_vicat_rows
            if row.get("dominant_nonviral_class", "").strip()
        })
        if "resolved_provirus_with_vicat_support" in vicat_patterns:
            vicat_provirus_status = "resolved_with_vicat_support"
        elif "localized_viral_cluster" in vicat_patterns:
            vicat_provirus_status = "potential_provirus"
        else:
            vicat_provirus_status = "none"
        viral_decision = adjudicate_viral_decision(
            discovery_status,
            record_type,
            strong_tools,
            qualified_tools,
            cellular,
            plasmid,
            tesorter_summary["sequence_interpretation"],
        )
        metadata = {
            "original_contig_name": original_name,
            "normalized_name": parent_id,
            "final_sequence_id": final_id,
            "record_type": record_type,
            "provirus_coordinates": coords,
            "original_length": region["original_length"],
            "refined_length": len(sequence),
            "discovery_status": discovery_status,
            "viral_decision": viral_decision,
            "viral_confidence": viral_confidence,
            "strong_tools": ",".join(sorted(strong_tools)) or "NA",
            "qualified_tools": ",".join(sorted(qualified_tools)) or "NA",
            "cellular_conflict_tools": ",".join(cellular) or "NA",
            "checkv_complete_viral_contig_support": str(
                checkv_complete_support
            ).lower(),
            "plasmid_conflict_tools": ",".join(plasmid) or "NA",
            "strict_taxonomy": strict_taxonomy,
            "strict_taxonomy_rank": strict_rank,
            "analysis_taxonomy": analysis_taxonomy,
            "analysis_taxonomy_rank": next((rank for rank in reversed(RANKS) if analysis.get(rank)), "domain"),
            "taxonomy_confidence": tax_confidence,
            "taxonomy_supporting_tools": ",".join(tax_tools) or "NA",
            "taxonomy_conflict": "true" if any(vote["reason"] == "rank_tie" for vote in vote_rows) else "false",
            "vcontact3_groups": ",".join(sorted(set(group_labels))) or "NA",
            "vcontact3_group_status": ",".join(group_statuses),
            "vicat_origin_pattern": ",".join(vicat_patterns) or "NA",
            "vicat_provirus_status": vicat_provirus_status,
            "vicat_evidence_scope": ",".join(vicat_scopes) or "NA",
            "vicat_viral_supported_loci": str(sum(
                int(float(row.get("viral_supported_loci", "0") or 0))
                for row in selected_vicat_rows
            )),
            "vicat_cellular_supported_loci": str(sum(
                int(float(row.get("cellular_supported_loci", "0") or 0))
                for row in selected_vicat_rows
            )),
            "vicat_nonviral_supported_classes": (
                ",".join(vicat_nonviral_classes) or "NA"
            ),
            "vicat_dominant_nonviral_class": (
                ",".join(vicat_dominant_nonviral_classes) or "NA"
            ),
            **tesorter_summary,
        }
        metadata_rows.append(metadata)
        reasons = []
        parent_vicat_viral = any(
            row.get("tool", "").strip().lower() == "vicat"
            and row.get("classification", "").strip().lower() == "virus"
            and row.get("evidence_scope", "parent_discovery").strip() == "parent_discovery"
            for row in all_rows
        )
        region_vicat_viral = any(
            row.get("classification", "").strip().lower() == "virus"
            for row in region_vicat_rows
        )
        if final_id != parent_id and parent_vicat_viral and not region_vicat_viral:
            reasons.append("provirus_boundary_vicat_discordance")
            metadata["vicat_provirus_status"] = "boundary_vicat_discordance"
        if viral_confidence == "provisional": reasons.append("single_qualified_viral_tool")
        if cellular: reasons.append("cellular_conflict")
        if plasmid: reasons.append("plasmid_conflict")
        if metadata["taxonomy_conflict"] == "true": reasons.append("taxonomy_conflict")
        if "ambiguous_vcontact3_groups" in group_statuses: reasons.append("vcontact3_group_ambiguity")
        if "incompatible_or_unknown_parent_context" in group_statuses: reasons.append("vcontact3_group_context_conflict")
        if tesorter_summary["sequence_interpretation"] == "viral_retroelement_conflict":
            reasons.append("viral_retroelement_conflict")
        elif tesorter_summary["sequence_interpretation"] == "viral_with_mobile_element_features":
            reasons.append("mobile_element_annotation")
        elif tesorter_summary["sequence_interpretation"] == "likely_retroelement":
            reasons.append("likely_retroelement")
        elif tesorter_summary["sequence_interpretation"] == "ambiguous_mobile_element":
            reasons.append("ambiguous_mobile_element")
        if strict_rank == "domain": reasons.append("taxonomy_unclassified")
        if reasons:
            review_rows.append({"final_sequence_id": final_id, "review_reasons": ",".join(reasons), **metadata})

        if viral_decision in PRIMARY_VIRAL_DECISIONS:
            normalized_fasta_rows.append((final_id, sequence))
        elif viral_decision in PROVISIONAL_VIRAL_DECISIONS:
            provisional_normalized_rows.append((final_id, sequence))
        original_output = original_output_id(original_name, record_type, coords)
        if original_output in used_original_output_ids:
            original_output = f"{original_output}|visum_{parent_id}"
        used_original_output_ids.add(original_output)
        if viral_decision in PRIMARY_VIRAL_DECISIONS:
            original_fasta_rows.append((original_output, sequence))
        elif viral_decision in PROVISIONAL_VIRAL_DECISIONS:
            provisional_original_rows.append((original_output, sequence))
        map_rows.append({
            "original_contig_name": original_name, "original_header": header["original_header"],
            "normalized_parent_id": parent_id, "final_sequence_id": final_id,
            "final_original_id": original_output, "record_type": record_type,
            "provirus_coordinates": coords,
        })

    prefix = args.output_prefix
    outputs = {
        "normalized_fasta": Path(f"{prefix}.final.normalized.fasta"),
        "original_fasta": Path(f"{prefix}.final.original_ids.fasta"),
        "provisional_normalized_fasta": Path(
            f"{prefix}.provisional.normalized.fasta"
        ),
        "provisional_original_fasta": Path(
            f"{prefix}.provisional.original_ids.fasta"
        ),
        "provisional_metadata": Path(f"{prefix}.provisional_metadata.tsv"),
        "metadata": Path(f"{prefix}.final_metadata.tsv"),
        "review": Path(f"{prefix}.review_queue.tsv"),
        "disposition": Path(f"{prefix}.sequence_disposition.tsv"),
        "map": Path(f"{prefix}.sequence_map.tsv"),
        "database_fasta": Path(f"{prefix}.database_candidates.fasta"),
        "database_tsv": Path(f"{prefix}.database_candidates.tsv"),
        "all_candidates_fasta": Path(f"{prefix}.all_candidates.fasta"),
        "review_candidates_fasta": Path(f"{prefix}.review_candidates.fasta"),
    }
    write_fasta(outputs["normalized_fasta"], normalized_fasta_rows)
    write_fasta(outputs["original_fasta"], original_fasta_rows)
    write_fasta(outputs["provisional_normalized_fasta"], provisional_normalized_rows)
    write_fasta(outputs["provisional_original_fasta"], provisional_original_rows)
    write_tsv(outputs["metadata"], METADATA_COLUMNS, metadata_rows)
    write_tsv(outputs["review"], ["review_reasons", *METADATA_COLUMNS], review_rows)
    write_tsv(outputs["map"], ["original_contig_name", "original_header", "normalized_parent_id", "final_sequence_id", "final_original_id", "record_type", "provirus_coordinates"], map_rows)
    write_fasta(outputs["database_fasta"], normalized_fasta_rows)
    primary_ids = {identifier for identifier, _ in normalized_fasta_rows}
    primary_metadata = [
        row for row in metadata_rows if row["final_sequence_id"] in primary_ids
    ]
    provisional_ids = {
        identifier for identifier, _ in provisional_normalized_rows
    }
    provisional_metadata = [
        row for row in metadata_rows if row["final_sequence_id"] in provisional_ids
    ]
    write_tsv(outputs["provisional_metadata"], METADATA_COLUMNS, provisional_metadata)
    write_tsv(outputs["database_tsv"], METADATA_COLUMNS, primary_metadata)
    write_fasta(outputs["all_candidates_fasta"], list(refined.items()))
    review_ids = {str(row["final_sequence_id"]) for row in review_rows}
    write_fasta(
        outputs["review_candidates_fasta"],
        [
            (identifier, sequence)
            for identifier, sequence in refined.items()
            if identifier in review_ids
        ],
    )

    final_by_parent: dict[str, list[str]] = defaultdict(list)
    provisional_by_parent: dict[str, list[str]] = defaultdict(list)
    decision_by_final = {
        str(row["final_sequence_id"]): str(row["viral_decision"])
        for row in metadata_rows
    }
    for row in map_rows:
        final_id = str(row["final_sequence_id"])
        if decision_by_final.get(final_id) in PRIMARY_VIRAL_DECISIONS:
            final_by_parent[str(row["normalized_parent_id"])].append(final_id)
        elif decision_by_final.get(final_id) in PROVISIONAL_VIRAL_DECISIONS:
            provisional_by_parent[str(row["normalized_parent_id"])].append(final_id)
    disposition_rows = []
    for sequence_id in normalized:
        gate_row = gate[sequence_id]
        final_ids = final_by_parent.get(sequence_id, [])
        provisional_ids_for_parent = provisional_by_parent.get(sequence_id, [])
        if final_ids:
            disposition = "retained_as_provirus" if any("|viral_region_" in item for item in final_ids) else "retained_contig"
        elif provisional_ids_for_parent:
            disposition = (
                "provisional_provirus_after_harmony"
                if any("|viral_region_" in item for item in provisional_ids_for_parent)
                else "provisional_viral_after_harmony"
            )
        elif gate_row["advance_to_refinement"].strip().lower() == "true":
            candidate_decisions = [
                decision_by_final.get(str(row["final_sequence_id"]), "")
                for row in map_rows
                if str(row["normalized_parent_id"]) == sequence_id
            ]
            if "likely_nonviral" in candidate_decisions:
                disposition = "likely_nonviral_after_harmony"
            elif "likely_retroelement" in candidate_decisions:
                disposition = "likely_retroelement_after_harmony"
            else:
                disposition = "ambiguous_review_after_harmony"
        elif gate_row["discovery_status"] == "likely_nonviral":
            disposition = "discovery_noncandidate_conflicting_origin"
        else:
            disposition = "discovery_noncandidate_insufficient_support"
        disposition_rows.append({"original_contig_name": headers[sequence_id]["original_id"], "normalized_name": sequence_id, "disposition": disposition, "final_sequence_ids": ",".join(final_ids) or "NA", "provisional_sequence_ids": ",".join(provisional_ids_for_parent) or "NA"})
    write_tsv(outputs["disposition"], ["original_contig_name", "normalized_name", "disposition", "final_sequence_ids", "provisional_sequence_ids"], disposition_rows)

    if args.audit_mode == "full":
        evidence_path = Path(f"{prefix}.evidence_audit.tsv.gz")
        audit_columns = [
            "final_sequence_id", "applied_scope", "evidence_application",
            "_source_file", *evidence_columns,
        ]
        with gzip.open(evidence_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=audit_columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader(); writer.writerows(evidence_audit_rows)
        taxonomy_path = Path(f"{prefix}.taxonomy_vote_audit.tsv.gz")
        tax_columns = ["sample_id", "final_sequence_id", "tool", "method_family", "rank", "taxon", "tier", "valid", "accepted", "reason"]
        with gzip.open(taxonomy_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=tax_columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader(); writer.writerows(taxonomy_audit_rows)
        outputs["evidence_audit"] = evidence_path
        outputs["taxonomy_audit"] = taxonomy_path

    manifest_path = Path(f"{prefix}.harmonizer_manifest.json")
    manifest = {
        "schema_version": "viharmony-0.6",
        "sample_id": args.sample_id,
        "input_type": args.input_type,
        "ictv_msl": args.ictv_msl.name,
        "ictv_msl_release": msl_release,
        "audit_mode": args.audit_mode,
        "tools_present": sorted({row.get("tool", "") for row in evidence_rows if row.get("tool")}),
        "policy": {
            "viral_confidence": "high=2+ strong; supported=1 strong or 2+ qualified; provisional=1 qualified",
            "checkv_complete_contig": "complete/high-quality, provirus=No, viral genes present, and zero host genes is strong intact-contig viral support",
            "taxonomy": "configured-ICTV-MSL-validated method-aware rank voting",
            "vcontact3_min_taxonomy_length": minimum_vcontact3_length,
            "vcontact3_project_groups": "ancestry-compatible, unambiguous groups only",
            "vicat_scope": "refined-region evidence supersedes parent-discovery evidence",
            "provirus_adjudication": "refinement selects coordinates; origin and mobile-element conflicts still control final retention",
            "tesorter_retrovirus": "explicit viral-like labels remain compatible; high-confidence multi-tool viral consensus overrides generic mobile-element conflict while retaining annotation",
            "final_output": "primary requires strong or multi-tool qualified viral support; single-tool qualified calls are emitted separately as provisional",
        },
        "outputs": {},
    }
    for label, path in outputs.items():
        manifest["outputs"][label] = {"file": path.name, "sha256": sha256(path)}
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"viHARMONY sample={args.sample_id} candidates={len(metadata_rows)} "
        f"retained={len(normalized_fasta_rows)} review={len(review_rows)}"
    )


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, ValueError, csv.Error) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
