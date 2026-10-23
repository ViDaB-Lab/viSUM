#!/usr/bin/env python3

"""Combine viSUM evidence into final sequences, decisions, and taxonomy calls."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import multiprocessing
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager, ExitStack
from dataclasses import dataclass, field
from functools import lru_cache
from heapq import merge
from itertools import chain, islice
from pathlib import Path
from typing import Iterable, TextIO

from evidence_schema import is_review_only_evidence
from evidence_store import EvidenceIndex, EvidenceTable
from ictv_reconciliation import Reconciler
import te_routing
import reporting_bundle


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
GENERIC_ORIGIN_CLASSIFIERS = {"deep6", "deepmicroclass2"}
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
    "parent_cellular_context_tools", "parent_plasmid_context_tools",
    "checkv_complete_viral_contig_support", "generic_cellular_conflict_overridden",
    "plasmid_conflict_tools", "strict_taxonomy", "strict_taxonomy_rank",
    "strict_taxonomy_confidence", "strict_taxonomy_supporting_tools",
    "strict_taxonomy_supporting_method_families", "strict_taxonomy_stop_reason",
    "strict_taxonomy_conflict_rank", "strict_taxonomy_conflicting_taxa",
    "exploratory_taxonomy", "exploratory_taxonomy_rank",
    "exploratory_taxonomy_confidence", "exploratory_taxonomy_supporting_tools",
    "exploratory_taxonomy_supporting_method_families",
    "exploratory_taxonomy_status", "strict_species", "exploratory_species",
    "species_assignment_status", "species_assignment_basis",
    "species_conflicting_taxa", "taxonomy_conflict", "vcontact3_groups",
    "vcontact3_group_status", "sequence_interpretation", "tesorter_status",
    "tesorter_evidence_strength", "tesorter_categories", "tesorter_orders",
    "tesorter_superfamilies", "tesorter_assignment_methods",
    "vicat_origin_pattern", "vicat_provirus_status", "vicat_evidence_scope",
    "vicat_viral_supported_loci", "vicat_cellular_supported_loci",
    "vicat_nonviral_supported_classes", "vicat_dominant_nonviral_class",
    "provirus_boundary_source", "provirus_boundary_status",
    "provirus_boundary_supporting_tools", "provirus_retention_basis",
    "viral_entity_interpretation",
    "rna_checkv_conflict_exception", "rna_checkv_exception_supporting_loci",
]

METADATA_COLUMNS += te_routing.METADATA_COLUMNS

STRENGTH_ORDER = {"": 0, "weak": 1, "qualified": 2, "strong": 3}

PRIMARY_VIRAL_DECISIONS = {"retained_viral", "retained_provirus"}
PROVISIONAL_VIRAL_DECISIONS = {"provisional_viral", "provisional_provirus"}

# Preserve the complete pre-cap consensus, not only individual tool predictions.
REPORTING_TAXONOMY_COLUMNS = [
    name for name in METADATA_COLUMNS
    if name.startswith(("strict_taxonomy", "exploratory_taxonomy", "species_"))
    or name in {"strict_species", "exploratory_species", "taxonomy_conflict"}
]
REPORTING_AUDIT_COLUMNS = [
    "sample_id", "final_sequence_id", "harmonizer_max_rank", "reporting_status",
    "reported_strict_taxonomy", "reported_strict_rank",
    "reported_exploratory_taxonomy", "reported_exploratory_rank",
    *["uncapped_" + name for name in REPORTING_TAXONOMY_COLUMNS],
]


@dataclass
class MslCatalog:
    """Indexed ICTV lineages used to validate and resolve partial predictions."""

    lineages: tuple[tuple[str, ...], ...]
    by_taxon: dict[tuple[str, str], frozenset[int]]

    def matching(self, lineage: dict[str, str]) -> frozenset[int]:
        named = [
            (rank, lineage.get(rank, ""))
            for rank in RANKS[1:]
            if lineage.get(rank)
        ]
        if not named:
            return frozenset()
        matches: set[int] | None = None
        for rank, taxon in named:
            indexed = self.by_taxon.get((rank, taxon), frozenset())
            matches = set(indexed) if matches is None else matches.intersection(indexed)
            if not matches:
                return frozenset()
        return frozenset(matches or set())


@dataclass
class TaxonomyDecision:
    strict: dict[str, str]
    exploratory: dict[str, str]
    strict_rank: str
    exploratory_rank: str
    strict_confidence: str
    exploratory_confidence: str
    strict_tools: list[str]
    exploratory_tools: list[str]
    strict_families: list[str]
    exploratory_families: list[str]
    strict_stop_reason: str
    conflict_rank: str
    conflicting_taxa: list[str]
    exploratory_status: str
    species_status: str
    species_basis: str
    species_conflicting_taxa: list[str]
    votes: list[dict[str, object]]


@dataclass(frozen=True)
class CandidateContext:
    sample_id: str
    input_type: str
    minimum_vcontact3_length: int
    pair_floor: bool
    has_locus_evidence: bool
    disable_rna_checkv_exception: bool
    full_audit: bool
    reconciliation_bundle: str | None = None
    reconciliation_sha256: str | None = None
    max_rank: str = "family"


@dataclass(frozen=True)
class CandidateTask:
    final_id: str
    length: int
    region: dict[str, str]
    parent_id: str
    header: dict[str, str]
    discovery_status: str
    evidence: list[dict[str, str]]
    groups: list[dict[str, str]]
    loci: list[dict[str, str]] | tuple


@dataclass
class CandidateResult:
    metadata: dict[str, object]
    review: dict[str, object] | None
    evidence_audit: list[dict[str, object]]
    taxonomy_audit: list[dict[str, object]]
    pair_floor_status: str
    rna_checkv_exception: bool
    reconciliation_audit: list[dict[str, object]] = field(default_factory=list)
    reporting_audit: dict[str, object] = field(default_factory=dict)
    routing_audit: dict[str, object] = field(default_factory=dict)


CANDIDATE_BATCH_SIZE = 16
_WORKER_CONTEXT: CandidateContext | None = None
_WORKER_MSL: MslCatalog | None = None


def _initialize_candidate_worker(context: CandidateContext, msl: MslCatalog) -> None:
    global _WORKER_CONTEXT, _WORKER_MSL
    # Only the small configuration/catalog are initialized in each worker.
    # Entire FASTAs, evidence tables, and output collections stay in the parent.
    _WORKER_CONTEXT, _WORKER_MSL = context, msl


def _evaluate_candidate_batch(batch: list[CandidateTask]) -> list[CandidateResult]:
    if _WORKER_CONTEXT is None or _WORKER_MSL is None:
        raise RuntimeError("viHARMONY candidate worker was not initialized")
    return [evaluate_candidate(task, _WORKER_CONTEXT, _WORKER_MSL) for task in batch]


def _ordered_candidate_batches(executor, tasks, batch_size, pending_limit):
    """At most pending_limit batches exist in flight, including finished ones."""
    tasks = iter(tasks)
    pending = deque()
    try:
        for _ in range(pending_limit):
            batch = list(islice(tasks, batch_size))
            if not batch:
                break
            pending.append(executor.submit(_evaluate_candidate_batch, batch))
        while pending:
            # Consume in submission order, not completion order. This also
            # bounds buffering when an early candidate takes longer to evaluate.
            yield from pending.popleft().result()
            batch = list(islice(tasks, batch_size))
            if batch:
                pending.append(executor.submit(_evaluate_candidate_batch, batch))
    finally:
        for future in pending:
            future.cancel()


@contextmanager
def candidate_results(tasks, context, msl, workers, batch_size=CANDIDATE_BATCH_SIZE):
    if workers == 1:
        yield (evaluate_candidate(task, context, msl) for task in tasks)
        return
    # Explicit spawn avoids inheriting the parent dataset and behaves the same
    # on Linux and Windows. Never eagerly enqueue a million-candidate iterable.
    executor = ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
        initializer=_initialize_candidate_worker, initargs=(context, msl),
    )
    results = _ordered_candidate_batches(executor, tasks, batch_size, 2 * workers)
    try:
        yield results
    finally:
        results.close()
        executor.shutdown(wait=True, cancel_futures=True)


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
    if row_record_type == "short_hallmark_region" and row_interval is None:
        # An unlocalized hallmark cannot corroborate arbitrary refined children.
        return final_record_type == "input_contig" and parent_id == final_parent
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


def iter_evidence_by_final(evidence_rows, final_to_parent, regions, refined):
    """Yield matches in original row/final order, without retaining row copies."""
    finals_by_parent = defaultdict(list)
    final_order = {}
    final_properties = {}
    for position, (final_id, parent_id) in enumerate(final_to_parent.items()):
        finals_by_parent[parent_id].append(final_id)
        final_order[final_id] = position
        final_properties[final_id] = (
            regions[final_id].get("record_type", "").strip(),
            parse_optional_interval(regions[final_id].get("coordinates", "")),
            len(refined[final_id]),
        )
    known_final_ids = set(final_to_parent)
    for row in evidence_rows:
        exact = row.get("sequence_id", "").strip()
        parent = row.get("parent_sequence_id", "").strip() or exact
        row_record_type = row.get("record_type", "").strip().lower()
        row_interval = parse_optional_interval(row.get("coordinates", ""))
        if exact in known_final_ids:
            # Exact-final matches take precedence, even over a different parent.
            candidates = (exact,)
        elif parent == exact:
            candidates = finals_by_parent.get(parent, ())
        else:
            # Both parent IDs can match the existing generic-evidence rule.
            # Distinct parent buckets are disjoint; merge in original final order.
            candidates = merge(
                finals_by_parent.get(parent, ()),
                finals_by_parent.get(exact, ()),
                key=final_order.__getitem__,
            )
        for final_id in candidates:
            record_type, interval, length = final_properties[final_id]
            if evidence_applies_to_final(
                exact, parent, row_record_type, row_interval,
                final_id, final_to_parent[final_id], record_type, interval,
                length, known_final_ids,
            ):
                yield final_id, row


def index_evidence_by_final(evidence_rows, final_to_parent, regions, refined):
    """Dictionary index retained for callers needing original row identities."""
    rows_by_final = defaultdict(list)
    for final_id, row in iter_evidence_by_final(evidence_rows, final_to_parent, regions, refined):
        rows_by_final[final_id].append(row)
    return rows_by_final


def rna_floor_qualified_tools(rows, qualified_tools):
    """Hallmark-only VS2 support must not exempt a Deep6/viCAT pair."""
    viral_vs2 = [row for row in rows if row.get('tool', '').lower() == 'virsorter2'
                 and row.get('classification', '').lower() == 'virus']
    if viral_vs2 and all(row.get('strength_basis') == 'virsorter2_hallmark_only'
                         for row in viral_vs2):
        return qualified_tools - {'virsorter2'}
    return set(qualified_tools)


def apply_rna_pair_homology_floor(decision, enabled, input_type, strong_tools,
                                 qualified_tools, loci, parent_id, coordinates=""):
    """Restrict only primary RNA Deep6/viCAT-only calls; preserve other decisions."""
    if not enabled or input_type != "rna" or strong_tools or qualified_tools != {"deep6", "vicat"} or decision not in PRIMARY_VIRAL_DECISIONS:
        return decision, "not_applicable"
    interval = tuple(map(int, coordinates.split("-"))) if coordinates else None
    for row in loci:
        if row["sequence_id"] != parent_id or row["locus_classification"] != "viral_supported":
            continue
        if interval:
            start, stop = map(int, row["coordinates"].split("-"))
            if not (interval[0] <= start <= stop <= interval[1]):
                continue
        if (float(row["protein_length"]) >= 100 and float(row["best_bitscore"]) >= 100
                and float(row["best_query_coverage"]) >= 50
                and float(row["best_subject_coverage"]) >= 50):
            return decision, "passed"
    return ("provisional_provirus" if decision == "retained_provirus" else "provisional_viral"), "demoted"


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

    generic_cellular_override = (
        len(strong_tools) >= 2
        and bool(cellular_tools)
        and set(cellular_tools).issubset(GENERIC_ORIGIN_CLASSIFIERS)
        and not plasmid_tools
    )
    conflicts = bool(cellular_tools or plasmid_tools)
    if conflicts:
        if generic_cellular_override:
            return "retained_provirus" if refined_region else "retained_viral"
        if (
            not refined_region
            and plasmid_tools
            and not cellular_tools
            and len(strong_tools) >= 2
        ):
            # Strong, independent viral support can coexist with an explicit
            # plasmid call for phage-plasmids and related hybrid replicons.
            return "retained_viral"
        if discovery_status == "ambiguous" and not strong_tools and len(qualified_tools) <= 1:
            return "likely_nonviral"
        return "ambiguous_review"
    if strong_tools or len(qualified_tools) >= 2:
        return "retained_provirus" if refined_region else "retained_viral"
    if qualified_tools:
        return "provisional_provirus" if refined_region else "provisional_viral"
    return "ambiguous_review"


def rna_checkv_exception_loci(input_type, decision, record_type, cellular, plasmid,
                            supported_tools, interpretation, loci, parent_id, length):
    """Balanced RNA-only exception; keep the original conflict in output metadata."""
    if (input_type != "rna" or decision != "ambiguous_review" or record_type != "input_contig"
        or set(cellular) != {"checkv"} or plasmid
        or not {"deep6", "virbot", "vicat"}.issubset(supported_tools)
        or interpretation in {"likely_retroelement", "viral_retroelement_conflict", "ambiguous_mobile_element"}):
        return []
    selected = [r for r in loci if r["sequence_id"] == parent_id]
    if any(r["locus_classification"] in {"cellular_supported", "ambiguous"} for r in selected):
        return []
    matches = []
    for row in selected:
        if row["locus_classification"] != "viral_supported":
            continue
        def number(key):
            value = row.get(key, "")
            return float(value) if value not in {"", "NA", None} else 0.0
        start, end = map(int, row["coordinates"].split("-"))
        if not 1 <= start <= end <= length:
            raise ValueError("RNA CheckV exception locus exceeds parent")
        if (number("protein_length") >= 100 and number("best_bitscore") >= 100
            and number("best_identity") >= 50
            and min(number("best_query_coverage"), number("best_subject_coverage")) >= 80
            and (end - start + 1) / length >= 0.70):
            matches.append(row.get("locus_id", row["coordinates"]))
    return sorted(set(matches))


def adjudicate_disputed_ct3_boundary(
    region: dict[str, str],
    rows: list[dict[str, str]],
    vicat_rows: list[dict[str, str]],
) -> tuple[bool, str]:
    """Decide whether a conflicting CT3/viCAT region remains primary.

    CT3-only boundaries remain governed by the ordinary evidence policy. This
    safeguard applies only when refinement explicitly records disagreement
    between a CT3-selected boundary and viCAT's projected region.
    """
    record_type = region.get("record_type", "").strip().lower()
    boundary_source = region.get("boundary_source", "").strip().lower()
    boundary_status = region.get("boundary_status", "").strip().lower()
    supporting_tools = {
        tool.strip().lower()
        for tool in region.get("supporting_boundary_tools", "").split(",")
        if tool.strip()
    }
    disputed_ct3_vicat = (
        record_type in {"provirus", "viral_region"}
        and boundary_source == "cenotetaker3"
        and boundary_status == "selected_boundary_conflict"
        and "vicat" in supporting_tools
    )
    if not disputed_ct3_vicat:
        return True, "not_disputed_ct3_vicat_boundary"

    ct3_strong = any(
        row.get("tool", "").strip().lower() == "cenotetaker3"
        and row.get("classification", "").strip().lower() == "virus"
        and row.get("evidence_strength", "").strip().lower() == "strong"
        and row.get("record_type", "").strip().lower()
            in {"provirus", "viral_region"}
        for row in rows
    )
    if ct3_strong:
        return True, "ct3_strong_multi_hallmark_rescue"

    viral_loci = sum(
        int(float(row.get("viral_supported_loci", "0") or 0))
        for row in vicat_rows
    )
    cellular_loci = sum(
        int(float(row.get("cellular_supported_loci", "0") or 0))
        for row in vicat_rows
    )
    informative_loci = viral_loci + cellular_loci
    viral_fraction = viral_loci / informative_loci if informative_loci else 0.0
    if viral_loci >= 2 and viral_fraction >= 0.60:
        return True, "vicat_regional_support_2_loci_fraction_ge_0.60"
    return False, "provisional_disputed_ct3_vicat_boundary"


def viral_entity_interpretation(
    viral_decision: str,
    record_type: str,
    parent_plasmid_tools: list[str],
    exact_plasmid_tools: list[str],
    has_vicat_plasmid_loci: bool,
    strong_tools: set[str],
) -> str:
    """Describe the retained entity without conflating context and origin."""
    refined_region = record_type in {"provirus", "viral_region"}
    primary = viral_decision in PRIMARY_VIRAL_DECISIONS
    provisional = viral_decision in PROVISIONAL_VIRAL_DECISIONS
    if refined_region and (parent_plasmid_tools or exact_plasmid_tools):
        if primary:
            return "plasmid_associated_provirus"
        if provisional:
            return "provisional_plasmid_associated_provirus"
        return "nonprimary_plasmid_associated_region"
    if (
        not refined_region
        and primary
        and len(strong_tools) >= 2
        and (exact_plasmid_tools or has_vicat_plasmid_loci)
    ):
        return "virus_plasmid_hybrid_candidate"
    if refined_region:
        if primary:
            return "provirus"
        if provisional:
            return "provisional_provirus"
    elif primary:
        return "viral_contig"
    elif provisional:
        return "provisional_viral_contig"
    return "nonprimary_or_review"


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
    parser.add_argument("--ictv-reconciliation", type=Path, default=argparse.SUPPRESS,
                        help="Pinned offline reconciliation bundle; taxonomy fields only (omitting preserves legacy behavior)")
    parser.add_argument("--evidence", nargs="*", default=[], type=Path)
    parser.add_argument("--vcontact3-groups", nargs="*", default=[], type=Path)
    parser.add_argument("--vcontact3-min-taxonomy-length", type=int, default=1000)
    parser.add_argument("--audit-mode", choices=("none", "compact", "full"), default="compact")
    parser.add_argument("--max-rank", choices=RANKS, default="family",
                        help="Maximum reported taxonomy rank (default: family); inference, retention and uncapped audits are unchanged")
    parser.add_argument("--workers", type=int, default=1,
                        help="Bounded candidate-evaluation processes (default: 1); decisions and output order are unchanged")
    parser.add_argument("--rna-pair-homology-floor", action="store_true",
                        help="Experimental RNA-only Deep6/viCAT primary-retention requirement")
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--disable-rna-checkv-exception", action="store_true",
                        help="Disable the audited RNA-only CheckV conflict exception for controlled comparisons")
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


def make_msl_catalog(lineages: Iterable[tuple[str, ...]]) -> MslCatalog:
    ordered = tuple(sorted(set(lineages)))
    by_taxon: dict[tuple[str, str], set[int]] = defaultdict(set)
    for index, lineage in enumerate(ordered):
        for rank, taxon in zip(RANKS[1:], lineage):
            if taxon:
                by_taxon[(rank, taxon)].add(index)
    return MslCatalog(
        ordered,
        {key: frozenset(indices) for key, indices in by_taxon.items()},
    )


def load_msl(path: Path) -> tuple[MslCatalog, str]:
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
    return make_msl_catalog(lineages), release


def as_msl_catalog(msl: MslCatalog | set[tuple[str, ...]]) -> MslCatalog:
    return msl if isinstance(msl, MslCatalog) else make_msl_catalog(msl)


def lineage_is_valid(
    lineage: dict[str, str], msl: MslCatalog | set[tuple[str, ...]]
) -> bool:
    return bool(as_msl_catalog(msl).matching(lineage))


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


def partition_refined_region_context(
    rows: list[dict[str, str]],
    final_id: str,
    record_type: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Keep parent nonviral calls as context after a partial-region refinement.

    A selected partial provirus/viral region no longer contains the parent flanks.
    Consequently, a cellular or plasmid call made only on the parent contig must
    remain auditable but cannot veto the extracted region. Exact region evidence
    (including projected viCAT evidence) remains active in adjudication.
    """
    if record_type.strip().lower() not in {"provirus", "viral_region"}:
        return rows, []
    decision_rows: list[dict[str, str]] = []
    context_rows: list[dict[str, str]] = []
    for row in rows:
        classification = row.get("classification", "").strip().lower()
        exact_region_evidence = row.get("sequence_id", "").strip() == final_id
        if classification in {"cellular", "plasmid"} and not exact_region_evidence:
            context_rows.append(row)
        else:
            decision_rows.append(row)
    return decision_rows, context_rows


def taxonomy_string(lineage: dict[str, str]) -> str:
    return ";".join(f"{PREFIXES[rank]}{lineage.get(rank) or 'unclassified'}" for rank in RANKS)


def apply_taxonomy_reporting_cap(
    metadata: dict[str, object], taxonomy: TaxonomyDecision, max_rank: str,
    sample_id: str,
) -> dict[str, object]:
    """Mask only final reporting fields, after every biological/review decision.

    Votes and the TaxonomyDecision are never mutated. An always-written audit
    preserves the uncapped consensus even with audit-mode=none or compact.
    """
    if max_rank not in RANKS:
        raise ValueError("viHARMONY maximum reporting rank must be one of: " + ", ".join(RANKS))
    audit = {"sample_id": sample_id, "final_sequence_id": metadata["final_sequence_id"],
             "harmonizer_max_rank": max_rank,
             **{"uncapped_" + name: metadata[name] for name in REPORTING_TAXONOMY_COLUMNS}}
    cutoff = RANKS.index(max_rank)
    changed = False
    for mode, lineage in (("strict", taxonomy.strict), ("exploratory", taxonomy.exploratory)):
        reported = {rank: value if RANKS.index(rank) <= cutoff else ""
                    for rank, value in lineage.items()}
        if reported == lineage:
            continue
        changed = True
        rank = next((r for r in reversed(RANKS) if reported.get(r)), "domain")
        metadata[mode + "_taxonomy"] = taxonomy_string(reported)
        metadata[mode + "_taxonomy_rank"] = rank
        # Same explicit rank/taxon support used by the inference engine. Do not
        # relabel species-level confidence as family-level confidence. When a
        # parent rank is filled only by the MSL path, say so explicitly.
        supporting = [vote for vote in taxonomy.votes
                      if vote["valid"] and vote["rank"] == rank
                      and vote["taxon"] == reported.get(rank)]
        tools = sorted({str(vote["tool"]) for vote in supporting})
        families = sorted({str(vote["method_family"]) for vote in supporting})
        tier = max((TIER_SCORE[str(vote["tier"])] for vote in supporting), default=0)
        metadata[mode + "_taxonomy_supporting_tools"] = ",".join(tools) or "NA"
        metadata[mode + "_taxonomy_supporting_method_families"] = ",".join(families) or "NA"
        if rank == "domain":
            confidence = "unclassified" if mode == "strict" else "ambiguous"
        elif not supporting:
            confidence = "lineage_inferred"
        elif mode == "strict":
            confidence = ("multi_method_consensus" if len(families) >= 2 else
                          "single_method_high" if tier == TIER_SCORE["high"] else "moderate_supported")
        else:
            confidence = {3: "high", 2: "moderate", 1: "provisional", 0: "ambiguous"}[tier]
        metadata[mode + "_taxonomy_confidence"] = confidence
        if mode == "strict":
            metadata["strict_taxonomy_stop_reason"] = "reporting_rank_cap"
    if cutoff < RANKS.index("species"):
        metadata["strict_species"] = "NA"
        metadata["exploratory_species"] = "NA"
        metadata["species_conflicting_taxa"] = "NA"
        if (taxonomy.strict.get("species") or taxonomy.exploratory.get("species")
                or taxonomy.species_conflicting_taxa):
            metadata["species_assignment_status"] = "withheld_by_reporting_cap"
            metadata["species_assignment_basis"] = "maximum_reported_rank=" + max_rank
    if taxonomy.conflict_rank in RANKS and RANKS.index(taxonomy.conflict_rank) > cutoff:
        # Preserve the conflict flag/rank and review reason, but put deeper taxon
        # names in the audit rather than suggesting reportable assignments.
        metadata["strict_taxonomy_conflicting_taxa"] = "NA"
    if changed and metadata["strict_taxonomy"] == metadata["exploratory_taxonomy"]:
        metadata["exploratory_taxonomy_status"] = "matches_strict"
    audit.update({"reporting_status": "capped" if changed else "unchanged",
                  "reported_strict_taxonomy": metadata["strict_taxonomy"],
                  "reported_strict_rank": metadata["strict_taxonomy_rank"],
                  "reported_exploratory_taxonomy": metadata["exploratory_taxonomy"],
                  "reported_exploratory_rank": metadata["exploratory_taxonomy_rank"]})
    return audit


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
    """Report project groups without inserting them into formal ICTV taxonomy."""
    exploratory = dict(strict)
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
        elif exploratory.get(rank):
            reason = "strict_rank_already_classified"
            statuses.add(reason)
        elif any(
            group_parent_context_agrees(group, rank, strict)
            for _, group in candidates
        ):
            reason = "reported_compatible_project_group"
            statuses.add(reason)
        else:
            reason = "incompatible_or_unknown_parent_context"
            statuses.add(reason)

        accepted_prediction = (
            unique_predictions[0]
            if reason == "reported_compatible_project_group"
            else ""
        )
        for prediction in unique_predictions:
            audit_rows.append({
                "tool": "vcontact3_group",
                "method_family": "gene_sharing_network",
                "rank": rank,
                "taxon": prediction,
                "tier": "auxiliary",
                "valid": reason == "reported_compatible_project_group",
                "accepted": prediction == accepted_prediction,
                "reason": reason,
            })

    if not labels:
        statuses.add("no_project_group")
    return exploratory, sorted(set(labels)), sorted(statuses), audit_rows


def decide_taxonomy(
    rows: list[dict[str, str]],
    msl: MslCatalog | set[tuple[str, ...]],
    vcontact3_minimum_length: int = 1000,
) -> TaxonomyDecision:
    """Resolve taxonomy globally against complete ICTV MSL paths.

    Strict species calls require either two independent method families or one
    high-tier call without any named alternative. Single moderate/provisional
    calls remain available only in the exploratory taxonomy.
    """
    catalog = as_msl_catalog(msl)
    votes: list[dict[str, object]] = []
    claims: list[dict[str, object]] = []
    for row in rows:
        tool = row.get("tool", "").strip().lower()
        if tool not in TAXONOMY_TOOLS or row.get("classification") != "virus":
            continue
        lineage = {rank: clean_taxon(row.get(column, ""), rank) for rank, column in RANK_COLUMNS.items()}
        lineage["domain"] = "Viruses"
        matching_paths = catalog.matching(lineage)
        valid = bool(matching_paths)
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
        claims.append({
            "tool": tool,
            "method_family": METHOD_FAMILY[tool],
            "tier": tier,
            "lineage": lineage,
            "matching_paths": matching_paths,
            "valid": valid,
        })
        for rank in RANKS[1:]:
            taxon = lineage.get(rank, "")
            if taxon:
                votes.append({
                    "tool": tool, "method_family": METHOD_FAMILY[tool],
                    "rank": rank, "taxon": taxon, "tier": tier,
                    "valid": valid, "accepted": False, "reason": reason,
                })

    valid_claims = [claim for claim in claims if claim["valid"]]
    candidate_indices = sorted({
        index for claim in valid_claims for index in claim["matching_paths"]
    })
    blank = {rank: "" for rank in RANKS}
    blank["domain"] = "Viruses"
    if not candidate_indices:
        return TaxonomyDecision(
            blank, dict(blank), "domain", "domain", "unclassified",
            "ambiguous", [], [], [], [], "no_valid_taxonomy_evidence",
            "NA", [], "unclassified", "no_taxonomy_evidence",
            "no_valid_taxonomy_evidence", [], votes,
        )

    def deepest_claim_rank(claim: dict[str, object]) -> int:
        lineage = claim["lineage"]
        return max(
            (RANKS.index(rank) for rank in RANKS[1:] if lineage.get(rank)),
            default=0,
        )

    def path_score(index: int) -> tuple[int, int, int, int]:
        supporters = [
            claim for claim in valid_claims if index in claim["matching_paths"]
        ]
        return (
            max(TIER_SCORE[str(claim["tier"])] for claim in supporters),
            len({str(claim["method_family"]) for claim in supporters}),
            sum(
                TIER_SCORE[str(claim["tier"])] * deepest_claim_rank(claim)
                for claim in supporters
            ),
            len({str(claim["tool"]) for claim in supporters}),
        )

    scores = {index: path_score(index) for index in candidate_indices}
    best_score = max(scores.values())
    best_indices = [index for index, score in scores.items() if score == best_score]
    best_supporters = [
        claim for claim in valid_claims
        if any(index in claim["matching_paths"] for index in best_indices)
    ]
    deepest_supported_index = max(
        (deepest_claim_rank(claim) for claim in best_supporters), default=0
    )

    common = dict(blank)
    ambiguous_rank = ""
    for offset, rank in enumerate(RANKS[1:]):
        if RANKS.index(rank) > deepest_supported_index:
            continue
        values = {catalog.lineages[index][offset] for index in best_indices}
        if len(values) == 1:
            common[rank] = next(iter(values))
        elif not ambiguous_rank:
            ambiguous_rank = rank

    strict = dict(common)
    conflict_rank = ""
    conflicting_taxa: list[str] = []
    for rank in RANKS[1:]:
        selected = strict.get(rank, "")
        if not selected:
            continue
        alternatives: dict[str, set[str]] = defaultdict(set)
        for claim in valid_claims:
            taxon = str(claim["lineage"].get(rank, ""))
            if taxon and taxon != selected:
                alternatives[taxon].add(str(claim["method_family"]))
        supported_alternatives = sorted(
            taxon for taxon, families in alternatives.items()
            if len(families) >= 2
        )
        if supported_alternatives:
            conflict_rank = rank
            conflicting_taxa = [selected, *supported_alternatives]
            for later_rank in RANKS[RANKS.index(rank):]:
                strict[later_rank] = ""
            break

    selected_species = common.get("species", "")
    species_supporters = [
        claim for claim in valid_claims
        if selected_species and claim["lineage"].get("species") == selected_species
    ]
    species_families = {
        str(claim["method_family"]) for claim in species_supporters
    }
    species_max_tier = max(
        (TIER_SCORE[str(claim["tier"])] for claim in species_supporters),
        default=0,
    )
    named_species = sorted({
        str(claim["lineage"].get("species", ""))
        for claim in valid_claims if claim["lineage"].get("species")
    })
    alternative_species = [
        taxon for taxon in named_species if taxon != selected_species
    ]
    species_status = "unresolved"
    species_basis = "insufficient_species_support"
    if selected_species and len(species_families) >= 2 and not conflict_rank:
        species_status = "strict_consensus"
        species_basis = "two_or_more_independent_method_families"
    elif (
        selected_species and species_max_tier == TIER_SCORE["high"]
        and not alternative_species and not conflict_rank
    ):
        species_status = "strict_single_high"
        species_basis = "single_high_call_without_named_alternative"
    else:
        strict["species"] = ""
        if conflict_rank and conflict_rank != "species":
            species_status = "withheld_lineage_conflict"
            species_basis = "ancestor_lineage_conflict"
        elif not selected_species and ambiguous_rank == "species":
            species_status = "ambiguous_msl_species"
            species_basis = "equally_supported_msl_species"
        elif not selected_species and ambiguous_rank:
            species_status = "withheld_lineage_conflict"
            species_basis = "ambiguous_ancestor_lineage"
        elif alternative_species or conflict_rank == "species":
            species_status = "withheld_species_conflict"
            species_basis = "named_alternative_species"
        elif selected_species and species_max_tier == TIER_SCORE["moderate"]:
            species_status = "exploratory_single_moderate"
            species_basis = "single_moderate_species_call"
        elif selected_species:
            species_status = "exploratory_single_provisional"
            species_basis = "single_provisional_species_call"

    strict_rank = next(
        (rank for rank in reversed(RANKS) if strict.get(rank)), "domain"
    )
    exploratory = dict(common)
    exploratory_rank = next(
        (rank for rank in reversed(RANKS) if exploratory.get(rank)), "domain"
    )

    def support_at(
        lineage: dict[str, str], rank: str
    ) -> tuple[list[str], list[str], int]:
        taxon = lineage.get(rank, "")
        supporting = [
            claim for claim in valid_claims
            if taxon and claim["lineage"].get(rank) == taxon
        ]
        return (
            sorted({str(claim["tool"]) for claim in supporting}),
            sorted({str(claim["method_family"]) for claim in supporting}),
            max(
                (TIER_SCORE[str(claim["tier"])] for claim in supporting),
                default=0,
            ),
        )

    strict_tools, strict_families, strict_tier = support_at(strict, strict_rank)
    exploratory_tools, exploratory_families, exploratory_tier = support_at(
        exploratory, exploratory_rank
    )
    if len(strict_families) >= 2:
        strict_confidence = "multi_method_consensus"
    elif strict_tier == TIER_SCORE["high"]:
        strict_confidence = "single_method_high"
    elif strict_rank == "domain":
        strict_confidence = "unclassified"
    else:
        strict_confidence = "moderate_supported"
    exploratory_confidence = {
        3: "high", 2: "moderate", 1: "provisional", 0: "ambiguous",
    }[exploratory_tier]

    if conflict_rank:
        stop_reason = "independent_method_conflict"
    elif strict_rank != exploratory_rank:
        stop_reason = "insufficient_rank_support"
    elif ambiguous_rank:
        stop_reason = "ambiguous_msl_paths"
    else:
        stop_reason = "deepest_supported_rank"
    if taxonomy_string(strict) == taxonomy_string(exploratory):
        exploratory_status = "matches_strict"
    elif conflict_rank:
        exploratory_status = "alternative_at_conflict"
    elif strict_rank == "domain":
        exploratory_status = "exploratory_only"
    else:
        exploratory_status = "extends_strict"

    for vote in votes:
        if not vote["valid"]:
            continue
        rank = str(vote["rank"])
        taxon = str(vote["taxon"])
        if strict.get(rank) == taxon:
            vote["accepted"] = True
            vote["reason"] = "selected_strict"
        elif exploratory.get(rank) == taxon:
            vote["reason"] = "selected_exploratory"
        elif rank == conflict_rank:
            vote["reason"] = "independent_method_conflict"
        else:
            vote["reason"] = "lower_global_path_support"

    return TaxonomyDecision(
        strict, exploratory, strict_rank, exploratory_rank,
        strict_confidence, exploratory_confidence, strict_tools,
        exploratory_tools, strict_families, exploratory_families,
        stop_reason, conflict_rank or "NA", conflicting_taxa,
        exploratory_status, species_status, species_basis,
        alternative_species, votes,
    )


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


@lru_cache(maxsize=1)
def load_reconciliation(path: str, expected_sha256: str) -> Reconciler:
    """Load once per process, keyed by the exact bundle identity, never online."""
    bundle = Path(path)
    if sha256(bundle) != expected_sha256:
        raise ValueError("Pinned reconciliation bundle changed before loading")
    return Reconciler.load(bundle)


def decide_candidate_taxonomy(rows, msl, minimum_length, context, audit):
    """The existing taxonomy engine gets copies; retention sees original rows."""
    if context.reconciliation_bundle is None:
        return decide_taxonomy(rows, msl, minimum_length)
    resolver = load_reconciliation(context.reconciliation_bundle, context.reconciliation_sha256)
    copied_rows = []
    for row in rows:
        tool = row.get("tool", "").strip().lower()
        if tool not in TAXONOMY_TOOLS or row.get("classification") != "virus":
            copied_rows.append(row)
            continue
        if tool == "vcontact3":
            try:
                length = int(float(row.get("length", "0") or 0))
            except ValueError:
                length = 0
            if length < minimum_length:
                copied_rows.append(row)
                continue
        original = {rank: clean_taxon(row.get(RANK_COLUMNS[rank], ""), rank) for rank in RANKS[1:]}
        resolution = resolver.resolve(original)
        copied = dict(row)
        if resolution.normalized:
            for rank in RANKS[1:]:
                copied[RANK_COLUMNS[rank]] = PREFIXES[rank] + resolution.normalized.get(rank, "unclassified")
        copied_rows.append(copied)
        if resolution.status != "no_named_taxonomy" and (context.full_audit or resolution.status != "unchanged"):
            audit.append({"sample_id": context.sample_id, "tool": tool,
                "source_sequence_id": row.get("sequence_id", ""),
                "source_file": row.get("_source_file", ""),
                "database_domain": row.get("vcontact3_database_domain", ""),
                "status": resolution.status, "target_release": resolution.target_release,
                "source_release": resolution.source_release,
                "original_taxonomy": json.dumps(resolution.original, sort_keys=True),
                "normalized_taxonomy": json.dumps(resolution.normalized, sort_keys=True),
                "unresolved": json.dumps(resolution.unresolved, sort_keys=True),
                "mappings": json.dumps(resolution.mappings, sort_keys=True)})
    return decide_taxonomy(copied_rows, msl, minimum_length)


@contextmanager
def reconciliation_audit_sink(path):
    """Stream provenance in the ordered parent, including in compact mode."""
    if path is None:
        yield None
        return
    columns = ["sample_id", "final_sequence_id", "tool", "source_sequence_id",
               "source_file", "database_domain", "status", "target_release",
               "source_release", "original_taxonomy", "normalized_taxonomy", "unresolved", "mappings"]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        yield writer


@contextmanager
def candidate_routing_sinks(prefix):
    """Stream audits/summaries; do not retain another full-dataset row collection."""
    with ExitStack() as stack:
        audit = stack.enter_context(gzip.open(
            Path(f"{prefix}.candidate_routing_audit.tsv.gz"), "wt", encoding="utf-8", newline=""))
        writers = {"all": csv.DictWriter(audit, fieldnames=te_routing.AUDIT_COLUMNS,
                                        delimiter="\t", lineterminator="\n")}
        for destination in ("TE_candidates", "viral_candidates"):
            handle = stack.enter_context(Path(f"{prefix}.{destination}.summary.tsv").open(
                "w", encoding="utf-8", newline=""))
            writers[destination] = csv.DictWriter(handle, fieldnames=te_routing.SUMMARY_COLUMNS,
                                                   delimiter="\t", lineterminator="\n")
        for writer in writers.values():
            writer.writeheader()
        yield writers


def evaluate_candidate(task: CandidateTask, context: CandidateContext, msl: MslCatalog) -> CandidateResult:
    """Evaluate one complete candidate; decision rules are unchanged."""
    final_id = task.final_id
    minimum_vcontact3_length = context.minimum_vcontact3_length
    pair_floor = context.pair_floor
    evidence_audit_rows = []
    taxonomy_audit_rows = []
    reconciliation_audit_rows = []
    review = None
    region = task.region
    record_type = region["record_type"].strip()
    parent_id = task.parent_id
    header = task.header
    all_rows = task.evidence
    rows, region_vicat_rows, suppress_parent_vicat = select_vicat_scope(
        all_rows, final_id, parent_id
    )
    rows, parent_context_rows = partition_refined_region_context(
        rows, final_id, record_type
    )
    parent_context_ids = {id(row) for row in parent_context_rows}
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
        if id(source_row) in parent_context_ids:
            audit["evidence_application"] = "parent_context_only"
        if is_review_only_evidence(source_row):
            audit["evidence_application"] = "review_only_no_consensus_vote"
        evidence_audit_rows.append(audit)
    tools_by_class: dict[str, set[str]] = defaultdict(set)
    strong_tools: set[str] = set()
    qualified_tools: set[str] = set()
    checkv_complete_support = checkv_supports_complete_viral_contig(
        rows, record_type
    )
    for row in rows:
        tool = row.get("tool", "").strip().lower()
        if is_review_only_evidence(row):
            continue
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

    taxonomy = decide_candidate_taxonomy(
        rows, msl, minimum_vcontact3_length, context, reconciliation_audit_rows
    )
    for vote in taxonomy.votes:
        vote["sample_id"] = context.sample_id
        vote["final_sequence_id"] = final_id
        taxonomy_audit_rows.append(vote)

    _, group_labels, group_statuses, group_audit = apply_vcontact3_groups(
        task.groups,
        taxonomy.strict,
        task.length,
        minimum_vcontact3_length,
        context.sample_id,
    )
    for vote in group_audit:
        vote["sample_id"] = context.sample_id
        vote["final_sequence_id"] = final_id
        taxonomy_audit_rows.append(vote)

    coords = region["coordinates"].strip() or "NA"
    discovery_status = task.discovery_status.strip()
    strict_taxonomy = taxonomy_string(taxonomy.strict)
    exploratory_taxonomy = taxonomy_string(taxonomy.exploratory)
    original_name = header["original_id"]
    cellular = sorted(tools_by_class.get("cellular", set()))
    plasmid = sorted(tools_by_class.get("plasmid", set()))
    parent_cellular_context = sorted({
        row.get("tool", "").strip().lower()
        for row in parent_context_rows
        if row.get("classification", "").strip().lower() == "cellular"
        and row.get("tool", "").strip()
    })
    parent_plasmid_context = sorted({
        row.get("tool", "").strip().lower()
        for row in parent_context_rows
        if row.get("classification", "").strip().lower() == "plasmid"
        and row.get("tool", "").strip()
    })
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
    has_vicat_plasmid_loci = any(
        item.upper().split(":", 1)[0] == "PLASMID"
        for item in vicat_nonviral_classes
    )
    reported_plasmid_conflicts = sorted(
        set(plasmid) | ({"vicat"} if has_vicat_plasmid_loci else set())
    )
    if "resolved_provirus_with_vicat_support" in vicat_patterns:
        vicat_provirus_status = "resolved_with_vicat_support"
    elif "localized_viral_cluster" in vicat_patterns:
        vicat_provirus_status = "potential_provirus"
    else:
        vicat_provirus_status = "none"
    interpretation = tesorter_summary["sequence_interpretation"]
    previous_decision = adjudicate_viral_decision(
        discovery_status, record_type, strong_tools, qualified_tools, cellular, plasmid, interpretation)
    viral_decision = adjudicate_viral_decision(
        discovery_status, record_type, strong_tools, qualified_tools, cellular, plasmid,
        "viral_candidate" if interpretation in te_routing.TE_VETO else interpretation)
    boundary_primary_eligible, provirus_retention_basis = adjudicate_disputed_ct3_boundary(
        region, rows, selected_vicat_rows)
    floor_qualified_tools = rna_floor_qualified_tools(rows, qualified_tools)
    if pair_floor and context.input_type == "rna" and not context.has_locus_evidence and not strong_tools and floor_qualified_tools == {"deep6", "vicat"}:
        raise ValueError("RNA pair homology floor requires viCAT locus evidence")

    def preserve_origin_safeguards(decision):
        exception_loci = []
        if not context.disable_rna_checkv_exception:
            # Keep the ACTUAL mobile-element interpretation in this exception.
            exception_loci = rna_checkv_exception_loci(
                context.input_type, decision, record_type, cellular, reported_plasmid_conflicts,
                strong_tools | qualified_tools, interpretation, task.loci, parent_id, task.length)
        if exception_loci:
            decision = "retained_viral"
        if decision == "retained_provirus" and not boundary_primary_eligible:
            decision = "provisional_provirus"
        decision, floor_status = apply_rna_pair_homology_floor(
            decision, pair_floor, context.input_type, strong_tools, floor_qualified_tools,
            task.loci, parent_id, coords if record_type != "input_contig" else "")
        return decision, floor_status, exception_loci

    ordinary_input = viral_decision
    viral_decision, pair_floor_status, checkv_exception_loci = preserve_origin_safeguards(ordinary_input)
    previous_decision, previous_floor_status, _ = (
        preserve_origin_safeguards(previous_decision) if previous_decision != ordinary_input
        else (viral_decision, pair_floor_status, checkv_exception_loci))
    previous_confidence = "provisional" if previous_floor_status == "demoted" else viral_confidence
    if pair_floor_status == "demoted":
        viral_confidence = "provisional"
    entity_interpretation = viral_entity_interpretation(
        viral_decision,
        record_type,
        parent_plasmid_context,
        reported_plasmid_conflicts,
        has_vicat_plasmid_loci,
        strong_tools,
    )
    generic_cellular_conflict_overridden = (
        viral_decision in PRIMARY_VIRAL_DECISIONS
        and len(strong_tools) >= 2
        and bool(cellular)
        and set(cellular).issubset(GENERIC_ORIGIN_CLASSIFIERS)
        and not plasmid
    )
    metadata = {
        "original_contig_name": original_name,
        "rna_checkv_conflict_exception": str(bool(checkv_exception_loci)).lower(),
        "rna_checkv_exception_supporting_loci": ",".join(checkv_exception_loci) or "NA",
        "normalized_name": parent_id,
        "final_sequence_id": final_id,
        "record_type": record_type,
        "provirus_coordinates": coords,
        "original_length": region["original_length"],
        "refined_length": task.length,
        "discovery_status": discovery_status,
        "viral_decision": viral_decision,
        "viral_confidence": viral_confidence,
        "strong_tools": ",".join(sorted(strong_tools)) or "NA",
        "qualified_tools": ",".join(sorted(qualified_tools)) or "NA",
        "cellular_conflict_tools": ",".join(cellular) or "NA",
        "parent_cellular_context_tools": (
            ",".join(parent_cellular_context) or "NA"
        ),
        "parent_plasmid_context_tools": (
            ",".join(parent_plasmid_context) or "NA"
        ),
        "checkv_complete_viral_contig_support": str(
            checkv_complete_support
        ).lower(),
        "generic_cellular_conflict_overridden": str(
            generic_cellular_conflict_overridden
        ).lower(),
        "plasmid_conflict_tools": (
            ",".join(reported_plasmid_conflicts) or "NA"
        ),
        "strict_taxonomy": strict_taxonomy,
        "strict_taxonomy_rank": taxonomy.strict_rank,
        "strict_taxonomy_confidence": taxonomy.strict_confidence,
        "strict_taxonomy_supporting_tools": (
            ",".join(taxonomy.strict_tools) or "NA"
        ),
        "strict_taxonomy_supporting_method_families": (
            ",".join(taxonomy.strict_families) or "NA"
        ),
        "strict_taxonomy_stop_reason": taxonomy.strict_stop_reason,
        "strict_taxonomy_conflict_rank": taxonomy.conflict_rank,
        "strict_taxonomy_conflicting_taxa": (
            ",".join(taxonomy.conflicting_taxa) or "NA"
        ),
        "exploratory_taxonomy": exploratory_taxonomy,
        "exploratory_taxonomy_rank": taxonomy.exploratory_rank,
        "exploratory_taxonomy_confidence": taxonomy.exploratory_confidence,
        "exploratory_taxonomy_supporting_tools": (
            ",".join(taxonomy.exploratory_tools) or "NA"
        ),
        "exploratory_taxonomy_supporting_method_families": (
            ",".join(taxonomy.exploratory_families) or "NA"
        ),
        "exploratory_taxonomy_status": taxonomy.exploratory_status,
        "strict_species": taxonomy.strict.get("species", "") or "NA",
        "exploratory_species": (
            taxonomy.exploratory.get("species", "") or "NA"
        ),
        "species_assignment_status": taxonomy.species_status,
        "species_assignment_basis": taxonomy.species_basis,
        "species_conflicting_taxa": (
            ",".join(taxonomy.species_conflicting_taxa) or "NA"
        ),
        "taxonomy_conflict": str(taxonomy.conflict_rank != "NA").lower(),
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
        "provirus_boundary_source": (
            region.get("boundary_source", "").strip() or "NA"
        ),
        "provirus_boundary_status": (
            region.get("boundary_status", "").strip() or "NA"
        ),
        "provirus_boundary_supporting_tools": (
            region.get("supporting_boundary_tools", "").strip() or "NA"
        ),
        "provirus_retention_basis": provirus_retention_basis,
        "viral_entity_interpretation": entity_interpretation,
        **tesorter_summary,
    }
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
    if not boundary_primary_eligible:
        reasons.append("disputed_ct3_vicat_boundary_provisional")
    if viral_confidence == "provisional": reasons.append("single_qualified_viral_tool")
    if generic_cellular_conflict_overridden:
        reasons.append("generic_cellular_conflict_overridden")
    elif cellular:
        reasons.append("cellular_conflict")
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
    if taxonomy.strict_rank == "domain": reasons.append("taxonomy_unclassified")
    if reasons:
        review = {"final_sequence_id": final_id, "review_reasons": ",".join(reasons), **metadata}

    reporting_audit = apply_taxonomy_reporting_cap(metadata, taxonomy, context.max_rank, context.sample_id)
    routing_audit = te_routing.describe_candidate(
        metadata, evidence_audit_rows, reporting_audit, taxonomy_audit_rows,
        context.sample_id, previous_decision, previous_confidence)
    metadata.update({name: routing_audit[name] for name in te_routing.METADATA_COLUMNS})
    if routing_audit["candidate_destination"] == "TE_candidates":
        metadata["viral_decision"] = "TE_candidate"
    if review is not None:
        review.update(metadata)

    # Compact/none never emit these audits. Do not ship unused audit copies
    # back through the process queue or retain them in the parent.
    return CandidateResult(
        metadata, review,
        evidence_audit_rows if context.full_audit else [],
        taxonomy_audit_rows if context.full_audit else [],
        pair_floor_status, bool(checkv_exception_loci),
        reconciliation_audit_rows,
        reporting_audit,
        routing_audit,
    )


def run(args: argparse.Namespace) -> None:
    workers = getattr(args, "workers", 1)
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError("viHARMONY workers must be a positive integer")
    max_rank = getattr(args, "max_rank", "family")
    if not isinstance(max_rank, str) or max_rank not in RANKS:
        raise ValueError("viHARMONY maximum reporting rank must be one of: " + ", ".join(RANKS))
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
    evidence_tables: list[EvidenceTable] = []
    tools_present: set[str] = set()
    evidence_columns: list[str] = []
    pair_loci_by_parent = EvidenceIndex()
    pair_floor = getattr(args, "rna_pair_homology_floor", False)
    for path in args.evidence:
        locus_file = path.name.endswith(".vicat_orf_evidence.tsv")
        rows = EvidenceTable(path, source_file=None if locus_file else path.name)
        columns = rows.fieldnames
        if locus_file:
            require(columns, {"sample_id", "sequence_id", "coordinates", "locus_classification",
                              "protein_length", "best_bitscore", "best_query_coverage", "best_subject_coverage"}, path)
            if any(row["sample_id"] != args.sample_id for row in rows):
                raise ValueError(f"Locus sample mismatch in {path}")
            for row in rows:
                # Keep exact ID matching and original locus order.
                pair_loci_by_parent.add(row["sequence_id"], row)
            continue
        require(columns, {"sample_id", "sequence_id", "parent_sequence_id", "tool", "classification"}, path)
        evidence_columns.extend(column for column in columns if column not in evidence_columns)
        for row in rows:
            if row["sample_id"] != args.sample_id:
                raise ValueError(f"Evidence sample mismatch in {path}: {row['sample_id']}")
            if row.get("tool"):
                tools_present.add(row["tool"])
        evidence_tables.append(rows)

    groups_by_sequence = EvidenceIndex()
    for path in args.vcontact3_groups:
        rows = EvidenceTable(path)
        columns = rows.fieldnames
        if "sequence_id" not in columns:
            continue
        for row in rows:
            groups_by_sequence.add(row["sequence_id"], row)

    msl, msl_release = load_msl(args.ictv_msl)
    reconciliation_bundle = getattr(args, "ictv_reconciliation", None)
    reconciliation_identity = None
    reconciliation_metadata = None
    if reconciliation_bundle is not None:
        reconciliation_bundle = reconciliation_bundle.resolve()
        reconciliation_identity = sha256(reconciliation_bundle)
        # Revalidate the manifest on every run, even when a previous sample in
        # this Python process populated the per-worker resolver cache.
        resolver = Reconciler.load(reconciliation_bundle)
        target_lineages = {tuple(resolver.lineages[nid].get(rank, "") for rank in RANKS[1:])
                           for nid, node in resolver.nodes.items() if node["rank"] == "species"}
        if target_lineages != set(msl.lineages):
            raise ValueError("Reconciliation bundle target differs from configured ICTV MSL")
        reconciliation_metadata = {"bundle": reconciliation_bundle.name,
            "sha256": reconciliation_identity, "target_release": resolver.target_release,
            "source_commit": resolver.bundle["source_commit"],
            "manifest_sha256": sha256(reconciliation_bundle.with_name(reconciliation_bundle.name + ".manifest.json")),
            "scope": "taxonomy-only copies; original evidence used for viral retention",
            "source_release_handling": "unknown per-tool release; exact current names retain current meaning",
            "audit_scope": "all named eligible calls" if args.audit_mode == "full" else "changed or unresolved named eligible calls"}
    rows_by_final = EvidenceIndex.from_pairs(iter_evidence_by_final(
        chain.from_iterable(evidence_tables), final_to_parent, regions, refined
    ))
    # Indexes own the tables still needed for decisions/reports. No flat list
    # of row views or per-row dictionaries remains resident in the parent.
    del evidence_tables

    metadata_rows: list[dict[str, object]] = []
    pair_floor_demoted_ids = []
    rna_checkv_exception_ids = []
    review_rows: list[dict[str, object]] = []
    evidence_audit_rows: list[dict[str, object]] = []
    taxonomy_audit_rows: list[dict[str, object]] = []
    normalized_fasta_rows: list[tuple[str, str]] = []
    original_fasta_rows: list[tuple[str, str]] = []
    provisional_normalized_rows: list[tuple[str, str]] = []
    provisional_original_rows: list[tuple[str, str]] = []
    used_original_output_ids: set[str] = set()
    map_rows: list[dict[str, object]] = []

    context = CandidateContext(
        args.sample_id, args.input_type, minimum_vcontact3_length, pair_floor,
        any(path.name.endswith(".vicat_orf_evidence.tsv") for path in args.evidence),
        getattr(args, "disable_rna_checkv_exception", False), args.audit_mode == "full",
        str(reconciliation_bundle) if reconciliation_bundle is not None else None,
        reconciliation_identity,
        max_rank,
    )
    tasks = (
        CandidateTask(
            final_id, len(sequence), regions[final_id], final_to_parent[final_id],
            headers[final_to_parent[final_id]], gate[final_to_parent[final_id]]["discovery_status"],
            rows_by_final.get(final_id, []), groups_by_sequence.get(final_id, []),
            pair_loci_by_parent.get(final_to_parent[final_id], ()),
        )
        for final_id, sequence in refined.items()
    )
    active_workers = min(workers, len(refined)) or 1
    batch_size = min(CANDIDATE_BATCH_SIZE, max(1, (len(refined) + active_workers * 4 - 1) // (active_workers * 4)))
    reconciliation_path = Path(f"{args.output_prefix}.taxonomy_reconciliation_audit.tsv.gz") if reconciliation_metadata else None
    reporting_path = Path(f"{args.output_prefix}.taxonomy_reporting_audit.tsv.gz")
    with gzip.open(reporting_path, "wt", encoding="utf-8", newline="") as reporting_handle, reconciliation_audit_sink(reconciliation_path) as reconciliation_writer, candidate_routing_sinks(args.output_prefix) as routing_writers, candidate_results(tasks, context, msl, active_workers, batch_size) as results:
        reporting_writer = csv.DictWriter(reporting_handle, fieldnames=REPORTING_AUDIT_COLUMNS, delimiter="\t", lineterminator="\n")
        reporting_writer.writeheader()
        for final_id, result in zip(refined, results, strict=True):
            metadata = result.metadata
            if metadata["final_sequence_id"] != final_id:
                raise RuntimeError("viHARMONY worker returned an out-of-order candidate")
            metadata_rows.append(metadata)
            routing_writers["all"].writerow(result.routing_audit)
            destination = result.routing_audit["candidate_destination"]
            if destination in ("TE_candidates", "viral_candidates"):
                routing_writers[destination].writerow({
                    column: result.routing_audit[column] for column in te_routing.SUMMARY_COLUMNS})
            reporting_writer.writerow(result.reporting_audit)
            if result.review is not None:
                review_rows.append(result.review)
            evidence_audit_rows.extend(result.evidence_audit)
            taxonomy_audit_rows.extend(result.taxonomy_audit)
            if reconciliation_writer is not None:
                for audit in result.reconciliation_audit:
                    reconciliation_writer.writerow({**audit, "final_sequence_id": final_id})
            if result.pair_floor_status == "demoted":
                pair_floor_demoted_ids.append(final_id)
            if result.rna_checkv_exception:
                rna_checkv_exception_ids.append(final_id)
            sequence = refined[final_id]
            parent_id = final_to_parent[final_id]
            header = headers[parent_id]
            original_name = metadata["original_contig_name"]
            record_type = metadata["record_type"]
            coords = metadata["provirus_coordinates"]
            viral_decision = metadata["viral_decision"]
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
        "taxonomy_reporting_audit": reporting_path,
        "TE_candidates_fasta": Path(f"{prefix}.TE_candidates.fasta"),
        "TE_candidates_summary": Path(f"{prefix}.TE_candidates.summary.tsv"),
        "viral_candidates_fasta": Path(f"{prefix}.viral_candidates.fasta"),
        "viral_candidates_summary": Path(f"{prefix}.viral_candidates.summary.tsv"),
        "candidate_routing_audit": Path(f"{prefix}.candidate_routing_audit.tsv.gz"),
    }
    for destination, label in (("TE_candidates", "TE_candidates"), ("viral_candidates", "viral_candidates")):
        selected_ids = {row["final_sequence_id"] for row in metadata_rows
                        if row["candidate_destination"] == destination}
        write_fasta(outputs[label + "_fasta"],
                    [(sid, sequence) for sid, sequence in refined.items() if sid in selected_ids])
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

    final_by_parent: dict[str, list[str]] = defaultdict(list)
    provisional_by_parent: dict[str, list[str]] = defaultdict(list)
    decisions_by_parent: dict[str, list[str]] = defaultdict(list)
    decision_by_final = {
        str(row["final_sequence_id"]): str(row["viral_decision"])
        for row in metadata_rows
    }
    for row in map_rows:
        final_id = str(row["final_sequence_id"])
        decisions_by_parent[str(row["normalized_parent_id"])].append(
            decision_by_final.get(final_id, "")
        )
        if decision_by_final.get(final_id) in PRIMARY_VIRAL_DECISIONS:
            final_by_parent[str(row["normalized_parent_id"])].append(final_id)
        elif decision_by_final.get(final_id) in PROVISIONAL_VIRAL_DECISIONS:
            provisional_by_parent[str(row["normalized_parent_id"])].append(final_id)
    TE_by_parent = defaultdict(list)
    for row in metadata_rows:
        if row["candidate_destination"] == "TE_candidates":
            TE_by_parent[row["normalized_name"]].append(row["final_sequence_id"])
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
        elif TE_by_parent.get(sequence_id):
            disposition = "TE_candidates_after_harmony"
        elif gate_row["advance_to_refinement"].strip().lower() == "true":
            candidate_decisions = decisions_by_parent.get(sequence_id, [])
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

    if reconciliation_metadata is not None:
        if sha256(reconciliation_bundle) != reconciliation_identity or sha256(reconciliation_bundle.with_name(reconciliation_bundle.name + ".manifest.json")) != reconciliation_metadata["manifest_sha256"]:
            raise ValueError("Pinned reconciliation bundle or manifest changed during the run")
        outputs["taxonomy_reconciliation_audit"] = reconciliation_path
    outputs = reporting_bundle.finalize(
        prefix, args.sample_id, metadata_rows, review_rows,
        rows_by_final, groups_by_sequence, outputs,
    )
    manifest_path = prefix.parent / "audit" / f"{prefix.name}.harmonizer_manifest.json"
    manifest = {
        "schema_version": "viharmony-1.0",
        "sample_id": args.sample_id,
        "input_type": args.input_type,
        "ictv_msl": args.ictv_msl.name,
        "ictv_msl_release": msl_release,
        "audit_mode": args.audit_mode,
        "reporting_bundle": {"schema_version": reporting_bundle.VERSION, "output_root": "..",
                             "missing_values": "NA", "review_fasta": False},
        "tools_present": sorted(tools_present),
        "policy": {
            "rna_checkv_conflict_exception": bool(args.input_type == "rna" and not getattr(args, "disable_rna_checkv_exception", False)),
            "rna_checkv_exception_ids": sorted(rna_checkv_exception_ids),
            "rna_checkv_exception_thresholds": "RNA intact only; sole CheckV cellular conflict; no plasmid/mobile conflict; Deep6+VirBot+viCAT qualified/strong; no cellular/ambiguous loci; viral protein>=100aa; bitscore>=100; identity>=50%; query and subject coverage>=80%; coding span>=70% of input; preserve CheckV annotation",
            "rna_pair_homology_floor": bool(pair_floor and args.input_type == "rna"),
            "virsorter2_hallmark_only_policy": "review-only; no discovery routing or primary-consensus votes; no invented scores/boundaries; does not exempt RNA homology floor",
            "rna_pair_homology_floor_demoted_ids": sorted(pair_floor_demoted_ids),
            "rna_pair_homology_floor_thresholds": "protein_length>=100; bitscore>=100; query_coverage>=50; subject_coverage>=50; one viral-supported locus; qualified Deep6+viCAT without strong support after ignoring hallmark-only VirSorter2 for exemption eligibility",
            "viral_confidence": "high=2+ strong; supported=1 strong or 2+ qualified; provisional=1 qualified",
            "checkv_complete_contig": "complete/high-quality, provirus=No, viral genes present, and zero host genes is strong intact-contig viral support",
            "generic_cellular_override": "2+ strong independent viral tools override Deep6/DeepMicroClass2-only cellular conflict; the conflict remains annotated",
            "taxonomy": "global configured-ICTV-MSL path resolution with conservative strict and explicit exploratory outputs",
            "harmonizer_max_rank": max_rank,
            "taxonomy_reporting_cap": "reporting only, after inference/retention/review; uncapped consensus preserved in taxonomy_reporting_audit in every audit mode; conflicts remain flagged",
            "strict_species": "2+ independent method families, or one high-tier species call with no named alternative",
            "exploratory_species": "single moderate/provisional calls remain exploratory and never populate strict species",
            "vcontact3_min_taxonomy_length": minimum_vcontact3_length,
            "vcontact3_project_groups": "ancestry-compatible, unambiguous groups only",
            "vicat_scope": "refined-region evidence supersedes parent-discovery evidence",
            "provirus_adjudication": "refinement selects coordinates; disputed CT3/viCAT boundaries require strong multi-hallmark CT3 evidence or at least 2 viCAT viral loci with regional viral fraction >=0.60; parent-only cellular/plasmid calls are retained as context",
            "plasmid_viral_entities": "parent plasmid evidence annotates extracted proviruses without vetoing them; strong intact-contig viral consensus plus exact plasmid evidence is retained as a hybrid candidate",
            "tesorter_retrovirus": "TE-routing-v1: strong coherent complete TE with compatible uncapped consensus exports to TE candidates; partial/mixed TE annotates ordinary eligibility; actual TE interpretation remains in RNA exception",
            "TE_routing_policy": te_routing.POLICY_VERSION,
            "viral_evidence_scope": "localized clusters, single/multiple loci, and record-level predictions are distinct; no whole-sequence viral-origin assertion",
            "final_output": "primary requires strong or multi-tool qualified viral support; single-tool qualified calls are emitted separately as provisional",
        },
        "outputs": {},
    }
    if reconciliation_metadata is not None:
        manifest["taxonomy_reconciliation"] = reconciliation_metadata
    for label, path in outputs.items():
        manifest["outputs"][label] = {"file": path.relative_to(prefix.parent).as_posix(), "sha256": sha256(path)}
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"viHARMONY sample={args.sample_id} workers={active_workers} candidates={len(metadata_rows)} "
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
