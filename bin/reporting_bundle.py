"""Reporting-only projections of completed viHARMONY decisions and evidence.

Never adjudicates candidates or edits evidence, metadata, sequences or taxonomy.
Reports stream in final-record order; raw tools remain available in the audits.
"""
from __future__ import annotations

import csv
import re
from contextlib import ExitStack
from pathlib import Path

VERSION = "reporting-bundle-v1"
RANKS = ("domain", "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species")
GROUP_RANKS = ("realm", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus")
SUMMARY_COLUMNS = [
    "sample_id", "final_sequence_id", "original_sequence_id", "normalized_parent_id",
    "record_type", "coordinates", "length_nt", "candidate_status", "viral_confidence",
    "strong_tools", "qualified_tools", "conflict_tools", "reported_taxonomy",
    "reported_taxonomy_rank", "taxonomy_support_status", "taxonomy_supporting_tools",
    "viral_evidence_scope", "TE_flag", "TE_architecture_status", "TE_evidence_scope",
    "tesorter_orders", "tesorter_superfamilies", "checkv_quality", "checkv_complete_call",
    "vcontact3_group_status", "review_flags",
]
TAXONOMY_COLUMNS = ["final_sequence_id", *RANKS]
QUALITY_COLUMNS = [
    "sample_id", "final_sequence_id", "checkv_sequence_id", "assessment_status",
    "assessment_scope", "assessment_coordinates", "checkv_record_length_nt",
    "checkv_quality", "miuvig_quality", "checkv_complete_call", "completeness_percent",
    "completeness_method", "contamination_percent", "provirus", "proviral_length_nt",
    "viral_genes", "host_genes", "complete_genome_prediction", "complete_genome_confidence",
    "topology", "aai_confidence", "hmm_completeness_lower_percent",
    "hmm_completeness_upper_percent", "warnings",
]
GROUP_COLUMNS = [
    "sample_id", "final_sequence_id", "database_domain", "database_version",
    "assignment_status", "rank", "raw_prediction", "prediction_type", "evidence_type",
    "network_support", "vog_support", "standardization_decision", "host_domain", "vog_host_domain",
]
REVIEW_COLUMNS = [*SUMMARY_COLUMNS, "review_category", "review_reason"]
QUALITY_MAP = {
    "checkv_record_length_nt": "length", "completeness_percent": "completeness",
    "contamination_percent": "contamination", "proviral_length_nt": "proviral_length",
    "hmm_completeness_lower_percent": "hmm_completeness_lower",
    "hmm_completeness_upper_percent": "hmm_completeness_upper",
}
LEGACY_LABELS = {
    "normalized_fasta", "original_fasta", "provisional_normalized_fasta",
    "provisional_original_fasta", "provisional_metadata", "database_fasta", "database_tsv",
    "all_candidates_fasta",
}


def cell(value):
    # Literal upstream categories such as "Not-determined" are not erased.
    return "NA" if value is None or str(value).strip() == "" else str(value)


def writer(stack, path, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = stack.enter_context(path.open("w", encoding="utf-8", newline=""))
    result = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
    result.writeheader()
    return result


def complete_call(quality, scope):
    if scope != "exact_candidate":
        return "NA"
    normalized = str(quality).strip().lower()
    if normalized == "complete":
        return "true"
    if normalized in {"high-quality", "medium-quality", "low-quality"}:
        return "false"
    return "NA"


def quality_rows(sample, metadata, evidence):
    final = metadata["final_sequence_id"]
    found = False
    for source in evidence:
        if source.get("tool", "").strip().lower() != "checkv":
            continue
        found = True
        source_id = source.get("sequence_id", "")
        # The current standardizer copies aggregate quality onto proviral rows.
        # An exact child ID therefore does NOT establish child-specific quality.
        if source.get("parent_sequence_id", "").strip() or source.get("record_type") == "provirus":
            scope = "proviral_aggregate_context"
        elif source_id == final and metadata["record_type"] == "input_contig":
            scope = "exact_candidate"
        else:
            scope = "parent_context"
        result = {key: cell(source.get(QUALITY_MAP.get(key, key))) for key in QUALITY_COLUMNS}
        result.update(sample_id=sample, final_sequence_id=final, checkv_sequence_id=cell(source_id),
                      assessment_status="available", assessment_scope=scope,
                      assessment_coordinates=(cell(source.get("coordinates")) if scope == "exact_candidate" else "NA"),
                      checkv_complete_call=complete_call(source.get("checkv_quality"), scope))
        yield result
    if not found:
        yield {**dict.fromkeys(QUALITY_COLUMNS, "NA"), "sample_id": sample,
               "final_sequence_id": final, "assessment_status": "no_matching_record"}


def group_rows(sample, final, sources):
    if not sources:
        yield {**dict.fromkeys(GROUP_COLUMNS, "NA"), "sample_id": sample,
               "final_sequence_id": final, "assignment_status": "no_matching_record"}
    for source in sources:
        base = {key: cell(source.get(key)) for key in GROUP_COLUMNS}
        base.update(sample_id=sample, final_sequence_id=final,
                    assignment_status=cell(source.get("candidate_assignment_status")))
        populated = False
        for rank in GROUP_RANKS:
            raw = str(source.get(rank + "_prediction", "")).strip()
            if not raw or raw.lower() in {"na", "nan", "none", "singleton", "unclassified", "-"}:
                continue
            populated = True
            kind = ("novel_group" if raw.startswith("novel_") else
                    "unplaced_group" if raw.startswith("unplaced_") else
                    "unresolved" if "|" in raw else "named_prediction")
            yield {**base, "rank": rank, "raw_prediction": raw, "prediction_type": kind,
                   "evidence_type": cell(source.get(rank + "_evidence")),
                   "network_support": cell(source.get(rank + "_network_support")),
                   "vog_support": cell(source.get(rank + "_vog_support"))}
        if not populated:
            yield {**base, **dict.fromkeys(("rank", "raw_prediction", "prediction_type",
                                           "evidence_type", "network_support", "vog_support"), "NA")}


def summary_row(sample, metadata, reasons, quality):
    final = metadata["final_sequence_id"]
    parent = metadata["normalized_name"]
    conflicts = []
    for kind, column in (("cellular", "cellular_conflict_tools"), ("plasmid", "plasmid_conflict_tools")):
        conflicts.extend(kind + ":" + tool for tool in str(metadata.get(column, "")).split(",")
                         if tool and tool != "NA")
    direct = [r for r in quality if r["assessment_scope"] == "exact_candidate"]
    qualities = {r["checkv_quality"] for r in direct}
    complete = {r["checkv_complete_call"] for r in direct}
    result = {key: cell(metadata.get(key)) for key in SUMMARY_COLUMNS}
    result.update(sample_id=sample, original_sequence_id=cell(metadata["original_contig_name"]),
                  normalized_parent_id=parent if final != parent else "NA",
                  coordinates=cell(metadata["provirus_coordinates"]), length_nt=cell(metadata["refined_length"]),
                  candidate_status=metadata["viral_decision"], conflict_tools=";".join(conflicts) or "NA",
                  reported_taxonomy=metadata["strict_taxonomy"],
                  reported_taxonomy_rank=metadata["strict_taxonomy_rank"],
                  taxonomy_support_status=metadata["strict_taxonomy_confidence"],
                  taxonomy_supporting_tools=metadata["strict_taxonomy_supporting_tools"],
                  checkv_quality=next(iter(qualities)) if len(qualities) == 1 else "NA",
                  checkv_complete_call=next(iter(complete)) if len(complete) == 1 else "NA",
                  review_flags=",".join(reasons) or "NA")
    return {key: cell(result[key]) for key in SUMMARY_COLUMNS}


def taxonomy_row(metadata):
    parts = metadata["strict_taxonomy"].split(";")
    if len(parts) != len(RANKS):
        raise ValueError("Reporting expected a complete nine-rank lineage")
    names = [part.split("__", 1)[-1] for part in parts]
    return {"final_sequence_id": metadata["final_sequence_id"],
            **{rank: "NA" if name.lower() in {"unclassified", "na", ""} else name
               for rank, name in zip(RANKS, names, strict=True)}}


def finalize(prefix, sample, metadata_rows, reviews, evidence_by_final, groups_by_final, outputs):
    """Write presentation tables, then relocate unmodified scientific artifacts."""
    prefix = Path(prefix)
    root, name = prefix.parent, prefix.name
    paths = dict(outputs)
    paths.update(taxonomy=root / (name + ".taxonomy.tsv"),
                 checkv_report=root / "quality" / (name + ".checkv.tsv"),
                 vcontact3_report=root / "groups" / (name + ".vcontact3.tsv"),
                 review_report=root / "review" / (name + ".review.tsv"))
    review_iter = iter(reviews)
    review = next(review_iter, None)
    reason_writers = {}
    with ExitStack() as stack:
        summaries = {destination: writer(stack, paths[destination + "_summary"], SUMMARY_COLUMNS)
                     for destination in ("viral_candidates", "TE_candidates")}
        tax = writer(stack, paths["taxonomy"], TAXONOMY_COLUMNS)
        checkv = writer(stack, paths["checkv_report"], QUALITY_COLUMNS)
        groups = writer(stack, paths["vcontact3_report"], GROUP_COLUMNS)
        review_table = writer(stack, paths["review_report"], REVIEW_COLUMNS)
        for metadata in metadata_rows:
            final = metadata["final_sequence_id"]
            reasons = []
            if review is not None and review["final_sequence_id"] == final:
                reasons = str(review["review_reasons"]).split(",")
                review = next(review_iter, None)
            quality = list(quality_rows(sample, metadata, evidence_by_final.get(final, ())))
            checkv.writerows(quality)
            groups.writerows(group_rows(sample, final, groups_by_final.get(final, ())))
            summary = summary_row(sample, metadata, reasons, quality)
            destination = metadata["candidate_destination"]
            if destination in summaries:
                summaries[destination].writerow(summary)
            if destination == "viral_candidates":
                tax.writerow(taxonomy_row(metadata))
            category = ("viral_candidate_flag" if destination == "viral_candidates" else
                        "TE_candidate_flag" if destination == "TE_candidates" else
                        "excluded_candidate_audit" if destination == "excluded" else "review_only")
            for reason in reasons:
                if not re.fullmatch(r"[A-Za-z0-9_]+", reason):
                    raise ValueError("Unsafe review-reason filename: " + reason)
                if reason not in reason_writers:
                    key = "review_reason_" + reason
                    paths[key] = root / "review/by_reason" / (name + ".review." + reason + ".tsv")
                    reason_writers[reason] = writer(stack, paths[key], REVIEW_COLUMNS)
                row = {**summary, "review_category": category, "review_reason": reason}
                review_table.writerow(row)
                reason_writers[reason].writerow(row)
        if review is not None:
            raise ValueError("Review records are not in final-candidate order")
    for label, source in list(paths.items()):
        if label in LEGACY_LABELS:
            destination = root / "audit/legacy" / source.name
        elif label in {"TE_candidates_fasta", "TE_candidates_summary"}:
            destination = root / "TE_candidates" / source.name
        elif label in {"metadata", "review", "disposition", "map", "candidate_routing_audit",
                       "taxonomy_reporting_audit", "taxonomy_reconciliation_audit", "evidence_audit", "taxonomy_audit"}:
            destination = root / "audit" / source.name
        else:
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source != destination:
            source.replace(destination)
        paths[label] = destination
    return paths
