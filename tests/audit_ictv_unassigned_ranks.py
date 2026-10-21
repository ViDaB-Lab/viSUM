#!/usr/bin/env python3
"""Audit missing formal ancestors in the pinned target, without pipeline work.

Injected higher labels below are deliberately synthetic test challenges, not
observed tool calls or claims that ICTV should assign these missing ranks.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from ictv_reconciliation import RANKS, Reconciler, sha256


def audit_release(resolver):
    rank_nodes = {rank: sorted((nid for nid, node in resolver.nodes.items() if node["rank"] == rank),
                              key=lambda nid: resolver.nodes[nid]["name"]) for rank in RANKS}
    candidate_cache = {}
    recovery_cases, conflict_cases, recoverable_taxa = 0, 0, 0
    rows, counts, gap_counts, top_level = [], Counter(), Counter(), []
    for nid, node in sorted(resolver.nodes.items()):
        if node["rank"] not in RANKS:
            continue
        expected = resolver.lineages[nid]
        gaps = [rank for rank in RANKS[:RANKS.index(node["rank"])] if rank not in expected]
        if not gaps:
            continue
        baseline = resolver.resolve(expected)
        if baseline.normalized != expected or baseline.status != "unchanged":
            raise AssertionError(f"Valid unassigned ancestry rejected: {node['name']}")
        recovered_gaps, conflicting_gaps = [], []
        for rank in gaps:
            # Challenge one gap at a time. Prefer a distractor compatible with
            # explicitly assigned upper ancestry, so synthetic labels do not
            # accidentally create a genuine conflict between upper claims.
            upper = tuple((r, expected[r]) for r in RANKS[:RANKS.index(rank)] if r in expected)
            cache_key = rank, upper
            if cache_key not in candidate_cache:
                candidate_cache[cache_key] = next((n for n in rank_nodes[rank]
                    if all(r not in resolver.lineages[n] or resolver.lineages[n][r] == name
                           for r, name in upper)), None)
            candidate = candidate_cache[cache_key]
            challenge = {**expected, rank: resolver.nodes[candidate or rank_nodes[rank][0]]["name"]}
            corrected = resolver.resolve(challenge)
            if corrected.original != challenge:
                raise AssertionError(f"Original challenge lost: {node['name']}")
            if candidate is None:
                if corrected.status != "incompatible" or corrected.normalized:
                    raise AssertionError(f"Assigned-ancestor conflict hidden: {node['name']}/{rank}")
                conflict_cases += 1
                conflicting_gaps.append(rank)
            else:
                reasons = {item["rank"]: item["reason"] for item in corrected.unresolved}
                if (corrected.normalized != expected or corrected.status != "partial"
                        or reasons != {rank: "target_higher_rank_unassigned"}):
                    raise AssertionError(f"Unsupported higher label not safely flagged: {node['name']}/{rank}")
                recovery_cases += 1
                recovered_gaps.append(rank)
        recoverable_taxa += bool(recovered_gaps)
        counts[node["rank"]] += 1
        gap_counts.update(gaps)
        row = {"target_node_id": nid, "rank": node["rank"], "name": node["name"],
               "missing_higher_ranks": gaps, "lineage": expected,
               "synthetic_recovery_tested_ranks": recovered_gaps,
               "synthetic_assigned_conflict_tested_ranks": conflicting_gaps}
        rows.append(row)
        if not any(rank in expected for rank in RANKS[:RANKS.index(node["rank"])]):
            top_level.append(row)
    return {"target_release": resolver.target_release, "prototype_only": True,
        "challenge_type": "one synthetic known label per missing higher rank; preserve compatible lower ancestry, but reject any explicit assigned-ancestor conflict",
        "taxa_with_missing_higher_formal_ranks": len(rows),
        "by_taxon_rank": dict(counts), "by_missing_rank": dict(gap_counts),
        "unchanged_valid_lineages_verified": len(rows),
        "synthetic_recovery_cases_verified": recovery_cases,
        "taxa_with_recoverable_challenges": recoverable_taxa,
        "synthetic_assigned_conflicts_verified": conflict_cases,
        "top_level_taxa_below_realm": top_level, "taxa": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Use a new audit output; never overwrite an earlier report")
    report = audit_release(Reconciler.load(args.bundle))
    report["bundle_sha256"] = sha256(args.bundle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k not in {"taxa", "top_level_taxa_below_realm"}}, indent=2))


if __name__ == "__main__":
    main()
