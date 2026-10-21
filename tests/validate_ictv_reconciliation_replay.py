#!/usr/bin/env python3
"""Compare taxonomy only on saved audits; never launch tools or edit inputs.

Retention is NOT recomputed here. Scientific source, source FASTAs, metadata,
and evidence hashes are checked before/after. The existing taxonomy engine is
called unchanged; this harness replaces only rank/name fields on copied rows.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
import csv
import gzip
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import run_viharmony as harmony
from ictv_reconciliation import RANKS, Reconciler, sha256


PROTECTED = ("bin/run_viharmony.py", "visum_nextflow.nf", "visum.config", "nextflow.config",
             "modules/local/viharmony.nf", "bin/standardize_vicat.py", "bin/standardize_vitap.py",
             "bin/standardize_vcontact3.py")


def read_tsv(path):
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_row(row, resolver):
    """Do not modify classification, evidence strength, confidence, or scores."""
    original = {rank: harmony.clean_taxon(row.get(harmony.RANK_COLUMNS[rank], ""), rank) for rank in RANKS}
    decision = resolver.resolve(original)
    copied = dict(row)
    if decision.normalized:
        for rank in RANKS:
            column = harmony.RANK_COLUMNS[rank]
            copied[column] = harmony.PREFIXES[rank] + decision.normalized.get(rank, "unclassified")
    non_taxonomy = set(row) - {harmony.RANK_COLUMNS[r] for r in RANKS}
    assert all(row[key] == copied[key] for key in non_taxonomy)
    assert original == {rank: harmony.clean_taxon(row.get(harmony.RANK_COLUMNS[rank], ""), rank) for rank in RANKS}
    return copied, decision


def validate(results_root, bundle, output_dir, samples=None, affected_only=False):
    results_root, output_dir = results_root.resolve(), output_dir.resolve()
    if output_dir == results_root or results_root in output_dir.parents:
        raise ValueError("Output must be outside the saved results tree")
    if output_dir.exists():
        raise ValueError("Use a new output directory; never overwrite a previous validation")
    resolver = Reconciler.load(bundle)
    catalog, release = harmony.load_msl(ROOT / "assets/ICTV_VMR_MSL41.csv")
    actual = {tuple(resolver.lineages[nid].get(rank, "") for rank in RANKS)
              for nid, node in resolver.nodes.items() if node["rank"] == "species"}
    if actual != set(catalog.lineages):
        raise ValueError("Reconciliation target and unchanged taxonomy engine reference differ")
    protected_before = {name: sha256(ROOT / name) for name in PROTECTED}
    input_before = {}
    summaries, changes, resolutions = [], [], []
    for folder in sorted(results_root.iterdir()):
        sample = folder.name
        if not folder.is_dir() or (samples and sample not in samples):
            continue
        directory = folder / "updated"
        audit = directory / (sample + ".evidence_audit.tsv.gz")
        metadata = directory / (sample + ".final_metadata.tsv")
        if not audit.exists() or not metadata.exists():
            continue
        manifest_path = directory / (sample + ".harmonizer_manifest.json")
        manifest = json.loads(manifest_path.read_text())
        minimum_length = manifest["policy"]["vcontact3_min_taxonomy_length"]
        for path in directory.iterdir():
            if path.is_file():
                input_before[str(path)] = sha256(path)
        by_final = defaultdict(list)
        for row in read_tsv(audit):
            if row.get("evidence_application") == "used_for_final_decision":
                by_final[row["final_sequence_id"]].append(row)
        counters, statuses, affected, recovered, tool_pairs = Counter(), Counter(), set(), set(), set()
        unassigned_by_tool, unassigned_by_rank = Counter(), Counter()
        records = read_tsv(metadata)
        affected_ids = None
        if affected_only:
            affected_ids = {r["final_sequence_id"] for r in read_tsv(directory / (sample + ".taxonomy_vote_audit.tsv.gz"))
                            if r["reason"] == "invalid_configured_ictv_msl"}
        for record in records:
            sid = record["final_sequence_id"]
            if affected_ids is not None and sid not in affected_ids:
                counters["unaffected_records_not_replayed"] += 1
                continue
            rows = by_final.get(sid, [])
            before = harmony.decide_taxonomy(rows, catalog, minimum_length)
            if (harmony.taxonomy_string(before.strict), before.species_status, before.strict_confidence) != (
                    record["strict_taxonomy"], record["species_assignment_status"], record["strict_taxonomy_confidence"]):
                raise AssertionError(f"Saved baseline could not be reproduced: {sample}/{sid}")
            counters["baseline_matched"] += 1
            copied_rows = []
            for row in rows:
                tool = row.get("tool", "")
                if tool not in harmony.TAXONOMY_TOOLS or row.get("classification") != "virus":
                    copied_rows.append(dict(row)); continue
                original = {rank: harmony.clean_taxon(row.get(harmony.RANK_COLUMNS[rank], ""), rank) for rank in RANKS}
                if not any(original.values()):
                    copied_rows.append(dict(row)); continue
                if tool == "vcontact3" and float(row.get("length") or 0) < minimum_length:
                    copied_rows.append(dict(row)); continue
                valid_before = bool(catalog.matching(original))
                new_row, resolution = normalize_row(row, resolver)
                copied_rows.append(new_row)
                statuses[resolution.status] += 1
                unassigned = [item for item in resolution.unresolved
                              if item["reason"] == "target_higher_rank_unassigned"]
                if unassigned:
                    counters["unassigned_higher_rank_rows"] += 1
                    unassigned_by_tool[tool] += 1
                    unassigned_by_rank.update(item["rank"] for item in unassigned)
                if not valid_before:
                    affected.add(sid)
                    tool_pairs.add((sid, tool))
                    normalized = {rank: harmony.clean_taxon(new_row.get(harmony.RANK_COLUMNS[rank], ""), rank) for rank in RANKS}
                    if resolution.normalized and catalog.matching(normalized):
                        recovered.add((sid, tool))
                    resolutions.append({"sample": sample, "id": sid, "tool": tool, **asdict(resolution)})
            after = harmony.decide_taxonomy(copied_rows, catalog, minimum_length)
            names_changed = harmony.taxonomy_string(before.strict) != harmony.taxonomy_string(after.strict)
            support_changed = (before.strict_tools, before.strict_confidence, before.species_status) != (after.strict_tools, after.strict_confidence, after.species_status)
            counters["changed_names_or_depth"] += names_changed
            counters["changed_support_or_species_status"] += support_changed
            if names_changed or support_changed:
                changes.append({"sample": sample, "id": sid, "original_id": record["original_contig_name"],
                    "saved_viral_decision": record["viral_decision"],
                    "before": {"taxonomy": before.strict, "tools": before.strict_tools, "confidence": before.strict_confidence, "species_status": before.species_status},
                    "after": {"taxonomy": after.strict, "tools": after.strict_tools, "confidence": after.strict_confidence, "species_status": after.species_status}})
        summaries.append({"sample": sample, "candidate_records": len(records), **dict(counters),
            "affected_sequences": len(affected), "rejected_sequence_tool_pairs": len(tool_pairs),
            "recovered_sequence_tool_pairs": len(recovered),
            "recovered_sequences": len({sid for sid, tool in recovered}),
            "recovered_by_tool": dict(Counter(tool for sid, tool in recovered)),
            "unassigned_higher_rank_by_tool": dict(unassigned_by_tool),
            "unassigned_higher_rank_by_rank": dict(unassigned_by_rank),
            "resolution_status_counts": dict(statuses)})
        print(f"Validated {sample}: {counters['baseline_matched']} saved taxonomy baselines; {len(recovered)} recovered sequence/tool pairs", flush=True)
    if not summaries:
        raise ValueError("No complete saved audit panels found")
    if samples and set(samples) != {s["sample"] for s in summaries}:
        raise ValueError("A requested saved panel is missing")
    if any(sha256(Path(path)) != digest for path, digest in input_before.items()):
        raise AssertionError("A saved input changed during validation")
    if any(sha256(ROOT / name) != digest for name, digest in protected_before.items()):
        raise AssertionError("A protected scientific/workflow source changed during validation")
    output_dir.mkdir(parents=True)
    report = {"prototype_only": True, "target_release": resolver.target_release,
        "source_commit": resolver.bundle["source_commit"], "bundle_sha256": sha256(bundle),
        "source_release_handling": "Unknown per-tool source release; exact current names use current meaning; obsolete names use historical relationships.",
        "unassigned_higher_rank_handling": "Preserve jointly compatible lower placements only when all surviving lower targets lack that higher rank; flag and audit omitted higher predictions; assigned-ancestor conflicts remain incompatible.",
        "affected_only": affected_only,
        "retention_recomputed": False, "retention_rules_modified": False,
        "input_and_fasta_checksums_unchanged": True, "protected_source_checksums_unchanged": True,
        "protected_source_sha256": protected_before, "input_sha256": input_before,
        "prototype_source_sha256": {name: sha256(ROOT / name) for name in (
            "bin/ictv_reconciliation.py", "bin/build_ictv_reconciliation_bundle.py",
            "tests/validate_ictv_reconciliation_replay.py", "tests/audit_ictv_unassigned_ranks.py")},
        "baseline_mismatch_count": 0, "samples": summaries}
    (output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for name, rows in (("taxonomy_changes.jsonl", changes), ("reconciliation_audit.jsonl", resolutions)):
        with (output_dir / name).open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--affected-only", action="store_true", help="Replay records with saved ICTV reference rejection; explicitly report skipped unaffected records")
    args = parser.parse_args()
    report = validate(args.results_root, args.bundle, args.output_dir, args.samples, args.affected_only)
    print(json.dumps({"report": str(args.output_dir / "report.json"), "samples": report["samples"]}, indent=2))


if __name__ == "__main__":
    main()
