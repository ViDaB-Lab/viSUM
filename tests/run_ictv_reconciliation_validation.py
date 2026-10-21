#!/usr/bin/env python3
"""Dependency-free validation of the standalone, self-contained Linux package.

Never launches Nextflow, searches, setup, or deployment. Reads only this
package's snapshots and writes a new validation directory.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
sys.path.insert(0, str(ROOT / "tests"))
from ictv_reconciliation import sha256
from audit_ictv_unassigned_ranks import audit_release
from ictv_reconciliation import Reconciler
from validate_ictv_reconciliation_replay import validate

TEST_MODULES = ("test_ictv_reconciliation.py", "test_ictv_reference.py",
                "test_viharmony.py", "test_viharmony_workflow_contracts.py")


def verify_package(manifest):
    for name, expected in manifest["files"].items():
        path = (ROOT / name).resolve()
        if ROOT not in path.parents or not path.is_file() or sha256(path) != expected:
            raise AssertionError(f"Package file missing, outside package, or changed: {name}")


def read_jsonl(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]


def canonical_rows(path):
    return Counter(json.dumps(row, sort_keys=True) for row in read_jsonl(path))


def result_hashes(report):
    # Preserve the sample/updated/filename identity, not Windows/Linux roots.
    result = {}
    for name, digest in report["input_sha256"].items():
        key = "/".join(name.replace("\\", "/").split("/")[-3:])
        if key in result:
            raise AssertionError(f"Non-unique saved-file identity: {key}")
        result[key] = digest
    return result


def review_examples(output, previous):
    changes = read_jsonl(output / "benchmark_replay/taxonomy_changes.jsonl")
    resolutions = read_jsonl(output / "benchmark_replay/reconciliation_audit.jsonl")
    old = {(r["sample"], r["id"]): r for r in read_jsonl(previous / "taxonomy_changes.jsonl")}
    flagged = {(r["sample"], r["id"]): r for r in resolutions
               if any(u["reason"] == "target_higher_rank_unassigned" for u in r["unresolved"])}
    selected, used = [], set()
    groups = (
        ("Unassigned higher rank: final result differs from earlier prototype",
         lambda r: (r["sample"], r["id"]) in flagged and
         r["after"] != old.get((r["sample"], r["id"]), {}).get("after")),
        ("Taxonomic names or depth changed versus unreconciled baseline",
         lambda r: r["before"]["taxonomy"] != r["after"]["taxonomy"]),
        ("Same taxonomy, different supporting tools/confidence",
         lambda r: r["before"]["taxonomy"] == r["after"]["taxonomy"]),
    )
    for label, predicate in groups:
        count = 0
        for row in changes:
            key = row["sample"], row["id"]
            if key not in used and predicate(row):
                selected.append({"review_category": label, **row,
                                 "flagged_evidence_example": flagged.get(key)})
                used.add(key)
                count += 1
                if count == 2:
                    break
    (output / "representative_taxonomy_changes.json").write_text(
        json.dumps(selected, indent=2) + "\n", encoding="utf-8")
    text = ["STANDALONE TAXONOMY REVIEW — NOT DEPLOYED",
            "Before = unchanged engine on original evidence; after = the same engine on reconciled taxonomy fields.",
            "Saved viral decisions are context only: retention was NOT recomputed.",
            "Reference compatibility/support is not proof of biological accuracy.",
            "Examples include nominal negative controls; recovered taxonomy does not validate their source labels or viral identity.", ""]
    for row in selected:
        text.extend([row["review_category"], f"{row['sample']} / {row['id']} / original: {row['original_id']}",
                     f"Saved viral decision (unchanged): {row['saved_viral_decision']}"])
        for stage in ("before", "after"):
            value = row[stage]
            lineage = "; ".join(f"{k}={v}" for k, v in value["taxonomy"].items() if v) or "unassigned"
            text.extend([f"{stage.upper()}: {lineage}",
                         f"  tools={value['tools']}; confidence={value['confidence']}; species status={value['species_status']}"])
        evidence = row["flagged_evidence_example"]
        if evidence:
            text.extend([f"Flagged source: {evidence['tool']}",
                         "Original: " + json.dumps(evidence["original"], sort_keys=True),
                         "Normalized: " + json.dumps(evidence["normalized"], sort_keys=True),
                         "Flags: " + json.dumps(evidence["unresolved"], sort_keys=True)])
        text.append("")
    (output / "representative_taxonomy_review.txt").write_text("\n".join(text) + "\n", encoding="utf-8")
    return len(selected)


def run(output):
    output = output.resolve()
    if output.exists():
        raise ValueError("Use a new output directory; never overwrite previous validation")
    if output == ROOT or (ROOT / "saved_results") == output or (ROOT / "saved_results") in output.parents:
        raise ValueError("Output cannot replace package source or saved results")
    output.mkdir(parents=True)
    status = {"status": "RUNNING", "prototype_only": True,
              "python": sys.version, "platform": platform.platform(), "retention_recomputed": False}
    try:
        manifest = json.loads((ROOT / "PACKAGE_MANIFEST.json").read_text(encoding="utf-8"))
        verify_package(manifest)
        status["package_manifest_sha256"] = sha256(ROOT / "PACKAGE_MANIFEST.json")
        print("PASS: package source, bundle, and saved-result checksums", flush=True)
        suite = unittest.TestSuite()
        for pattern in TEST_MODULES:
            suite.addTests(unittest.defaultTestLoader.discover(str(ROOT / "tests"), pattern=pattern))
        with (output / "regression.log").open("w", encoding="utf-8") as log:
            with redirect_stdout(log), redirect_stderr(log):
                result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
        status["regression"] = {"tests_run": result.testsRun, "failures": len(result.failures),
                                "errors": len(result.errors), "skipped": len(result.skipped)}
        if not result.wasSuccessful() or result.skipped or result.testsRun != manifest["expected_unittest_count"]:
            raise AssertionError("Regression tests failed, skipped, or did not run the expected count; see regression.log")
        print(f"PASS: {result.testsRun} dependency-free regression tests (no skips)", flush=True)
        bundle = ROOT / manifest["bundle"]
        release_report = audit_release(Reconciler.load(bundle))
        release_report["bundle_sha256"] = sha256(bundle)
        (output / "release_audit.json").write_text(json.dumps(release_report, indent=2) + "\n", encoding="utf-8")
        expected = ROOT / "expected_corrected"
        if release_report != json.loads((expected / "release_audit.json").read_text(encoding="utf-8")):
            raise AssertionError("Release-wide audit differs from validated local snapshot")
        print("PASS: complete unassigned-rank audit matches validated local snapshot", flush=True)
        report = validate(ROOT / "saved_results", bundle, output / "benchmark_replay",
                          samples=manifest["samples"], affected_only=True)
        reference = json.loads((expected / "report.json").read_text(encoding="utf-8"))
        for name in ("samples", "target_release", "bundle_sha256", "baseline_mismatch_count",
                     "protected_source_sha256", "input_and_fasta_checksums_unchanged",
                     "protected_source_checksums_unchanged", "retention_recomputed", "retention_rules_modified"):
            if report[name] != reference[name]:
                raise AssertionError(f"Saved-evidence comparison differs from local validation: {name}")
        if result_hashes(report) != result_hashes(reference):
            raise AssertionError("Saved result contents differ from local validation")
        for name in ("taxonomy_changes.jsonl", "reconciliation_audit.jsonl"):
            if canonical_rows(output / "benchmark_replay" / name) != canonical_rows(expected / name):
                raise AssertionError(f"Detailed taxonomy outputs differ from local validation: {name}")
        examples = review_examples(output, ROOT / "expected_previous")
        verify_package(manifest)
        status.update(status="PASS", exact_local_snapshot_match=True,
            input_and_package_checksums_unchanged=True,
            saved_baselines_reproduced=sum(s.get("baseline_matched", 0) for s in report["samples"]),
            recovered_sequence_tool_pairs=sum(s["recovered_sequence_tool_pairs"] for s in report["samples"]),
            representative_examples=examples)
        print("PASS: benchmark summaries and detailed taxonomy outputs exactly match local validation", flush=True)
        print(f"Review: {output / 'representative_taxonomy_review.txt'}", flush=True)
    except Exception:
        status["status"] = "FAIL"
        status["error"] = traceback.format_exc()
        (output / "validation_status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise
    (output / "validation_status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    print(f"Final status: {status['status']}; report: {output / 'validation_status.json'}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    run(args.output_dir)


if __name__ == "__main__":
    main()
