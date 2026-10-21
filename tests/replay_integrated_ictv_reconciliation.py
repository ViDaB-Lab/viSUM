#!/usr/bin/env python3
"""Replay the COMPLETE harmonizer only, against saved real benchmark inputs.

No discovery, database prep, searches, Nextflow, or student data are launched.
Requires an immutable pre-integration script, published baseline, and original
staged inputs. Verifies baseline reproduction before comparing integration.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import run_viharmony as harmony
from ictv_reconciliation import Reconciler
from validate_ictv_reconciliation_replay import normalize_row

FROZEN_SHA256 = "564cc199e58836829d1bb6cf7eb440835ccbc928554b16a2fd9616d0b11cd9a9"
TAXONOMY_FIELDS = {name for name in harmony.METADATA_COLUMNS
    if name.startswith(("strict_taxonomy", "exploratory_taxonomy", "species_", "vcontact3_group"))
    or name in {"strict_species", "exploratory_species", "taxonomy_conflict"}}
IDENTICAL_OUTPUTS = ("final.normalized.fasta", "final.original_ids.fasta",
    "provisional.normalized.fasta", "provisional.original_ids.fasta",
    "database_candidates.fasta", "all_candidates.fasta", "sequence_map.tsv", "sequence_disposition.tsv")
TAXONOMY_REVIEW_REASONS = {"taxonomy_conflict", "taxonomy_unclassified",
                         "vcontact3_group_ambiguity", "vcontact3_group_context_conflict"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reporting_artifact(path):
    """Read either the frozen flat layout or the current reporting layout."""
    if path.exists():
        return path
    for folder in ("audit", "audit/legacy", "TE_candidates", "quality", "groups", "review"):
        candidate = path.parent / folder / path.name
        if candidate.exists():
            return candidate
    return path


def table(path):
    path = reporting_artifact(path)
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def by_id(path):
    rows = table(path)
    result = {row["final_sequence_id"]: row for row in rows}
    if len(result) != len(rows):
        raise AssertionError(f"Duplicate candidate IDs: {path}")
    return result


def content(path):
    path = reporting_artifact(path)
    return gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()


def verify_published_table(saved, baseline):
    """Older saved tables lack two already-existing RNA audit annotations."""
    left, right = table(saved), table(baseline)
    if len(left) != len(right):
        return False
    with saved.open(encoding="utf-8-sig") as handle:
        old_columns = next(csv.reader(handle, delimiter="\t"))
    with baseline.open(encoding="utf-8-sig") as handle:
        current_columns = next(csv.reader(handle, delimiter="\t"))
    allowed_extra = {"rna_checkv_conflict_exception", "rna_checkv_exception_supporting_loci"}
    if set(old_columns) - set(current_columns) or set(current_columns) - set(old_columns) - allowed_extra:
        return False
    return all(all(a[name] == b[name] for name in old_columns) for a, b in zip(left, right))


def verify_expected_taxonomy(baseline, integrated, sample, bundle, minimum_length):
    """Independent standalone normalization + unchanged voting for EVERY ID."""
    from collections import defaultdict
    by_final = defaultdict(list)
    for row in table(baseline / (sample + ".evidence_audit.tsv.gz")):
        if row["evidence_application"] == "used_for_final_decision":
            by_final[row["final_sequence_id"]].append(row)
    resolver = Reconciler.load(bundle)
    catalog, _ = harmony.load_msl(ROOT / "assets/ICTV_VMR_MSL41.csv")
    actual = by_id(integrated / (sample + ".final_metadata.tsv"))
    for sid, record in actual.items():
        copied = []
        for row in by_final[sid]:
            if row["tool"] not in harmony.TAXONOMY_TOOLS or row["classification"] != "virus":
                copied.append(row)
            elif row["tool"] == "vcontact3" and float(row.get("length") or 0) < minimum_length:
                copied.append(row)
            else:
                copied.append(normalize_row(row, resolver)[0])
        expected = harmony.decide_taxonomy(copied, catalog, minimum_length)
        if (record["strict_taxonomy"] != harmony.taxonomy_string(expected.strict)
                or record["strict_taxonomy_confidence"] != expected.strict_confidence
                or record["species_assignment_status"] != expected.species_status
                or record["strict_taxonomy_supporting_tools"] != (",".join(expected.strict_tools) or "NA")):
            raise AssertionError(f"Full integrated taxonomy differs from standalone expectation: {sample}/{sid}")
    return len(actual)


def launch(script, arguments, directory, sample, extra=()):
    directory.mkdir()
    command = [sys.executable, "-B", str(script), *map(str, arguments), *map(str, extra),
               "--output-prefix", str(directory / sample)]
    # This gate isolates reconciliation against the historical reporting depth.
    # The separate reporting-cap replay validates the new family default.
    if Path(script).resolve() == (ROOT / "bin/run_viharmony.py").resolve():
        command.extend(["--max-rank", "species"])
    started = time.perf_counter()
    with (directory / "replay.log").open("w", encoding="utf-8") as log:
        process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if process.returncode:
        raise AssertionError(f"Harmonizer replay failed: {directory / 'replay.log'}")
    return time.perf_counter() - started


def compare_worker_outputs(serial, parallel, sample):
    """Compare every scientific output; gzip headers may contain timestamps."""
    manifests = []
    for directory in (serial, parallel):
        manifest = json.loads(reporting_artifact(directory / (sample + ".harmonizer_manifest.json")).read_text())
        for record in manifest["outputs"].values():
            path = directory / record["file"]
            if sha(path) != record["sha256"]:
                raise AssertionError(f"Output does not match its manifest: {path}")
        manifests.append(manifest)
    left, right = manifests
    if left["outputs"].keys() != right["outputs"].keys():
        raise AssertionError("Worker output sets differ")
    for label, record in left["outputs"].items():
        other = right["outputs"][label]
        if record["file"] != other["file"] or content(serial / record["file"]) != content(parallel / other["file"]):
            raise AssertionError(f"Worker scientific outputs differ: {sample}/{label}")
        # Only gzip timestamp-dependent digests may differ, after content check.
        record.pop("sha256")
        other.pop("sha256")
    if left != right:
        raise AssertionError(f"Worker policy/provenance manifests differ: {sample}")


def compare(baseline, integrated, sample):
    before = by_id(baseline / (sample + ".final_metadata.tsv"))
    after = by_id(integrated / (sample + ".final_metadata.tsv"))
    if before.keys() != after.keys():
        raise AssertionError("Candidate population changed")
    changes = []
    for sid, row in before.items():
        new = after[sid]
        changed = {name for name in row if row[name] != new[name]}
        if changed - TAXONOMY_FIELDS:
            raise AssertionError(f"Non-taxonomy metadata changed: {sample}/{sid}/{changed - TAXONOMY_FIELDS}")
        if changed:
            changes.append({"sample": sample, "id": sid, "fields": sorted(changed),
                            "before": {k: row[k] for k in sorted(changed)},
                            "after": {k: new[k] for k in sorted(changed)}})
    for suffix in IDENTICAL_OUTPUTS:
        if reporting_artifact(baseline / (sample + "." + suffix)).read_bytes() != reporting_artifact(integrated / (sample + "." + suffix)).read_bytes():
            raise AssertionError(f"Retained/provisional sequences, mappings, or disposition changed: {sample}/{suffix}")
    if content(baseline / (sample + ".evidence_audit.tsv.gz")) != content(integrated / (sample + ".evidence_audit.tsv.gz")):
        raise AssertionError("Original evidence audit changed")
    for suffix in ("database_candidates.tsv", "provisional_metadata.tsv"):
        left, right = by_id(baseline / (sample + "." + suffix)), by_id(integrated / (sample + "." + suffix))
        if left.keys() != right.keys():
            raise AssertionError(f"Retained/provisional membership changed: {suffix}")
        for sid, row in left.items():
            if any(row[k] != right[sid][k] for k in row if k not in TAXONOMY_FIELDS):
                raise AssertionError(f"Retention metadata changed: {suffix}/{sid}")
    left_review = by_id(baseline / (sample + ".review_queue.tsv"))
    right_review = by_id(integrated / (sample + ".review_queue.tsv"))
    for sid in left_review.keys() | right_review.keys():
        reasons = []
        for review in (left_review, right_review):
            reasons.append(set(review.get(sid, {}).get("review_reasons", "").split(",")) - TAXONOMY_REVIEW_REASONS - {""})
        if reasons[0] != reasons[1]:
            raise AssertionError(f"Non-taxonomy review reasons changed: {sample}/{sid}")
    return len(before), changes


def run(args):
    parallel_workers = getattr(args, "parallel_workers", 1)
    if parallel_workers < 1:
        raise ValueError("parallel-workers must be positive")
    out = args.output.resolve()
    if out.exists():
        raise ValueError("Use a new output directory")
    for source in (args.saved_results, args.input_root, args.fasta_root, args.frozen_script.parent):
        if out == source.resolve() or source.resolve() in out.parents:
            raise ValueError("Output must be outside every source input tree")
    if sha(args.frozen_script) != FROZEN_SHA256:
        raise ValueError("Frozen script is not the verified pre-integration harmonizer")
    out.mkdir(parents=True)
    source_hashes, panels, changes = {}, [], []
    status = {"status": "RUNNING", "scope": "complete harmonizer; saved fixed downstream inputs; no searches",
              "frozen_script_sha256": sha(args.frozen_script), "integrated_script_sha256": sha(ROOT / "bin/run_viharmony.py")}
    try:
        for panel in sorted(args.saved_results.iterdir()):
            sample = panel.name
            if not (panel / "updated" / (sample + ".final_metadata.tsv")).is_file():
                continue
            source = args.input_root / (sample + "_results")
            fasta = args.fasta_root / (sample + "_results")
            saved = panel / "updated"
            manifest = json.loads((saved / (sample + ".harmonizer_manifest.json")).read_text())
            evidence = []
            for tool in ("genomad", "virsorter2", "cenotetaker3", "deep6", "deepmicroclass2",
                         "virbot", "gianthunter", "vicat", "checkv", "tesorter", "vitap", "vcontact3"):
                path = source / tool / (sample + "." + tool + "_evidence.tsv")
                if tool == "virsorter2":
                    path = panel / (sample + ".virsorter2_evidence.tsv")
                if path.is_file():
                    evidence.append(path)
            evidence.append(source / "vicat" / (sample + ".vicat_refined_evidence.tsv"))
            if manifest["input_type"] == "rna":
                evidence.append(source / "vicat" / (sample + ".vicat_orf_evidence.tsv"))
            paths = {
                "normalized-fasta": fasta / "prep" / (sample + ".normalized.fasta"),
                "refined-fasta": fasta / "refinement" / (sample + ".refined_candidates.fasta"),
                "header-map": source / "prep" / (sample + ".header_map.tsv"),
                "discovery-gate": source / "discovery_gate" / (sample + ".discovery_gate.tsv"),
                "region-map": source / "refinement" / (sample + ".provirus_region_map.tsv"),
                "ictv-msl": ROOT / "assets/ICTV_VMR_MSL41.csv",
            }
            groups = source / "vcontact3" / (sample + ".vcontact3_group_membership.tsv")
            inputs = [*paths.values(), *evidence, groups, args.bundle,
                      args.bundle.with_name(args.bundle.name + ".manifest.json"), *[p for p in saved.iterdir() if p.is_file()]]
            for path in inputs:
                source_hashes[str(path.resolve())] = sha(path)
            arguments = ["--sample-id", sample, "--input-type", manifest["input_type"],
                         "--audit-mode", "full", "--workers", "1",
                         "--vcontact3-min-taxonomy-length", str(manifest["policy"]["vcontact3_min_taxonomy_length"])]
            for name, path in paths.items():
                arguments.extend(["--" + name, str(path.resolve())])
            arguments.extend(["--evidence", *map(lambda p: str(p.resolve()), evidence), "--vcontact3-groups", str(groups.resolve())])
            if manifest["policy"]["rna_pair_homology_floor"]:
                arguments.append("--rna-pair-homology-floor")
            historical_exception = manifest["policy"].get("rna_checkv_conflict_exception", False)
            if manifest["input_type"] == "rna" and not historical_exception:
                arguments.append("--disable-rna-checkv-exception")
            dest = out / sample
            dest.mkdir()
            launch(args.frozen_script, arguments, dest / "baseline", sample)
            # Every existing scientific output must reproduce its published
            # contents (gzip timestamps alone may differ).
            for record in manifest["outputs"].values():
                name = record["file"]
                identical = content(saved / name) == content(dest / "baseline" / name)
                if not identical and name.endswith(("database_candidates.tsv", "final_metadata.tsv", "provisional_metadata.tsv", "review_queue.tsv")):
                    identical = verify_published_table(saved / name, dest / "baseline" / name)
                if not identical:
                    raise AssertionError(f"Frozen full baseline does not reproduce published output: {sample}/{name}")
            print(f"PASS frozen published baseline: {sample}", flush=True)
            serial_seconds = launch(ROOT / "bin/run_viharmony.py", arguments, dest / "integrated", sample,
                   ("--ictv-reconciliation", args.bundle.resolve()))
            parallel_seconds = None
            if parallel_workers > 1:
                parallel_arguments = list(arguments)
                parallel_arguments[parallel_arguments.index("--workers") + 1] = str(parallel_workers)
                parallel_seconds = launch(ROOT / "bin/run_viharmony.py", parallel_arguments,
                    dest / "integrated_parallel", sample, ("--ictv-reconciliation", args.bundle.resolve()))
                compare_worker_outputs(dest / "integrated", dest / "integrated_parallel", sample)
            count, panel_changes = compare(dest / "baseline", dest / "integrated", sample)
            expected_count = verify_expected_taxonomy(dest / "baseline", dest / "integrated", sample,
                args.bundle, manifest["policy"]["vcontact3_min_taxonomy_length"])
            current_default_differences = []
            if manifest["input_type"] == "rna" and not historical_exception:
                current_arguments = [arg for arg in arguments if arg != "--disable-rna-checkv-exception"]
                launch(args.frozen_script, current_arguments, dest / "baseline_current_defaults", sample)
                launch(ROOT / "bin/run_viharmony.py", current_arguments, dest / "integrated_current_defaults", sample,
                       ("--ictv-reconciliation", args.bundle.resolve()))
                compare(dest / "baseline_current_defaults", dest / "integrated_current_defaults", sample)
                if parallel_workers > 1:
                    current_parallel = list(current_arguments)
                    current_parallel[current_parallel.index("--workers") + 1] = str(parallel_workers)
                    launch(ROOT / "bin/run_viharmony.py", current_parallel,
                        dest / "integrated_current_defaults_parallel", sample,
                        ("--ictv-reconciliation", args.bundle.resolve()))
                    compare_worker_outputs(dest / "integrated_current_defaults",
                        dest / "integrated_current_defaults_parallel", sample)
                historical = by_id(dest / "baseline" / (sample + ".final_metadata.tsv"))
                current = by_id(dest / "baseline_current_defaults" / (sample + ".final_metadata.tsv"))
                integrated_historical = by_id(dest / "integrated" / (sample + ".final_metadata.tsv"))
                integrated_current = by_id(dest / "integrated_current_defaults" / (sample + ".final_metadata.tsv"))
                for sid, row in current.items():
                    if row["viral_decision"] != historical[sid]["viral_decision"]:
                        current_default_differences.append({"id": sid, "historical": historical[sid]["viral_decision"],
                                                           "preexisting_current_default": row["viral_decision"]})
                    # Existing origin-only CheckV exception must not change
                    # expected strict taxonomy when current defaults are used.
                    for key in ("strict_taxonomy", "strict_taxonomy_confidence", "species_assignment_status", "strict_taxonomy_supporting_tools"):
                        if integrated_current[sid][key] != integrated_historical[sid][key]:
                            raise AssertionError(f"Strict taxonomy unexpectedly depends on CheckV exception: {sample}/{sid}/{key}")
            changes.extend(panel_changes)
            panels.append({"sample": sample, "candidate_records": count,
                "unchanged_retention_and_origin_metadata": True, "unchanged_sequence_outputs": list(IDENTICAL_OUTPUTS),
                "unchanged_original_evidence_audit": True, "changed_taxonomy_records": len(panel_changes),
                "standalone_expected_taxonomy_verified": expected_count,
                "current_default_retention_verified": True,
                "historical_checkv_exception_enabled": historical_exception,
                "parallel_workers": parallel_workers,
                "parallel_outputs_identical": parallel_workers > 1,
                "integrated_serial_seconds": serial_seconds,
                "integrated_parallel_seconds": parallel_seconds,
                "preexisting_current_default_vs_historical_retention_changes": current_default_differences})
            print(f"PASS complete integrated replay: {sample}; {count} candidates; retention unchanged", flush=True)
        if not panels:
            raise ValueError("No complete saved panels found")
        expected_path = args.expected_taxonomy / "taxonomy_changes.jsonl"
        with expected_path.open(encoding="utf-8") as handle:
            expected = [json.loads(line) for line in handle if line.strip()]
        for item in expected:
            row = by_id(out / item["sample"] / "integrated" / (item["sample"] + ".final_metadata.tsv"))[item["id"]]
            decision = item["after"]
            if (row["strict_taxonomy"] != harmony.taxonomy_string(decision["taxonomy"])
                    or row["strict_taxonomy_confidence"] != decision["confidence"]
                    or row["species_assignment_status"] != decision["species_status"]
                    or row["strict_taxonomy_supporting_tools"] != (",".join(decision["tools"]) or "NA")):
                raise AssertionError(f"Integrated taxonomy differs from validated prototype: {item['sample']}/{item['id']}")
        if any(sha(Path(path)) != digest for path, digest in source_hashes.items()):
            raise AssertionError("A saved source input changed")
        if sha(ROOT / "bin/run_viharmony.py") != status["integrated_script_sha256"] or sha(args.frozen_script) != FROZEN_SHA256:
            raise AssertionError("A replayed scientific source changed during validation")
        status.update(status="PASS", retention_recomputed=True,
            retention_decisions_identical=True, retained_and_provisional_sequences_identical=True,
            saved_inputs_unchanged=True, candidate_records=sum(p["candidate_records"] for p in panels),
            published_table_schema_note="Older published tables omit two pre-existing RNA CheckV annotation columns; all their existing fields reproduced exactly, all published sequence outputs reproduced byte-for-byte.",
            standalone_expected_taxonomy_verified=sum(p["standalone_expected_taxonomy_verified"] for p in panels),
            expected_changed_taxonomy_records_verified=len(expected), samples=panels, input_sha256=source_hashes)
        status["parallel_workers"] = parallel_workers
        status["parallel_outputs_identical"] = parallel_workers > 1
        if {(r["sample"], r["id"]) for r in changes} != {(r["sample"], r["id"]) for r in expected}:
            raise AssertionError("Changed taxonomy IDs differ from the validated standalone prototype")
    except Exception:
        import traceback
        status.update(status="FAIL", error=traceback.format_exc(), samples=panels)
        raise
    finally:
        (out / "report.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        with (out / "taxonomy_changes.jsonl").open("w", encoding="utf-8") as handle:
            for row in changes:
                handle.write(json.dumps(row, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("saved-results", "input-root", "fasta-root", "frozen-script", "bundle", "expected-taxonomy", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--parallel-workers", type=int, default=1)
    run(parser.parse_args())
