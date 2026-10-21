#!/usr/bin/env python3
"""Compare reporting ceilings with an independent pre-cap saved-input replay.

Executes Python viHARMONY only. Never calls Nextflow, searches, or student paths.
The saved root is the earlier integrated-reconciliation package's saved/ folder.
"""
from __future__ import annotations
import argparse
import ast
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import run_viharmony as harmony


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def content(path):
    data = path.read_bytes()
    return gzip.decompress(data) if path.suffix == ".gz" else data


def records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def launch(script, arguments, directory, sample, rank=None, workers=1):
    directory.mkdir(parents=True)
    command = [sys.executable, "-B", str(script), *arguments,
               "--workers", str(workers), "--output-prefix", str(directory / sample)]
    if rank is not None:
        command.extend(["--max-rank", rank])
    result = subprocess.run(command, text=True, capture_output=True, timeout=600)
    (directory / "run.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Harmonizer failed: {directory}\n{result.stdout}{result.stderr}")


def verify_manifest(directory, sample):
    manifest = json.loads((directory / (sample + ".harmonizer_manifest.json")).read_text())
    for item in manifest["outputs"].values():
        if sha(directory / item["file"]) != item["sha256"]:
            raise AssertionError("Output checksum mismatch: " + item["file"])
    return manifest


def policy_functions_unchanged(before, after):
    old, new = ast.parse(before.read_text(encoding="utf-8")), ast.parse(after.read_text(encoding="utf-8"))
    functions = lambda tree: {n.name: ast.dump(n, include_attributes=False) for n in tree.body
                              if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
    left, right = functions(old), functions(new)
    permitted = {"parse_args", "evaluate_candidate", "run"}
    checked = [name for name in left if name not in permitted]
    for name in checked:
        if left[name] != right[name]:
            raise AssertionError("Pre-existing policy/helper function changed: " + name)
    return len(checked)


def verify_cap(baseline, capped, sample, rank):
    old = records(baseline / (sample + ".final_metadata.tsv"))
    new = records(capped / (sample + ".final_metadata.tsv"))
    audits = records(capped / (sample + ".taxonomy_reporting_audit.tsv.gz"))
    if not len(old) == len(new) == len(audits):
        raise AssertionError("Candidate population changed")
    changed = 0
    for original, reported, audit in zip(old, new, audits, strict=True):
        for name in original:
            if name not in harmony.REPORTING_TAXONOMY_COLUMNS and original[name] != reported[name]:
                raise AssertionError(f"Non-taxonomy field changed: {sample}/{name}")
        for name in harmony.REPORTING_TAXONOMY_COLUMNS:
            if audit["uncapped_" + name] != original[name]:
                raise AssertionError(f"Uncapped consensus lost: {sample}/{name}")
        for mode in ("strict", "exploratory"):
            values = original[mode + "_taxonomy"].split(";")
            expected = ";".join(value if i <= harmony.RANKS.index(rank)
                                else harmony.PREFIXES[harmony.RANKS[i]] + "unclassified"
                                for i, value in enumerate(values))
            if reported[mode + "_taxonomy"] != expected:
                raise AssertionError("Reported lineage differs from independent rank mask")
        if audit["harmonizer_max_rank"] != rank:
            raise AssertionError("Incorrect audit rank ceiling")
        changed += audit["reporting_status"] == "capped"
    # Includes primary/provisional/all/review FASTAs: presentation-only changes
    # must not change any sequence or identifier in the existing output contract.
    for path in baseline.glob("*.fasta"):
        if path.read_bytes() != (capped / path.name).read_bytes():
            raise AssertionError("FASTA changed: " + path.name)
    for suffix in ("sequence_map.tsv", "sequence_disposition.tsv", "evidence_audit.tsv.gz",
                   "taxonomy_vote_audit.tsv.gz", "taxonomy_reconciliation_audit.tsv.gz"):
        if content(baseline / (sample + "." + suffix)) != content(capped / (sample + "." + suffix)):
            raise AssertionError("Original evidence/decisions/uncapped votes changed: " + suffix)
    for suffix in ("review_queue.tsv", "database_candidates.tsv", "provisional_metadata.tsv"):
        left, right = records(baseline / (sample + "." + suffix)), records(capped / (sample + "." + suffix))
        if len(left) != len(right):
            raise AssertionError("Output membership changed: " + suffix)
        for a, b in zip(left, right, strict=True):
            for name in a:
                if name not in harmony.REPORTING_TAXONOMY_COLUMNS and a[name] != b[name]:
                    raise AssertionError(f"Decision/review reason changed: {sample}/{suffix}/{name}")
    return len(old), changed


def run(args):
    saved = args.saved_root.resolve()
    out = args.output.resolve()
    if out.exists():
        raise ValueError("Choose a new output directory; existing results are never overwritten")
    if saved == out or saved in out.parents:
        raise ValueError("Output must be separate from saved inputs")
    out.mkdir(parents=True)
    baseline_script = args.baseline_script.resolve()
    script = ROOT / "bin/run_viharmony.py"
    protected = {str(p.resolve()): sha(p) for p in saved.rglob("*") if p.is_file()}
    report = {"status": "RUNNING", "scope": "saved benchmark inputs; Python harmonizer only; no searches",
              "baseline_script_sha256": sha(baseline_script), "updated_script_sha256": sha(script),
              "unchanged_preexisting_policy_functions": policy_functions_unchanged(baseline_script, script),
              "samples": []}
    try:
        bundle = ROOT / "assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz"
        for panel in sorted((saved / "published").iterdir()):
            sample = panel.name
            manifest_path = panel / "updated" / (sample + ".harmonizer_manifest.json")
            if not manifest_path.is_file():
                continue
            manifest = json.loads(manifest_path.read_text())
            source = saved / "inputs" / (sample + "_results")
            fasta = saved / "fasta" / (sample + "_results")
            evidence = []
            for tool in ("genomad", "virsorter2", "cenotetaker3", "deep6", "deepmicroclass2",
                         "virbot", "gianthunter", "vicat", "checkv", "tesorter", "vitap", "vcontact3"):
                path = panel / (sample + ".virsorter2_evidence.tsv") if tool == "virsorter2" else source / tool / (sample + "." + tool + "_evidence.tsv")
                if path.is_file():
                    evidence.append(path)
            evidence.append(source / "vicat" / (sample + ".vicat_refined_evidence.tsv"))
            if manifest["input_type"] == "rna":
                evidence.append(source / "vicat" / (sample + ".vicat_orf_evidence.tsv"))
            arguments = ["--sample-id", sample, "--input-type", manifest["input_type"], "--audit-mode", "full",
                         "--ictv-reconciliation", str(bundle), "--vcontact3-min-taxonomy-length",
                         str(manifest["policy"]["vcontact3_min_taxonomy_length"])]
            paths = {"normalized-fasta": fasta / "prep" / (sample + ".normalized.fasta"),
                     "refined-fasta": fasta / "refinement" / (sample + ".refined_candidates.fasta"),
                     "header-map": source / "prep" / (sample + ".header_map.tsv"),
                     "discovery-gate": source / "discovery_gate" / (sample + ".discovery_gate.tsv"),
                     "region-map": source / "refinement" / (sample + ".provirus_region_map.tsv"),
                     "ictv-msl": ROOT / "assets/ICTV_VMR_MSL41.csv"}
            for name, path in paths.items():
                arguments.extend(["--" + name, str(path)])
            arguments.extend(["--evidence", *map(str, evidence), "--vcontact3-groups",
                              str(source / "vcontact3" / (sample + ".vcontact3_group_membership.tsv"))])
            if manifest["policy"]["rna_pair_homology_floor"]:
                arguments.append("--rna-pair-homology-floor")
            dest = out / sample
            baseline = dest / "pre_cap"
            launch(baseline_script, arguments, baseline, sample)
            old_manifest = verify_manifest(baseline, sample)
            species = dest / "species"
            launch(script, arguments, species, sample, "species")
            verify_manifest(species, sample)
            for item in old_manifest["outputs"].values():
                if content(baseline / item["file"]) != content(species / item["file"]):
                    raise AssertionError("Species setting does not reproduce pre-cap output: " + item["file"])
            item = {"sample": sample, "species_restores_pre_cap_outputs": True}
            for rank in ("family", "genus"):
                directory = dest / rank
                launch(script, arguments, directory, sample, rank)
                verify_manifest(directory, sample)
                count, changed = verify_cap(baseline, directory, sample, rank)
                item.update({"candidate_records": count, rank + "_capped_candidates": changed})
            if sample in ("BMARK_DNA_CELLULAR", "BMARK_RNA_VIRAL"):
                parallel = dest / "family_parallel"
                launch(script, arguments, parallel, sample, "family", args.workers)
                left, right = verify_manifest(dest / "family", sample), verify_manifest(parallel, sample)
                if left["outputs"].keys() != right["outputs"].keys():
                    raise AssertionError("Parallel output set differs")
                for entry in left["outputs"].values():
                    if content(dest / "family" / entry["file"]) != content(parallel / entry["file"]):
                        raise AssertionError("Parallel output differs: " + entry["file"])
                item["serial_parallel_outputs_identical"] = True
                item["parallel_workers_requested"] = args.workers
            report["samples"].append(item)
            print("PASS " + sample + ": retention, uncapped audits, family/genus masks, species parity", flush=True)
        if not report["samples"]:
            raise AssertionError("No saved panels found")
        if any(sha(Path(name)) != digest for name, digest in protected.items()):
            raise AssertionError("Saved inputs changed")
        report.update(status="PASS", saved_inputs_unchanged=True, saved_files_verified=len(protected),
                      candidate_records=sum(item["candidate_records"] for item in report["samples"]),
                      viral_decisions_and_sequences_identical=True, review_membership_and_reasons_identical=True,
                      original_evidence_and_uncapped_consensus_preserved=True)
    except Exception as exc:
        report.update(status="FAIL", error=repr(exc))
        raise
    finally:
        (out / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--saved-root", type=Path, required=True)
    parser.add_argument("--baseline-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, choices=range(1, 33), default=8)
    run(parser.parse_args())
