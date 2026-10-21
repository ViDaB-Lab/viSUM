#!/usr/bin/env python3
"""Private saved-input reporting-cap validation; never deploys or runs searches."""
from __future__ import annotations
import argparse
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import traceback
import unittest

ROOT = Path(__file__).resolve().parents[1]
STAGES = (("first_family", None, "family"), ("unchanged_resume", None, "family"),
          ("genus", "genus", "genus"), ("species", "species", "species"),
          ("restore_family", "family", "family"))
TESTS = ("test_viharmony_reporting_cap.py", "test_rank_cap_linux_validation.py")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(path):
    if path.name.endswith("_audit.tsv.gz"):
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = sorted(reader.fieldnames or [])
            rows = Counter(json.dumps(row, sort_keys=True, ensure_ascii=True) for row in reader)
        value = json.dumps([columns, sorted(rows.items())], separators=(",", ":"))
        return hashlib.sha256(value.encode()).hexdigest()
    return sha(path)


def scientific_snapshot(directory, sample):
    manifest = json.loads((directory / (sample + ".harmonizer_manifest.json")).read_text())
    digests = {}
    for label, record in manifest["outputs"].items():
        path = (directory / record["file"]).resolve()
        if path.parent != directory.resolve() or sha(path) != record["sha256"]:
            raise AssertionError("Output path/checksum mismatch: " + record["file"])
        digests[label] = fingerprint(path)
        record.pop("sha256")
    return dict(manifest=manifest, scientific_digests=digests)


def verify_expected(directory, sample, rank, expected):
    actual = scientific_snapshot(directory, sample)
    wanted = expected[rank][sample]
    if actual != wanted:
        changed = [k for k in actual["scientific_digests"]
                   if actual["scientific_digests"][k] != wanted["scientific_digests"].get(k)]
        raise AssertionError(f"Saved expected outputs differ: {sample}/{rank}/{changed}; policy or output sets may differ")
    return dict(sample=sample, rank=rank, outputs_verified=len(actual["scientific_digests"]),
                retained_sequences_decisions_reviews_and_uncapped_audits_verified=True)


def verify_publication(work, published, sample):
    """Require real published copies to equal the CURRENT task, byte for byte.

    A stale but internally consistent old manifest is not sufficient. This
    verifier is read-only: it never repairs/copies outputs into published/.
    """
    task = scientific_snapshot(work, sample)
    copied = scientific_snapshot(published, sample)
    if copied != task:
        raise AssertionError("Published output bundle is stale or differs from current task: " + sample)
    names = [sample + ".harmonizer_manifest.json"]
    names.extend(record["file"] for record in task["manifest"]["outputs"].values())
    for name in names:
        if sha(work / name) != sha(published / name):
            raise AssertionError("Published copy differs byte-for-byte from current task: " + name)
    return dict(sample=sample, published_files_verified=len(names),
                manifest_sha256=sha(published / names[0]), exact_current_task_copies=True)


def verify_package(project=ROOT):
    project = project.resolve()
    manifest = json.loads((project / "RANK_CAP_PACKAGE.json").read_text())
    for name, digest in manifest["files"].items():
        path = (project / name).resolve()
        if project not in path.parents or not path.is_file() or sha(path) != digest:
            raise AssertionError("Package source/input checksum mismatch: " + name)
    if not re.search(r"harmonizer_max_rank\s*=\s*'family'", (project / "visum.config").read_text()):
        raise AssertionError("Real configuration no longer defaults to family")
    return manifest


def direct_arguments(spec, project=ROOT):
    args = ["--sample-id", spec["sample"], "--input-type", spec["type"],
            "--audit-mode", "full", "--ictv-msl", str(project / "assets/ICTV_VMR_MSL41.csv"),
            "--ictv-reconciliation", str(project / spec["bundle"]),
            "--vcontact3-min-taxonomy-length", str(spec["minimum_length"])]
    for flag, key in (("normalized-fasta", "normalized"), ("refined-fasta", "refined"),
                      ("header-map", "header"), ("discovery-gate", "gate"), ("region-map", "regions")):
        args.extend(["--" + flag, str(project / spec["folder"] / spec[key])])
    args.extend(["--evidence", *[str(project / spec["folder"] / f) for f in spec["evidence"]],
                 "--vcontact3-groups", str(project / spec["folder"] / spec["groups"])])
    if spec["type"] == "rna" and spec["pair_floor"]:
        args.append("--rna-pair-homology-floor")
    return args


def native_command(nextflow, native, config, parameters, stage, rank, session=None, project=ROOT):
    command = [str(nextflow), "-log", str(native / (stage + ".nextflow.log")), "-C", str(config),
               "run", str(project / "native_validation.nf"), "-params-file", str(parameters),
               "-work-dir", str(native / "work"), "-ansi-log", "false",
               "-with-trace", str(native / (stage + ".trace.tsv"))]
    if rank is not None:
        command.extend(["--harmonizer_max_rank", rank])
    if session is not None:
        if not re.fullmatch(r"[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}", session):
            raise ValueError("Invalid explicit resume session")
        command.extend(["-resume", session])
    return command


def verify_trace(rows, stage, samples):
    expected = {(name, sample) for name in ("STAGE_SAVED_INPUTS", "VIHARMONY") for sample in samples}
    seen = set()
    for row in rows:
        match = re.fullmatch(r"(?:\w+:)*(STAGE_SAVED_INPUTS|VIHARMONY) \(([^)]+)\)", row["name"])
        if not match or (match[1], match[2]) not in expected or (match[1], match[2]) in seen:
            raise AssertionError("Unexpected, duplicated or search task: " + str(row))
        name, sample = match.groups()
        wanted = "COMPLETED" if stage == "first_family" or (name == "VIHARMONY" and stage in {"genus", "species"}) else "CACHED"
        if row["status"] != wanted:
            raise AssertionError(f"Incorrect cache result: {stage}/{name}/{sample}: {row['status']} != {wanted}")
        seen.add((name, sample))
    if seen != expected:
        raise AssertionError("Missing staging or harmonizer task")


def verify_task_command(path, rank, workers):
    commands = [line for line in path.read_text().replace("\\\n", " ").splitlines()
                if line.strip().startswith("python3 ")]
    if len(commands) != 1:
        raise AssertionError("Expected one real harmonizer command")
    args = shlex.split(commands[0])
    if (args[args.index("--max-rank") + 1] != rank
            or args[args.index("--workers") + 1] != str(workers)):
        raise AssertionError("Incorrect actual rank/worker argument")


def private_environment(native, installed_home):
    framework = installed_home.resolve() / "framework/26.04.6"
    if not (framework / "nextflow-26.04.6-one.jar").is_file():
        raise ValueError("Offline Nextflow 26.04.6 runtime missing: " + str(framework))
    shutil.copytree(framework, native / "nextflow-home/framework/26.04.6")
    env = dict(os.environ)
    for name in ("NXF_HOME", "NXF_PLUGINS_DIR", "NXF_WORK", "NXF_CACHE_DIR", "NXF_CLOUDCACHE_PATH",
                 "NXF_DISABLE_CACHE", "NXF_PARAMS_FILE", "NXF_OPTS", "NXF_VER"):
        env.pop(name, None)
    env.update(NXF_HOME=str(native / "nextflow-home"), NXF_PLUGINS_DIR=str(native / "plugins"),
               NXF_OFFLINE="true", NXF_VER="26.04.6", NXF_DISABLE_CHECK_LATEST="true",
               NXF_SYNTAX_PARSER="v2", PYTHONDONTWRITEBYTECODE="1", TMPDIR="/tmp", TEMP="/tmp", TMP="/tmp")
    return env


def native_check(output, nextflow, installed_home, specs, expected, workers):
    if platform.system() != "Linux":
        raise ValueError("Native cache validation must run on Linux")
    native = output / "native"
    native.mkdir()
    bin_dir = native / "bin"
    bin_dir.mkdir()
    wrapper = bin_dir / "python3"
    wrapper.write_text("#!/bin/sh\nexec " + shlex.quote(sys.executable) + ' "$@"\n')
    wrapper.chmod(0o755)
    quote = lambda value: "'" + str(value).replace("\\", "/").replace("'", "\\'") + "'"
    config = native / "native.config"
    config.write_text("includeConfig " + quote(ROOT / "visum.config") + "\n"
        "process.executor = 'local'\nprocess.errorStrategy = 'terminate'\nprocess.time = null\n"
        "conda.enabled = false\ndocker.enabled = false\nsingularity.enabled = false\n"
        "report.enabled = false\ntimeline.enabled = false\ndag.enabled = false\n"
        "trace.fields = 'task_id,hash,name,status,exit,workdir'\n"
        "env.PATH = " + quote(bin_dir) + " + ':' + System.getenv('PATH')\n")
    payload = dict(saved_specs=[dict(s, files=[str(ROOT / p) for p in s["files"]]) for s in specs],
        harmonizer_cpus=workers, max_cpus=workers, harmonizer_memory="32 GB", harmonizer_time=None,
        max_memory="64 GB", outdir=str(native / "published"),
        rna_pair_homology_floor=next(s["pair_floor"] for s in specs if s["type"] == "rna"),
        ictv_csv=str(ROOT / "assets/ICTV_VMR_MSL41.csv"), ictv_reconciliation=str(ROOT / specs[0]["bundle"]))
    # Deliberately omit harmonizer_max_rank: the first run must use the actual default.
    parameters = native / "params.json"
    parameters.write_text(json.dumps(payload, indent=2) + "\n")
    env = private_environment(native, installed_home)
    session, stages, first_hashes = None, [], {}
    for stage, override, rank in STAGES:
        command = native_command(nextflow, native, config, parameters, stage, override, session)
        (native / (stage + ".command.txt")).write_text(shlex.join(command) + "\n")
        with (native / (stage + ".console.log")).open("w") as log:
            result = subprocess.run(command, cwd=native, env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise AssertionError("Nextflow failed; see " + str(native / (stage + ".console.log")))
        if session is None:
            matches = re.findall(r"Session UUID:\s*([0-9a-f-]{36})", (native / (stage + ".nextflow.log")).read_text())
            if len(set(matches)) != 1:
                raise AssertionError("Cannot identify the original native cache session")
            session = matches[0]
            if not (native / ".nextflow/cache" / session / "db").is_dir():
                raise AssertionError("Original cache is not in the private launch directory")
        with (native / (stage + ".trace.tsv")).open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        verify_trace(rows, stage, [s["sample"] for s in specs])
        checks, publication_checks = [], []
        receipts = native / (stage + "-receipts")
        receipts.mkdir()
        for row in rows:
            work = Path(row["workdir"]).resolve()
            if (native / "work").resolve() not in work.parents:
                raise AssertionError("Task work directory escaped the private test")
            if "VIHARMONY (" not in row["name"]:
                continue
            sample = row["name"].rsplit("(", 1)[1].rstrip(")")
            verify_task_command(work / ".command.sh", rank, workers)
            if stage == "first_family":
                first_hashes[sample] = row["hash"]
            if stage in {"unchanged_resume", "restore_family"} and row["hash"] != first_hashes[sample]:
                raise AssertionError("Did not return to the original family task hash")
            checks.append(verify_expected(work, sample, rank, expected))
            published = native / "published" / (sample + "_results/viharmony")
            publication_checks.append(verify_publication(work, published, sample))
            verify_expected(published, sample, rank, expected)
            shutil.copyfile(published / (sample + ".harmonizer_manifest.json"),
                            receipts / (sample + ".published.harmonizer_manifest.json"))
            for suffix in ("harmonizer_manifest.json",):
                shutil.copyfile(work / (sample + "." + suffix), receipts / (sample + "." + suffix))
            shutil.copyfile(work / ".command.sh", receipts / (sample + ".command.sh"))
        stages.append(dict(stage=stage, rank=rank, tasks=rows, checks=checks,
                           publication_checks=publication_checks, original_session=session))
        (native / (stage + ".verified.json")).write_text(json.dumps(stages[-1], indent=2) + "\n")
        print("PASS native Nextflow: " + stage, flush=True)
    return dict(status="PASS", original_session=session, stages=stages,
                scope="actual VIHARMONY module plus saved-input copy process; task AND published outputs verified after all five stages",
                published_outputs_verified=True,
                expensive_search_modules_present=False, full_visum_entrypoint_executed=False)


def run(args):
    output = args.output_dir.resolve()
    if output.exists() or ROOT == output or ROOT in output.parents:
        raise ValueError("Use a new output directory outside the extracted package")
    output.mkdir(parents=True)
    status = dict(status="RUNNING", python=sys.version, platform=platform.platform(), deployed=False,
                  student_runs_accessed=False, no_search_programs_launched=True)
    try:
        package = verify_package()
        status["package_manifest_sha256"] = sha(ROOT / "RANK_CAP_PACKAGE.json")
        specs = json.loads((ROOT / "saved_specs.json").read_text())
        expected = json.loads((ROOT / "expected_outputs.json").read_text())
        suite = unittest.TestSuite()
        for pattern in TESTS:
            suite.addTests(unittest.defaultTestLoader.discover(str(ROOT / "tests"), pattern=pattern))
        with (output / "regression.log").open("w") as log, redirect_stdout(log), redirect_stderr(log):
            result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
        status["regression"] = dict(tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors), skipped=len(result.skipped))
        if not result.wasSuccessful() or result.skipped or result.testsRun != package["test_count"]:
            raise AssertionError("Regression checks incomplete or failed; see regression.log")
        print(f"PASS {result.testsRun} dependency-free regression tests", flush=True)
        if args.local_only:
            checks = []
            for rank in ("family", "genus", "species"):
                for spec in specs:
                    directory = output / rank / spec["sample"]
                    directory.mkdir(parents=True)
                    command = [sys.executable, "-B", str(ROOT / "bin/run_viharmony.py"),
                               *direct_arguments(spec), "--workers", str(args.workers),
                               "--output-prefix", str(directory / spec["sample"])]
                    if rank != "family":
                        command.extend(["--max-rank", rank])
                    with (directory / "run.log").open("w") as log:
                        process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
                    if process.returncode:
                        raise AssertionError("Local extracted-package harmonizer failed: " + str(directory))
                    checks.append(verify_expected(directory, spec["sample"], rank, expected))
                    print(f"PASS extracted-package Python replay: {spec['sample']}/{rank}", flush=True)
            status.update(status="LOCAL_ONLY_NATIVE_PENDING", local_checks=checks, native_nextflow={"status": "NOT_RUN"})
        else:
            if not args.nextflow or not args.nextflow_home:
                raise ValueError("Native validation requires --nextflow and --nextflow-home")
            status["native_nextflow"] = native_check(output, args.nextflow.resolve(), args.nextflow_home, specs, expected, args.workers)
            status["status"] = "PASS"
        verify_package()
        status.update(package_and_saved_inputs_unchanged=True, candidate_records=package["candidate_records"],
                      biological_comparison="Exact expected outputs from independently verified pre-cap replay; audit row/column order and gzip timestamps only normalized")
    except Exception:
        status.update(status="FAIL", error=traceback.format_exc())
        raise
    finally:
        if status["status"] == "FAIL" and (output / "native/work").is_dir():
            work_root = (output / "native/work").resolve()
            diagnostics = output / "native/failure-receipts"
            diagnostics.mkdir(exist_ok=True)
            for name in (".command.sh", ".command.err", ".command.out", ".exitcode"):
                for path in work_root.rglob(name):
                    if work_root in path.resolve().parents:
                        target = diagnostics / ("_".join(path.relative_to(work_root).parts))
                        shutil.copyfile(path, target)
        (output / "validation_status.json").write_text(json.dumps(status, indent=2) + "\n")
    print("Final status: " + status["status"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, choices=range(2, 9), default=8)
    parser.add_argument("--nextflow", type=Path)
    parser.add_argument("--nextflow-home", type=Path)
    parser.add_argument("--local-only", action="store_true")
    run(parser.parse_args())
