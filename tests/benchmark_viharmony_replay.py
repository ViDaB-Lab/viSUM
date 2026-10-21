"""Time frozen and indexed viHARMONY on identical archived pipeline inputs.

Standalone, standard-library-only developer benchmark; does not invoke Nextflow
or discovery tools. Each replay uses a fresh process and separate output folder.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gc
import gzip
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import tempfile
import time
import types

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = "d4b90f3c235eae41b03a24ea4fbfd62b5febdb1c"
sys.path.insert(0, str(ROOT / "bin"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_for(version):
    if version == "reference":
        return subprocess.check_output(["git", "show", f"{REFERENCE}:bin/run_viharmony.py"], cwd=ROOT)
    return (ROOT / "bin/run_viharmony.py").read_bytes()


def worker(job, version, destination):
    source = source_for(version)
    if version == "indexed":
        # Spawned worker functions/classes must belong to an importable module.
        module = importlib.import_module("run_viharmony")
    else:
        module = types.ModuleType("_timed_viharmony")
        sys.modules[module.__name__] = module
        exec(compile(source, f"<{version} viHARMONY>", "exec"), module.__dict__)
    args = argparse.Namespace(**job["arguments"])
    for key in ("normalized_fasta", "refined_fasta", "header_map", "discovery_gate", "region_map", "ictv_msl"):
        setattr(args, key, Path(getattr(args, key)))
    for key in ("evidence", "vcontact3_groups"):
        setattr(args, key, [Path(value) for value in getattr(args, key)])
    destination.mkdir(parents=True, exist_ok=False)
    args.output_prefix = destination / args.sample_id
    gc.collect()
    with contextlib.redirect_stdout(io.StringIO()):
        wall_start, cpu_start = time.perf_counter(), time.process_time()
        module.run(args)
        wall, cpu = time.perf_counter() - wall_start, time.process_time() - cpu_start
    return dict(wall_seconds=wall, cpu_seconds=cpu,
                source_sha256=hashlib.sha256(source).hexdigest())


def output_signature(directory, sample):
    result = {}
    for path in sorted(directory.glob(sample + ".*")):
        data = path.read_bytes()
        if path.name.endswith(".harmonizer_manifest.json"):
            manifest = json.loads(data)
            for output in manifest["outputs"].values():
                raw = (directory / output["file"]).read_bytes()
                if hashlib.sha256(raw).hexdigest() != output["sha256"]:
                    raise AssertionError(f"Invalid manifest hash: {output['file']}")
                if output["file"].endswith(".gz"):
                    output["sha256"] = hashlib.sha256(gzip.decompress(raw)).hexdigest()
            data = json.dumps(manifest, sort_keys=True).encode()
        elif path.suffix == ".gz":
            data = gzip.decompress(data)
        result[path.name] = hashlib.sha256(data).hexdigest()
    if len(result) < 14:
        raise AssertionError(f"Incomplete output set: {directory}")
    return result


def fasta_count(path):
    with path.open("rb") as handle:
        return sum(line.startswith(b">") for line in handle)


def tsv_count(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def prepare_job(directory, audit_mode):
    sample = directory.name.removesuffix("_results")
    archived = json.loads((directory / "viharmony" / f"{sample}.harmonizer_manifest.json").read_text())
    input_type = archived["input_type"]
    shared = ["genomad", "virsorter2", "cenotetaker3", "vicat", "checkv", "tesorter", "vitap", "vcontact3"]
    specific = ["deep6", "virbot"] if input_type == "rna" else ["deepmicroclass2", "gianthunter"]
    evidence = [directory / tool / f"{sample}.{tool}_evidence.tsv" for tool in shared + specific]
    evidence.extend(directory / "vicat" / f"{sample}.{suffix}.tsv" for suffix in ["vicat_context", "vicat_refined_evidence"])
    if input_type == "rna":
        evidence.append(directory / "vicat" / f"{sample}.vicat_orf_evidence.tsv")
    evidence.sort(key=lambda path: path.name)
    args = dict(
        sample_id=sample, input_type=input_type,
        normalized_fasta=directory / "prep" / f"{sample}.normalized.fasta",
        header_map=directory / "prep" / f"{sample}.header_map.tsv",
        discovery_gate=directory / "discovery_gate" / f"{sample}.discovery_gate.tsv",
        refined_fasta=directory / "refinement" / f"{sample}.refined_candidates.fasta",
        region_map=directory / "refinement" / f"{sample}.provirus_region_map.tsv",
        ictv_msl=ROOT / "assets/ICTV_VMR_MSL41.csv",
        evidence=evidence,
        vcontact3_groups=[directory / "vcontact3" / f"{sample}.vcontact3_group_membership.tsv"],
        audit_mode=audit_mode,
        rna_pair_homology_floor=archived.get("policy", {}).get("rna_pair_homology_floor", False),
        vcontact3_min_taxonomy_length=archived.get("policy", {}).get("vcontact3_min_taxonomy_length", 1000),
    )
    inputs = [value for value in args.values() if isinstance(value, Path)] + evidence + args["vcontact3_groups"]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    job = dict(arguments=args, inputs=[dict(path=str(path), sha256=digest(path), bytes=path.stat().st_size) for path in inputs],
               normalized_sequences=fasta_count(args["normalized_fasta"]), candidates=fasta_count(args["refined_fasta"]),
               evidence_rows=sum(tsv_count(path) for path in evidence if not path.name.endswith(".vicat_orf_evidence.tsv")),
               locus_rows=sum(tsv_count(path) for path in evidence if path.name.endswith(".vicat_orf_evidence.tsv")))
    return job


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--audit-mode", choices=["compact", "full", "none"], default="compact")
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--version", choices=["reference", "indexed"])
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(json.loads(args.worker.read_text()), args.version, args.output)))
        return
    if not args.results or args.output is None or args.repeats < 1:
        parser.error("Provide --results, a new --output directory, and positive --repeats")
    jobs = [prepare_job(directory, args.audit_mode) for root in args.results for directory in sorted(root.glob("*_results"))]
    if not jobs:
        parser.error("No archived sample results found")
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(reference_commit=REFERENCE, python=sys.version, platform=platform.platform(),
                  audit_mode=args.audit_mode, repeats=args.repeats, warmups_per_version=1,
                  timing_scope="run(args), including all input loading, decision making, output writing and hashing; excludes interpreter/module startup and comparison",
                  ordering="One warm-up per version, then alternating reference/indexed order in fresh sequential subprocesses; fixed PYTHONHASHSEED=0; OS file cache not flushed",
                  source_sha256={v: hashlib.sha256(source_for(v)).hexdigest() for v in ("reference", "indexed")}, samples=[])
    environment = dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1")
    for job in jobs:
        sample = job["arguments"]["sample_id"]
        sample_dir = args.output / sample
        sample_dir.mkdir()
        job_path = sample_dir / "inputs.json"
        job_path.write_text(json.dumps(job, indent=2, default=str) + "\n")
        record = {key: value for key, value in job.items() if key != "arguments"}
        record.update(sample=sample, arguments=job["arguments"], timings={v: [] for v in ("reference", "indexed")})
        expected = None
        for trial in range(args.repeats + 1):
            versions = ("reference", "indexed") if trial % 2 == 0 else ("indexed", "reference")
            for version in versions:
                # Keep one complete output pair; discard only our own timed
                # scratch outputs after their hashes and contents are checked.
                scratch = contextlib.nullcontext(sample_dir) if trial == 0 else tempfile.TemporaryDirectory(prefix="timed-", dir=sample_dir)
                with scratch as temporary_root:
                    destination = Path(temporary_root) / f"{trial}_{version}"
                    output = subprocess.check_output([sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(job_path),
                                                      "--version", version, "--output", str(destination)], env=environment, text=True)
                    measurement = json.loads(output)
                    assert measurement["source_sha256"] == report["source_sha256"][version], "Source changed during timing"
                    signature = output_signature(destination, sample)
                if expected is None:
                    expected = signature
                if signature != expected:
                    raise AssertionError(f"Output mismatch: {sample}, trial {trial}, {version}")
                if trial:
                    record["timings"][version].append(measurement)
        for input_record in job["inputs"]:
            assert digest(Path(input_record["path"])) == input_record["sha256"], "Input changed during replay"
        record["output_files_compared_per_run"] = len(expected)
        record["output_signature"] = expected
        record["output_equality"] = "PASS"
        record["median_wall_seconds"] = {v: statistics.median(t["wall_seconds"] for t in trials) for v, trials in record["timings"].items()}
        record["speedup"] = record["median_wall_seconds"]["reference"] / record["median_wall_seconds"]["indexed"]
        report["samples"].append(record)
        (args.output / "report.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
        print(json.dumps({key: record[key] for key in ("sample", "normalized_sequences", "candidates", "evidence_rows", "median_wall_seconds", "speedup", "output_equality")}), flush=True)
    print(f"Report: {args.output / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
