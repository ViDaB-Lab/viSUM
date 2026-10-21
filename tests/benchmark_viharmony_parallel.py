"""Developer-only wall-time/RSS comparison using saved replay inputs.

Requires psutil in the test environment, not in the viHARMONY runtime.
"""

import argparse
import contextlib
import gc
import hashlib
import importlib.util
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

import benchmark_viharmony_replay as replay


def execute(job_path, destination, workers, baseline):
    if baseline:
        spec = importlib.util.spec_from_file_location("_indexed_baseline", baseline)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        source = baseline
    else:
        import run_viharmony as module
        source = Path(module.__file__)
    job = json.loads(job_path.read_text())
    args = argparse.Namespace(**job["arguments"])
    for key in ("normalized_fasta", "refined_fasta", "header_map", "discovery_gate", "region_map", "ictv_msl"):
        setattr(args, key, Path(getattr(args, key)))
    for key in ("evidence", "vcontact3_groups"):
        setattr(args, key, [Path(path) for path in getattr(args, key)])
    destination.mkdir(parents=True, exist_ok=False)
    args.output_prefix = destination / args.sample_id
    args.workers = workers
    gc.collect()
    with contextlib.redirect_stdout(io.StringIO()):
        start = time.perf_counter()
        module.run(args)
        elapsed = time.perf_counter() - start
    print(json.dumps(dict(wall_seconds=elapsed, source_sha256=replay.digest(source))))


def monitored_run(command):
    import psutil
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    watched = psutil.Process(process.pid)
    peak_rss = peak_processes = memory_errors = 0
    try:
        while process.poll() is None:
            try:
                tree = [watched, *watched.children(recursive=True)]
            except psutil.NoSuchProcess:
                break
            rss = 0
            for item in tree:
                try:
                    rss += item.memory_info().rss
                except psutil.NoSuchProcess:
                    pass
                except psutil.AccessDenied:
                    memory_errors += 1
            peak_rss = max(peak_rss, rss)
            peak_processes = max(peak_processes, len(tree))
            time.sleep(0.1)
        stdout, stderr = process.communicate()
        if process.returncode:
            raise RuntimeError(f"Replay failed ({process.returncode}): {stderr}")
        result = json.loads(stdout)
        result.update(peak_sampled_tree_rss_bytes=peak_rss, peak_observed_processes=peak_processes,
                      memory_sampling_access_errors=memory_errors)
        return result
    finally:
        if process.poll() is None:
            try:
                for item in watched.children(recursive=True):
                    try:
                        item.terminate()
                    except psutil.NoSuchProcess:
                        pass
            except psutil.NoSuchProcess:
                pass
            process.terminate()
            process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", nargs="+", type=int, default=[1, 2, 4, 12])
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, choices=[0, 1], default=1)
    parser.add_argument("--child-job", type=Path)
    args = parser.parse_args()
    if args.child_job:
        execute(args.child_job, args.output, args.workers[0], args.baseline)
        return
    import psutil
    if not args.jobs or min(args.workers) < 1 or args.repeats < 1:
        parser.error("Provide saved --jobs and positive worker/repeat counts")
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(platform=platform.platform(), python=sys.version, logical_cpus=psutil.cpu_count(),
                  physical_cpus=psutil.cpu_count(logical=False), repeats=args.repeats, warmups=args.warmups,
                  timing_scope="run(args), including worker startup/shutdown and all input/output; excluding interpreter/module startup",
                  memory_scope="Peak sampled sum of parent and descendant RSS at 100 ms intervals; shared pages may be counted multiple times; not a whole-dataset memory forecast",
                  current_source_sha256=replay.digest(replay.ROOT / "bin/run_viharmony.py"), samples=[])
    variants = [(f"workers_{workers}", workers, None) for workers in args.workers]
    if args.baseline:
        variants.insert(0, ("indexed_baseline", 1, args.baseline))
        report["baseline_sha256"] = replay.digest(args.baseline)
    for job_path in args.jobs:
        job = json.loads(job_path.read_text())
        sample = job["arguments"]["sample_id"]
        sample_dir = args.output / sample
        sample_dir.mkdir()
        expected = replay.output_signature(job_path.parent / "0_reference", sample)
        for record in job["inputs"]:
            assert replay.digest(Path(record["path"])) == record["sha256"], "Saved input changed"
        runs = {name: [] for name, _, _ in variants}
        for trial in range(args.repeats + args.warmups):
            order = variants[trial % len(variants):] + variants[:trial % len(variants)]
            for name, workers, baseline in order:
                scratch = contextlib.nullcontext(sample_dir) if trial == 0 else tempfile.TemporaryDirectory(prefix="timed-", dir=sample_dir)
                with scratch as temporary_root:
                    destination = Path(temporary_root) / name
                    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--child-job", str(job_path),
                               "--output", str(destination), "--workers", str(workers)]
                    if baseline:
                        command += ["--baseline", str(baseline)]
                    measurement = monitored_run(command)
                    expected_source = report["baseline_sha256"] if baseline else report["current_source_sha256"]
                    assert measurement["source_sha256"] == expected_source, "Source changed during benchmark"
                    assert replay.output_signature(destination, sample) == expected, f"Output difference in {sample}/{name}"
                if trial >= args.warmups:
                    runs[name].append(measurement)
                print(f"PASS {sample} {name} trial={trial} seconds={measurement['wall_seconds']:.3f} peak_tree_MiB={measurement['peak_sampled_tree_rss_bytes'] / 2**20:.1f}", flush=True)
        summary = dict(sample=sample, inputs=job["normalized_sequences"], candidates=job["candidates"], job=str(job_path),
                       output_equality="PASS", runs=runs,
                       median_wall_seconds={name: statistics.median(r["wall_seconds"] for r in values) for name, values in runs.items()},
                       maximum_sampled_tree_rss_mib={name: max(r["peak_sampled_tree_rss_bytes"] for r in values) / 2**20 for name, values in runs.items()})
        report["samples"].append(summary)
        (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({k: v for k, v in summary.items() if k != "runs"}), flush=True)


if __name__ == "__main__":
    main()
