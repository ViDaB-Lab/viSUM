#!/usr/bin/env python3
"""Finish pending viHARMONY tasks into a fresh directory using recorded inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import replay_viharmony_task as replay


def fingerprints(task, args):
    paths = [task / ".command.sh", task / ".command.run", *replay.input_paths(args, task)]
    return {str(path): replay.digest(path) for path in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tasks", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--python", type=Path)
    options = parser.parse_args()
    if options.workers < 1:
        parser.error("--workers must be a positive integer")
    sources = replay.source_hashes()
    own_hash = replay.digest(Path(__file__).resolve())
    bundle_manifest = replay.ROOT / "bundle_manifest.json"
    if bundle_manifest.is_file():
        expected = json.loads(bundle_manifest.read_text()).get("additional_files", {}).get("tests/finish_viharmony_tasks.py")
        if expected != own_hash:
            raise ValueError("Finisher does not match its bundle checksum")
    jobs, samples = [], set()
    for path in options.tasks:
        task = path.resolve()
        args = replay.task_arguments(task)
        if args.sample_id in samples:
            raise ValueError(f"Multiple task directories refer to sample {args.sample_id}")
        samples.add(args.sample_id)
        interpreter = replay.task_python(task, options.python)
        version = subprocess.check_output([str(interpreter), "-c", "import sys; print(sys.version); assert sys.version_info >= (3, 11)"], text=True).strip()
        jobs.append((task, args, interpreter, version, fingerprints(task, args)))
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = dict(status="RUNNING", workers=options.workers, source_sha256=sources,
                  finisher_sha256=own_hash, runs=[],
                  mode="Finish recorded viHARMONY tasks outside Nextflow; original tasks/results/cache are not written",
                  validation="Output manifest/file checksums and unchanged inputs; original results may be incomplete, so no claim of per-sample baseline equality",
                  memory_scope="Peak sampled sum of parent/descendant RSS every 200 ms; shared pages can be counted multiple times")
    report_path = output / "finish_report.json"
    def save():
        report_path.write_text(json.dumps(report, indent=2) + "\n")
    save()
    print(f"Finish report: {report_path}", flush=True)
    try:
        for task, args, interpreter, version, inputs in jobs:
            if fingerprints(task, args) != inputs or replay.source_hashes() != sources:
                raise ValueError("Recorded inputs or test sources changed before execution")
            destination = output / args.sample_id
            print(f"Starting {args.sample_id}: workers={options.workers}; outputs={destination}", flush=True)
            measurement = replay.execute(replay.command(args, interpreter, destination, options.workers), task, destination)
            signature = replay.signature(destination, args.sample_id)
            if (fingerprints(task, args) != inputs or replay.source_hashes() != sources
                    or replay.digest(Path(__file__).resolve()) != own_hash):
                raise ValueError("Recorded inputs or test sources changed during execution")
            measurement.update(sample=args.sample_id, task=str(task), python=version, interpreter=str(interpreter),
                               input_sha256=inputs, output_signature=signature, output_integrity="PASS",
                               outputs=str(destination), baseline_comparison="NOT_PERFORMED")
            report["runs"].append(measurement)
            save()
            print(f"FINISHED {args.sample_id}: {measurement['wall_seconds']:.2f}s; output checks PASS", flush=True)
        report["status"] = "FINISHED"
    except BaseException as error:
        report.update(status="FAIL", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        save()
    print(f"Final files are ready in {output}. Finish report: {report_path}", flush=True)


if __name__ == "__main__":
    main()
