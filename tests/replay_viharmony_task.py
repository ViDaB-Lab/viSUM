#!/usr/bin/env python3
"""Replay one completed viHARMONY task without invoking Nextflow or other tools."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import run_viharmony as harmony

PATH_ARGUMENTS = ("normalized_fasta", "header_map", "discovery_gate", "refined_fasta", "region_map", "ictv_msl")


def digest(path, compressed=False):
    opener = gzip.open if compressed else open
    result = hashlib.sha256()
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def output_manifest(directory, sample):
    flat = directory / f"{sample}.harmonizer_manifest.json"
    return flat if flat.is_file() else directory / "audit" / flat.name


def signature(directory, sample, allow_inputs=False):
    manifest_path = output_manifest(directory, sample)
    manifest = json.loads(manifest_path.read_text())
    result = {}
    for output in manifest["outputs"].values():
        name = output["file"]
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or not relative.name.startswith(sample + "."):
            raise ValueError(f"Unexpected output filename: {name}")
        path = (directory / name).resolve()
        if directory.resolve() not in path.parents:
            raise ValueError(f"Output escapes bundle: {name}")
        if digest(path) != output["sha256"]:
            raise ValueError(f"Output does not match its manifest: {path}")
        result[name] = digest(path, compressed=True) if name.endswith(".gz") else output["sha256"]
        output["sha256"] = result[name]  # Ignore only gzip container timestamps.
    manifest_name = manifest_path.relative_to(directory).as_posix()
    actual = {path.relative_to(directory).as_posix() for path in directory.rglob(sample + ".*") if path.is_file()}
    if not allow_inputs and actual != set(result) | {manifest_name}:
        raise ValueError(f"Unexpected or missing output files in {directory}")
    result[manifest_name] = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    return result


def find_task(work_root, published, sample):
    expected_hash = digest(published)
    candidates = sorted([*work_root.glob(f"*/*/{sample}.harmonizer_manifest.json"),
                         *work_root.glob(f"*/*/audit/{sample}.harmonizer_manifest.json")])
    for manifest in candidates:
        task = manifest.parent.parent if manifest.parent.name == "audit" else manifest.parent
        exitcode = task / ".exitcode"
        if (digest(manifest) == expected_hash and exitcode.is_file()
                and exitcode.read_text().strip() == "0" and (task / ".command.sh").is_file()):
            return task
    raise ValueError("No successful saved task matches the published manifest. No replay launched.")


def task_arguments(task):
    # Parse the recorded command as data; never execute its shell script.
    source = (task / ".command.sh").read_text().replace("\\\n", " ")
    tokens = shlex.split(source, comments=True)
    if (len(tokens) < 3 or not re.fullmatch(r"python(?:\d+(?:\.\d+)*)?", Path(tokens[0]).name)
            or Path(tokens[1]).name != "run_viharmony.py"):
        raise ValueError("Saved command is not a single recognized viHARMONY Python invocation")
    previous = sys.argv
    try:
        sys.argv = [str(ROOT / "bin/run_viharmony.py"), *tokens[2:]]
        return harmony.parse_args()
    finally:
        sys.argv = previous


def task_python(task, override=None):
    if override:
        interpreter = override.resolve()
    else:
        source = (task / ".command.run").read_text()
        match = re.search(r'/bin/activate\s+(\"[^\"]+\"|\'[^\']+\'|/[^\s;]+)', source)
        if not match:
            raise ValueError("Cannot locate the saved Conda environment; supply --python /path/to/its/bin/python")
        interpreter = Path(shlex.split(match.group(1))[0]) / "bin/python"
    if not interpreter.is_file():
        raise FileNotFoundError(f"Saved Python interpreter is unavailable: {interpreter}")
    return interpreter


def input_paths(args, task):
    paths = [getattr(args, name) for name in PATH_ARGUMENTS] + args.evidence + args.vcontact3_groups
    paths = list(dict.fromkeys(path if path.is_absolute() else task / path for path in paths))
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Missing or broken staged input: {path}")
    return paths


def command(args, interpreter, destination, workers):
    values = vars(args).copy()
    values.update(output_prefix=destination / args.sample_id, workers=workers)
    result = [str(interpreter), "-B", str(ROOT / "bin/run_viharmony.py")]
    for name, value in values.items():
        flag = "--" + name.replace("_", "-")
        if isinstance(value, bool):
            if value:
                result.append(flag)
        elif isinstance(value, list):
            result.extend([flag, *map(str, value)])
        else:
            result.extend([flag, str(value)])
    return result


def process_tree_rss(pid):
    """Linux summed resident memory, including children started by any thread."""
    if not sys.platform.startswith("linux"):
        return None, None, 0
    pending, seen, resident, errors = [pid], set(), 0, 0
    page_size = os.sysconf("SC_PAGE_SIZE")
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        try:
            resident += int(Path(f"/proc/{current}/statm").read_text().split()[1]) * page_size
            for children in Path(f"/proc/{current}/task").glob("*/children"):
                pending.extend(map(int, children.read_text().split()))
        except (FileNotFoundError, ProcessLookupError):
            pass
        except PermissionError:
            errors += 1
    return resident, len(seen), errors


def stop_replay(process):
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    elif process.poll() is None:
        process.terminate()
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGKILL)
        else:
            process.kill()
        process.wait()


def execute(argv, task, destination):
    destination.mkdir()
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0")
    peak, processes, errors = None, None, 0
    started = time.perf_counter()
    with (destination / "stdout.log").open("x") as stdout, (destination / "stderr.log").open("x") as stderr:
        process = subprocess.Popen(argv, cwd=task, stdout=stdout, stderr=stderr, env=environment,
                                   start_new_session=os.name == "posix")
        try:
            while process.poll() is None:
                rss, count, denied = process_tree_rss(process.pid)
                if rss is not None:
                    peak, processes = max(peak or 0, rss), max(processes or 0, count)
                errors += denied
                time.sleep(0.2)
            code = process.wait()
        except BaseException:
            stop_replay(process)
            raise
        if code:
            stop_replay(process)
            raise RuntimeError(f"viHARMONY exited {code}; inspect {destination / 'stderr.log'}")
    return dict(wall_seconds=time.perf_counter() - started, peak_sampled_tree_rss_bytes=peak,
                peak_observed_processes=processes, memory_sampling_access_errors=errors,
                command=argv, stdout=str(destination / "stdout.log"), stderr=str(destination / "stderr.log"))


def source_hashes():
    files = ["bin/run_viharmony.py", "bin/evidence_schema.py", "tests/replay_viharmony_task.py",
             "bin/te_routing.py", "bin/ictv_reconciliation.py", "bin/reporting_bundle.py"]
    hashes = {name: digest(ROOT / name) for name in files}
    manifest = ROOT / "bundle_manifest.json"
    if manifest.is_file() and hashes != json.loads(manifest.read_text())["files"]:
        raise ValueError("Test bundle files do not match their recorded checksums")
    return hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--work-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", nargs="+", type=int, default=[1, 12])
    parser.add_argument("--python", type=Path)
    options = parser.parse_args()
    if min(options.workers) < 1 or len(set(options.workers)) != len(options.workers):
        parser.error("Worker counts must be positive and unique")
    results, work_root = options.results.resolve(), options.work_root.resolve()
    if not results.name.endswith("_results"):
        parser.error("--results must name the sample's *_results directory")
    sample = results.name[:-8]
    published = results / "viharmony"
    manifest_path = output_manifest(published, sample)
    manifest = json.loads(manifest_path.read_text())
    task = find_task(work_root, manifest_path, sample)
    args = task_arguments(task)
    for field in ("sample_id", "input_type", "audit_mode"):
        if getattr(args, field) != manifest[field]:
            raise ValueError(f"Saved command and published manifest disagree on {field}")
    if (args.vcontact3_min_taxonomy_length != manifest["policy"]["vcontact3_min_taxonomy_length"]
            or args.rna_pair_homology_floor != manifest["policy"]["rna_pair_homology_floor"]
            or (args.input_type == "rna" and not args.disable_rna_checkv_exception)
                != manifest["policy"]["rna_checkv_conflict_exception"]):
        raise ValueError("Saved command and published manifest disagree on decision settings")
    expected = signature(published, sample)
    if signature(task, sample, allow_inputs=True) != expected:
        raise ValueError("Saved task outputs differ from the published baseline")
    interpreter = task_python(task, options.python)
    version = subprocess.check_output([str(interpreter), "-c", "import sys; print(sys.version); assert sys.version_info >= (3, 11)"], text=True).strip()
    inputs = {str(path): digest(path) for path in input_paths(args, task)}
    sources = source_hashes()
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = dict(status="RUNNING", sample=sample, task=str(task), published=str(published), python=version,
                  interpreter=str(interpreter), source_sha256=sources, input_sha256=inputs, baseline_signature=expected,
                  timing_scope="Elapsed subprocess time including interpreter/pool startup, input loading, output writing and shutdown",
                  memory_scope="Linux peak sampled sum of parent/descendant RSS every 200 ms; shared pages can be counted multiple times; unavailable on other platforms",
                  runs=[])
    report_path = output / "report.json"
    def save():
        report_path.write_text(json.dumps(report, indent=2) + "\n")
    save()
    print(f"Matching completed task: {task}\nPython: {interpreter}\nReport: {report_path}", flush=True)
    try:
        for workers in options.workers:
            destination = output / f"workers_{workers}"
            print(f"Starting {workers} worker(s); logs: {destination}", flush=True)
            measurement = execute(command(args, interpreter, destination, workers), task, destination)
            actual = signature(destination, sample)
            measurement.update(workers=workers, output_equality="PASS" if actual == expected else "FAIL",
                               differing_files=sorted(name for name in set(expected) | set(actual) if expected.get(name) != actual.get(name)))
            report["runs"].append(measurement)
            if source_hashes() != sources or any(digest(Path(path)) != sha for path, sha in inputs.items()):
                raise ValueError("A test source or input changed during replay")
            save()
            if actual != expected:
                raise ValueError(f"Output differences: {measurement['differing_files']}")
            print(f"PASS workers={workers} seconds={measurement['wall_seconds']:.2f} peak_tree_rss_bytes={measurement['peak_sampled_tree_rss_bytes']}", flush=True)
        if signature(published, sample) != expected or signature(task, sample, allow_inputs=True) != expected:
            raise ValueError("Original outputs changed during the replay")
        report["status"] = "PASS"
    except BaseException as error:
        report.update(status="FAIL", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        save()
    print(f"All outputs match the original results. Report: {report_path}", flush=True)


if __name__ == "__main__":
    main()
