"""Checks for the standalone server replay; no discovery tools are invoked."""

import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import replay_viharmony_task as replay
import test_viharmony as existing


def test_parse_saved_command_preserves_order_and_settings(tmp_path):
    options = ["--sample-id", "sample", "--input-type", "rna", "--audit-mode", "full",
               "--normalized-fasta", "a file.fasta", "--header-map", "header.tsv",
               "--discovery-gate", "gate.tsv", "--refined-fasta", "refined.fasta",
               "--region-map", "regions.tsv", "--ictv-msl", "msl.csv", "--output-prefix", "sample",
               "--evidence", "second.tsv", "first.tsv", "--rna-pair-homology-floor",
               "--disable-rna-checkv-exception", "--vcontact3-min-taxonomy-length", "1234"]
    tokens = ["python3", "/old/bin/run_viharmony.py", *options]
    (tmp_path / ".command.sh").write_text("#!/bin/bash -ue\n# a comment\n" + " \\\n  ".join(map(shlex.quote, tokens)) + "\n")
    parsed = replay.task_arguments(tmp_path)
    assert parsed.evidence == [Path("second.tsv"), Path("first.tsv")]
    assert parsed.normalized_fasta == Path("a file.fasta")
    assert parsed.rna_pair_homology_floor and parsed.disable_rna_checkv_exception
    assert parsed.audit_mode == "full" and parsed.vcontact3_min_taxonomy_length == 1234
    rebuilt = replay.command(parsed, Path(sys.executable), tmp_path, 12)
    assert rebuilt[rebuilt.index("--workers") + 1] == "12"
    assert rebuilt[rebuilt.index("--output-prefix") + 1] == str(tmp_path / "sample")
    assert rebuilt[rebuilt.index("--evidence") + 1:rebuilt.index("--evidence") + 3] == ["second.tsv", "first.tsv"]


def test_unrecognized_shell_content_is_rejected(tmp_path):
    (tmp_path / ".command.sh").write_text("touch must-not-exist; python3 /old/bin/run_viharmony.py --sample-id sample")
    with pytest.raises(ValueError, match="recognized"):
        replay.task_arguments(tmp_path)
    assert not (tmp_path / "must-not-exist").exists()


@pytest.mark.parametrize("quote", ["", "'", '"'])
def test_saved_conda_interpreter_is_found(tmp_path, quote):
    prefix = Path("/home/alex/viSUM/.conda/env-123")
    (tmp_path / ".command.run").write_text(f"source $(conda info --json | awk '...')/bin/activate {quote}{prefix.as_posix()}{quote}\n")
    with patch.object(Path, "is_file", return_value=True):
        assert replay.task_python(tmp_path) == prefix / "bin/python"


def test_task_selection_requires_matching_manifest_and_success(tmp_path):
    published = tmp_path / "published.json"
    published.write_text("published")
    for name, content, exitcode in [("old", "old", "0"), ("failed", "published", "1"), ("match", "published", "0")]:
        task = tmp_path / "work" / "aa" / name
        task.mkdir(parents=True)
        (task / "sample.harmonizer_manifest.json").write_text(content)
        (task / ".exitcode").write_text(exitcode)
        (task / ".command.sh").touch()
    assert replay.find_task(tmp_path / "work", published, "sample").name == "match"
    (tmp_path / "work/aa/match/.exitcode").write_text("1")
    with pytest.raises(ValueError, match="No successful"):
        replay.find_task(tmp_path / "work", published, "sample")


def test_signature_checks_file_content_and_allows_only_explicit_task_inputs(tmp_path):
    content = tmp_path / "sample.data.tsv"
    content.write_text("original\n")
    manifest = dict(outputs={"metadata": dict(file=content.name, sha256=replay.digest(content))})
    (tmp_path / "sample.harmonizer_manifest.json").write_text(json.dumps(manifest))
    expected = replay.signature(tmp_path, "sample")
    (tmp_path / "sample.input.fasta").write_text(">input\nAAAA\n")
    assert replay.signature(tmp_path, "sample", allow_inputs=True) == expected
    with pytest.raises(ValueError, match="Unexpected or missing"):
        replay.signature(tmp_path, "sample")
    content.write_text("tampered\n")
    with pytest.raises(ValueError, match="does not match"):
        replay.signature(tmp_path, "sample", allow_inputs=True)


def test_linux_memory_sampler_includes_children_from_all_threads_once(monkeypatch):
    records = {
        "/proc/10/statm": "10 3", "/proc/20/statm": "10 7", "/proc/30/statm": "20 11",
        "/proc/10/task/10/children": "20", "/proc/10/task/11/children": "20 30",
    }
    def read(path):
        return records[path.as_posix()]
    def children(path, pattern):
        assert pattern == "*/children"
        return [Path("/proc/10/task/10/children"), Path("/proc/10/task/11/children")] if path.as_posix() == "/proc/10/task" else []
    monkeypatch.setattr(replay.sys, "platform", "linux")
    monkeypatch.setattr(replay.os, "sysconf", lambda key: 4096, raising=False)
    with patch.object(Path, "read_text", read), patch.object(Path, "glob", children):
        assert replay.process_tree_rss(10) == (21 * 4096, 3, 0)


def test_complete_replay_cli_matches_fasta_tables_and_gzip_audits(tmp_path):
    original_run = existing.run_viharmony.run
    exercised = []
    def exercise(args):
        original_run(args)
        original = args.output_prefix.parent
        task = tmp_path / "work/aa/completed"
        published = tmp_path / "sample_results/viharmony"
        task.mkdir(parents=True)
        published.mkdir(parents=True)
        manifest_path = replay.output_manifest(original, "sample")
        manifest = json.loads(manifest_path.read_text())
        for name in [manifest_path.relative_to(original).as_posix(), *[entry["file"] for entry in manifest["outputs"].values()]]:
            (task / name).parent.mkdir(parents=True, exist_ok=True)
            (published / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original / name, task / name)
            shutil.copyfile(original / name, published / name)
        (task / "sample.staged_input.fasta").write_text(">staged\nAAAA\n")
        recorded = replay.command(args, Path("python3"), original, 1)
        recorded[1:3] = ["/old/install/bin/run_viharmony.py"]  # Remove -B; obsolete installation is never executed.
        (task / ".command.sh").write_text("# Saved command\n" + " \\\n  ".join(map(shlex.quote, recorded)) + "\n")
        (task / ".exitcode").write_text("0\n")
        before = replay.signature(published, "sample")
        result = subprocess.run([sys.executable, "-B", str(Path(replay.__file__)),
                                 "--results", str(published.parent), "--work-root", str(tmp_path / "work"),
                                 "--output", str(tmp_path / "replay"), "--python", sys.executable,
                                 "--workers", "1", "2"], capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        report = json.loads((tmp_path / "replay/report.json").read_text())
        assert report["status"] == "PASS"
        assert [run["output_equality"] for run in report["runs"]] == ["PASS", "PASS"]
        assert replay.signature(published, "sample") == before
        assert replay.signature(task, "sample", allow_inputs=True) == before
        exercised.append(True)
    with patch.object(existing.run_viharmony, "run", side_effect=exercise):
        method = "test_combines_parent_evidence_region_taxonomy_and_original_ids"
        getattr(existing.ViharmonyTests(method), method)()
    assert exercised == [True]


def test_finisher_accepts_incomplete_original_task_and_preserves_it(tmp_path):
    original_run = existing.run_viharmony.run
    exercised = []
    def exercise(args):
        original_run(args)
        task = tmp_path / "task"
        task.mkdir()
        recorded = replay.command(args, Path("python3"), args.output_prefix.parent, 1)
        recorded[1:3] = ["/old/install/bin/run_viharmony.py"]
        (task / ".command.sh").write_text("# Saved task\n" + shlex.join(recorded) + "\n")
        (task / ".command.run").write_text("# environment supplied by --python in this fixture\n")
        (task / ".command.out").write_text("old process still running\n")
        # No .exitcode or output manifest: the original task is incomplete.
        before = {file.name: replay.digest(file) for file in task.iterdir()}
        finisher = Path(replay.__file__).with_name("finish_viharmony_tasks.py")
        result = subprocess.run([sys.executable, "-B", str(finisher), "--tasks", str(task),
                                 "--output", str(tmp_path / "finished"), "--workers", "2", "--python", sys.executable],
                                capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        report = json.loads((tmp_path / "finished/finish_report.json").read_text())
        assert report["status"] == "FINISHED"
        assert report["runs"][0]["output_integrity"] == "PASS"
        assert report["runs"][0]["baseline_comparison"] == "NOT_PERFORMED"
        assert replay.signature(tmp_path / "finished/sample", "sample") == replay.signature(args.output_prefix.parent, "sample", allow_inputs=True)
        assert {file.name: replay.digest(file) for file in task.iterdir()} == before
        exercised.append(True)
    with patch.object(existing.run_viharmony, "run", side_effect=exercise):
        method = "test_combines_parent_evidence_region_taxonomy_and_original_ids"
        getattr(existing.ViharmonyTests(method), method)()
    assert exercised == [True]
