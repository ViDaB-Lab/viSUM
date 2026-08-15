#!/usr/bin/env python3
"""Run VITAP assignment while enforcing one CPU budget for both DIAMOND calls."""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace


UPSTREAM_THREAD_LITERAL = '"--threads", "10", "--quiet"'
PATCHED_THREAD_LITERAL = '"--threads", threads, "--quiet"'


def patch_assignment_source(source: str) -> str:
    """Replace VITAP 1.12's one hard-coded UniRef90 thread count."""
    occurrences = source.count(UPSTREAM_THREAD_LITERAL)
    if occurrences != 1:
        raise RuntimeError(
            "VITAP 1.12 source compatibility check failed: expected exactly one "
            "hard-coded UniRef90 DIAMOND thread argument, found "
            f"{occurrences}. The upstream implementation may have changed."
        )
    return source.replace(UPSTREAM_THREAD_LITERAL, PATCHED_THREAD_LITERAL, 1)


def load_patched_assignment():
    source_path = None
    spec = importlib.util.find_spec("VITAP_assignment")
    if spec is not None and spec.origin is not None:
        source_path = Path(spec.origin).resolve()

    # Bioconda installs VITAP's companion modules beside the VITAP executable
    # rather than necessarily exposing them as importable site-packages.
    vitap_executable = shutil.which("VITAP")
    if source_path is None and vitap_executable:
        candidate = Path(vitap_executable).resolve().parent / "VITAP_assignment.py"
        if candidate.is_file():
            source_path = candidate

    if source_path is None:
        raise RuntimeError("Could not locate the installed VITAP_assignment module.")

    source = source_path.read_text(encoding="utf-8")
    patched_source = patch_assignment_source(source)

    # VITAP_assignment imports its companion module as a top-level module.
    source_parent = str(source_path.parent)
    if source_parent not in sys.path:
        sys.path.insert(0, source_parent)

    temporary_directory = tempfile.TemporaryDirectory(prefix="visum_vitap_")
    patched_path = Path(temporary_directory.name) / "VITAP_assignment.py"
    patched_path.write_text(patched_source, encoding="utf-8")

    patched_spec = importlib.util.spec_from_file_location(
        "visum_patched_VITAP_assignment",
        patched_path,
    )
    if patched_spec is None or patched_spec.loader is None:
        temporary_directory.cleanup()
        raise RuntimeError("Could not load the patched VITAP assignment module.")

    module = importlib.util.module_from_spec(patched_spec)
    patched_spec.loader.exec_module(module)
    return module.assignment, temporary_directory


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--low-conf", action="store_true")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    return args


def main() -> None:
    args = parse_args()
    assignment, temporary_directory = load_patched_assignment()
    try:
        assignment(
            SimpleNamespace(
                fasta=args.fasta,
                db=args.db,
                out=args.out,
                cpu=args.threads,
                low_conf=args.low_conf,
            )
        )
    finally:
        temporary_directory.cleanup()


if __name__ == "__main__":
    main()
