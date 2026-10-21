#!/usr/bin/env python3
"""Package an explicit standalone source/data allowlist; no Git/server writes."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = (
    "bin/ictv_reconciliation.py", "bin/build_ictv_reconciliation_bundle.py",
    "bin/run_viharmony.py", "bin/evidence_schema.py", "bin/reconcile_benchmark_truth.py",
    "bin/replay_viharmony_taxonomy.py", "bin/standardize_vicat.py",
    "bin/standardize_vitap.py", "bin/standardize_vcontact3.py",
    "visum_nextflow.nf", "visum.config", "nextflow.config", "modules/local/viharmony.nf",
    "assets/ICTV_VMR_MSL41.csv", "assets/ictv_reconciliation/README.md",
    "assets/ictv_reconciliation/LICENSE.txt",
    "assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz",
    "assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz.manifest.json",
    "tests/test_ictv_reconciliation.py", "tests/test_ictv_reference.py",
    "tests/test_viharmony.py", "tests/test_viharmony_workflow_contracts.py",
    "tests/audit_ictv_unassigned_ranks.py", "tests/validate_ictv_reconciliation_replay.py",
    "tests/run_ictv_reconciliation_validation.py",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def package(results, archive):
    archive = archive.resolve()
    sidecar = archive.with_name(archive.name + ".sha256")
    if archive.exists() or sidecar.exists():
        raise ValueError("Use a new archive path; never overwrite a previous package")
    reference_dir = ROOT / "outputs/ictv_reconciliation_prototype/unassigned_rank_correction_verified"
    reference = json.loads((reference_dir / "benchmark_replay/report.json").read_text(encoding="utf-8"))
    expected_files = {}
    for name, value in reference["input_sha256"].items():
        relative = "/".join(name.replace("\\", "/").split("/")[-3:])
        if relative in expected_files:
            raise AssertionError("Duplicate saved-file identity")
        expected_files[relative] = value
    protected = {name: digest(ROOT / name) for name in reference["protected_source_sha256"]}
    if protected != reference["protected_source_sha256"]:
        raise AssertionError("Protected source differs from validated snapshot")
    for name, value in reference["prototype_source_sha256"].items():
        if digest(ROOT / name) != value:
            raise AssertionError(f"Prototype source differs from validated snapshot: {name}")
    with tempfile.TemporaryDirectory(prefix="visum-ictv-package-") as temporary:
        stage = Path(temporary) / "viSUM-reconciliation"
        stage.mkdir()
        for name in SOURCE_FILES:
            target = stage / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / name, target)
        for name, expected in expected_files.items():
            source = results / name
            if digest(source) != expected:
                raise AssertionError(f"Saved input differs from validated snapshot: {name}")
            target = stage / "saved_results" / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        reports = {"expected_corrected/release_audit.json": reference_dir / "release_audit.json"}
        for name in ("report.json", "taxonomy_changes.jsonl", "reconciliation_audit.jsonl"):
            reports["expected_corrected/" + name] = reference_dir / "benchmark_replay" / name
            reports["expected_previous/" + name] = ROOT / "outputs/ictv_reconciliation_prototype/validated_snapshot" / name
        for name, source in reports.items():
            target = stage / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        manifest = {"prototype_only": True, "target_release": reference["target_release"],
            "bundle": "assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz",
            "expected_unittest_count": 74,
            "test_scope": "dependency-free unittest cases only, not the full 156-test local pytest suite",
            "samples": [row["sample"] for row in reference["samples"]],
            "files": {path.relative_to(stage).as_posix(): digest(path)
                      for path in sorted(stage.rglob("*")) if path.is_file()}}
        (stage / "PACKAGE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "w:gz", compresslevel=6, format=tarfile.PAX_FORMAT) as tar:
            root_info = tarfile.TarInfo("viSUM-reconciliation")
            root_info.type = tarfile.DIRTYPE
            root_info.mode = 0o755
            tar.addfile(root_info)
            for path in sorted(stage.rglob("*")):
                if path.is_dir():
                    info = tarfile.TarInfo("viSUM-reconciliation/" + path.relative_to(stage).as_posix())
                    info.type = tarfile.DIRTYPE
                    info.mode = 0o755
                    tar.addfile(info)
                if path.is_file():
                    info = tar.gettarinfo(str(path), arcname="viSUM-reconciliation/" + path.relative_to(stage).as_posix())
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mode = 0o644
                    with path.open("rb") as handle:
                        tar.addfile(info, handle)
    sidecar.write_text(digest(archive) + "  " + archive.name + "\n", encoding="ascii")
    print(json.dumps({"archive": str(archive), "sha256": digest(archive),
                      "bytes": archive.stat().st_size, "saved_files": len(expected_files)}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    args = parser.parse_args()
    package(args.results_root, args.archive)


if __name__ == "__main__":
    main()
