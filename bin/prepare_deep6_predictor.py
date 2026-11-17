#!/usr/bin/env python3
"""Prepare a task-local Deep6 predictor with exact length-model selection.

The pinned external installation and pretrained weights are never modified.
Compatible with Deep6's existing Python 3.6 environment.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

POLICY = "exact_length_prefix_v1"
SOURCE_SHA256 = {
    "deep6.py": "eb24e6693f16275e9281c15d61bd243d51272403f3045e231c0eddaead33d3a9",
    "deep6_functions.py": "7f9535a84351d4c1d262acd7b774c6404b38fad78d8f707c2993452c22008f90",
}
OLD = b"and clength in x and"
NEW = b"and x.startswith('model_' + clength + '_') and"
LENGTHS = (250, 500, 1000, 1500)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def select_models(directory):
    root = Path(directory).resolve(strict=True)
    result = {}
    for length in LENGTHS:
        matches = sorted(p for p in root.glob("model_{}_*.h5".format(length)) if p.is_file())
        if len(matches) != 1:
            raise ValueError("Expected exactly one model_{}_*.h5 in {}; found {}."
                             .format(length, root, len(matches)))
        path = matches[0]
        with path.open("rb") as handle:
            if handle.read(8) != b"\x89HDF\r\n\x1a\n":
                raise ValueError("Invalid HDF5 model: {}".format(path))
        result[str(length)] = path
    return result


def outside(path, roots):
    resolved = Path(path).resolve()
    for root in roots:
        if os.path.commonpath([str(resolved), str(root)]) == str(root):
            raise ValueError("Refusing to write inside the installed code/models: {}".format(path))
    return resolved


def write_identical_or_new(path, content):
    path = Path(path)
    if path.is_symlink():
        raise ValueError("Refusing a symlink output: {}".format(path))
    if path.exists():
        if not path.is_file() or path.read_bytes() != content:
            raise ValueError("Refusing to overwrite different existing content: {}".format(path))
    else:
        with path.open("xb") as handle:
            handle.write(content)


def prepare(installation, models, output_dir, report):
    installation = Path(installation).resolve(strict=True)
    models = Path(models).resolve(strict=True)
    source = installation / "Master"
    original = {}
    for name, expected in SOURCE_SHA256.items():
        path = source / name
        if digest(path) != expected:
            raise ValueError("Unexpected Deep6 source; refusing an unverified patch: {}".format(path))
        original[name] = path.read_bytes()
    if original["deep6.py"].count(OLD) != 1:
        raise ValueError("Expected one pinned model-selection expression.")
    selected = select_models(models)
    if Path(output_dir).is_symlink():
        raise ValueError("Task-local output directory must not be a symlink.")
    output = outside(output_dir, [installation, models])
    report = outside(report, [installation, models])
    output.mkdir(parents=True, exist_ok=True)
    fixed = original["deep6.py"].replace(OLD, NEW, 1)
    write_identical_or_new(output / "deep6.py", fixed)
    write_identical_or_new(output / "deep6_functions.py", original["deep6_functions.py"])
    provenance = dict(policy=POLICY, original_installation=str(installation),
        source_sha256=SOURCE_SHA256, task_predictor_sha256=hashlib.sha256(fixed).hexdigest(),
        model_directory=str(models),
        models={length: dict(path=str(path), filename=path.name,
                            sha256=digest(path), bytes=path.stat().st_size)
                for length, path in selected.items()})
    write_identical_or_new(report, (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--installation", required=True, type=Path)
    parser.add_argument("--models", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    try:
        prepare(args.installation, args.models, args.output_dir, args.report)
    except (OSError, ValueError) as error:
        raise SystemExit("ERROR: Deep6 model selection: " + str(error))


if __name__ == "__main__":
    main()
