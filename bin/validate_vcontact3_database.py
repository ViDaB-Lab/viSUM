#!/usr/bin/env python3
"""Validate an official vConTACT3 database manifest and its runtime files."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Iterable


REQUIRED_REFSEQ_FILES = (
    "Genomes Report",
    "gene2genome",
    "proteins",
    "nucleotides",
)
REQUIRED_DOMAINS = ("prokaryotes", "eukaryotes")
REQUIRED_VOG_FILES = ("rank_exclusivity", "host_domain")


def fail(message: str) -> "NoReturn":
    raise ValueError(message)


def newest_manifest(database_path: Path) -> Path:
    if database_path.is_file():
        if database_path.suffix != ".json":
            fail(f"database file is not a JSON manifest: {database_path}")
        return database_path

    if not database_path.is_dir():
        fail(f"database path does not exist or is not readable: {database_path}")

    manifests = [
        path
        for path in database_path.glob("*.json")
        if path.stem.isdigit() and len(path.stem) == 3
    ]
    if not manifests:
        fail(f"no three-digit vConTACT3 database manifest was found in {database_path}")
    return max(manifests, key=lambda path: int(path.stem))


def strings(value: Any) -> Iterable[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from strings(child)


def require_nonempty_file(manifest_dir: Path, relative_path: Any, label: str) -> None:
    if not isinstance(relative_path, str) or not relative_path.strip():
        fail(f"manifest entry '{label}' is missing or is not a path")
    target = manifest_dir / relative_path
    if not target.is_file() or target.stat().st_size == 0:
        fail(f"required database file for '{label}' is missing or empty: {target}")


def require_mmseqs_prefix(manifest_dir: Path, relative_prefix: Any) -> None:
    if not isinstance(relative_prefix, str) or not relative_prefix.strip():
        fail("manifest entry 'VOGDB.mmseqs_db' is missing or is not a path")
    prefix = manifest_dir / relative_prefix
    matches = [
        path
        for path in prefix.parent.glob(f"{prefix.name}*")
        if path.is_file() and path.stat().st_size > 0
    ]
    if not matches:
        fail(f"VOG MMseqs2 database prefix has no non-empty files: {prefix}")


def validate(database_path: Path) -> tuple[Path, str, tuple[str, ...], int]:
    manifest = newest_manifest(database_path.resolve())
    if manifest.stat().st_size == 0:
        fail(f"database manifest is empty: {manifest}")

    try:
        data = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"could not parse database manifest {manifest}: {exc}")

    version = manifest.stem
    release = data.get(version)
    if not isinstance(release, dict):
        fail(f"manifest does not contain its expected release key '{version}'")

    refseq = release.get("RefSeq")
    if not isinstance(refseq, dict):
        fail("manifest does not contain a RefSeq block")
    domains = refseq.get("domains")
    if not isinstance(domains, dict):
        fail("manifest RefSeq block does not contain domains")

    missing_domains = [domain for domain in REQUIRED_DOMAINS if domain not in domains]
    if missing_domains:
        fail(f"manifest is missing required domains: {', '.join(missing_domains)}")

    checked_files: set[Path] = set()
    for domain_name, domain in domains.items():
        if not isinstance(domain, dict):
            fail(f"RefSeq domain '{domain_name}' is not an object")
        for key in REQUIRED_REFSEQ_FILES:
            require_nonempty_file(manifest.parent, domain.get(key), f"{domain_name}.{key}")
            checked_files.add(manifest.parent / domain[key])

        identities = domain.get("identities")
        if not isinstance(identities, (dict, list)):
            fail(f"RefSeq domain '{domain_name}' has no identities block")
        identity_paths = tuple(strings(identities))
        if not identity_paths:
            fail(f"RefSeq domain '{domain_name}' has an empty identities block")
        for index, relative_path in enumerate(identity_paths, start=1):
            require_nonempty_file(
                manifest.parent,
                relative_path,
                f"{domain_name}.identities[{index}]",
            )
            checked_files.add(manifest.parent / relative_path)

    vogdb = release.get("VOGDB")
    if not isinstance(vogdb, dict):
        fail("manifest does not contain the VOGDB block required by database v232+")
    require_mmseqs_prefix(manifest.parent, vogdb.get("mmseqs_db"))
    for key in REQUIRED_VOG_FILES:
        require_nonempty_file(manifest.parent, vogdb.get(key), f"VOGDB.{key}")
        checked_files.add(manifest.parent / vogdb[key])

    return manifest, version, tuple(sorted(domains)), len(checked_files)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("database_path", type=Path)
    parser.add_argument(
        "--field",
        choices=("version", "manifest", "domains", "checked_files", "summary"),
        default="summary",
    )
    args = parser.parse_args()

    try:
        manifest, version, domains, checked_files = validate(args.database_path)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    values = {
        "version": version,
        "manifest": str(manifest),
        "domains": ",".join(domains),
        "checked_files": str(checked_files),
        "summary": (
            f"version={version}\tmanifest={manifest}\t"
            f"domains={','.join(domains)}\tchecked_files={checked_files}"
        ),
    }
    print(values[args.field])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
