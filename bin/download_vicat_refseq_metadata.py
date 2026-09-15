#!/usr/bin/env python3
"""Download RefSeq feature tables and assembly reports for viCAT references."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


REQUIRED_MANIFEST_COLUMNS = {"assembly_accession", "cellular_group", "ftp_path"}
REQUIRED_FEATURE_COLUMNS = {
    "feature",
    "seq_type",
    "genomic_accession",
    "product_accession",
}
SAFE_COMPONENT = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class Assembly:
    accession: str
    group: str
    ftp_path: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=120)
    return parser.parse_args()


def read_manifest(path: Path) -> list[Assembly]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest is missing columns: {sorted(missing)}")
        assemblies = []
        seen = set()
        for row in reader:
            accession = row["assembly_accession"].strip()
            group = row["cellular_group"].strip()
            ftp_path = row["ftp_path"].strip().rstrip("/")
            if not accession or not group or not ftp_path:
                raise ValueError("Assembly accession, cellular group, and FTP path cannot be blank")
            if not SAFE_COMPONENT.fullmatch(accession) or not SAFE_COMPONENT.fullmatch(group):
                raise ValueError(f"Unsafe assembly or group value: {accession!r}, {group!r}")
            if accession in seen:
                raise ValueError(f"Duplicate assembly accession: {accession}")
            seen.add(accession)
            assemblies.append(Assembly(accession, group, ftp_path))
    if not assemblies:
        raise ValueError("Manifest contains no assemblies")
    return assemblies


def source_urls(assembly: Assembly) -> tuple[str, str]:
    basename = assembly.ftp_path.rsplit("/", 1)[-1]
    return (
        f"{assembly.ftp_path}/{basename}_feature_table.txt.gz",
        f"{assembly.ftp_path}/{basename}_assembly_report.txt",
    )


def validate_feature_table(path: Path) -> None:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = [field.lstrip("# ").strip() for field in next(csv.reader(handle, delimiter="\t"))]
    missing = REQUIRED_FEATURE_COLUMNS - set(header)
    if missing:
        raise ValueError(f"Feature table is missing columns: {sorted(missing)}")


def validate_report(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Assembly report is empty")


def valid_existing(path: Path, validator) -> bool:
    try:
        validator(path)
    except (OSError, EOFError, StopIteration, UnicodeError, ValueError):
        return False
    return True


def download(url: str, destination: Path, validator, retries: int, timeout: int) -> str:
    if valid_existing(destination, validator):
        return "cached"
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    partial.unlink(missing_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "viSUM-viCAT-builder/1"})
            with urllib.request.urlopen(request, timeout=timeout) as response, partial.open("wb") as output:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
            validator(partial)
            os.replace(partial, destination)
            return "downloaded"
        except (OSError, EOFError, StopIteration, UnicodeError, ValueError, urllib.error.URLError) as error:
            last_error = error
            partial.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(2 ** (attempt - 1), 10))
    raise RuntimeError(f"Failed after {retries} attempt(s): {url}: {last_error}")


def download_assembly(assembly: Assembly, output_dir: Path, retries: int, timeout: int) -> dict[str, str]:
    target = output_dir / assembly.group / assembly.accession
    feature = target / f"{assembly.accession}.feature_table.txt.gz"
    report = target / f"{assembly.accession}.assembly_report.txt"
    feature_url, report_url = source_urls(assembly)
    feature_status = download(feature_url, feature, validate_feature_table, retries, timeout)
    report_status = download(report_url, report, validate_report, retries, timeout)
    return {
        "cellular_group": assembly.group,
        "assembly_accession": assembly.accession,
        "feature_table": str(feature),
        "assembly_report": str(report),
        "feature_status": feature_status,
        "report_status": report_status,
        "error": "",
    }


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "cellular_group", "assembly_accession", "feature_table", "assembly_report",
        "feature_status", "report_status", "error",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["cellular_group"], row["assembly_accession"])))


def main() -> None:
    args = parse_args()
    if args.workers < 1 or args.retries < 1 or args.timeout < 1:
        raise SystemExit("--workers, --retries, and --timeout must be positive")
    try:
        assemblies = read_manifest(args.manifest)
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(download_assembly, assembly, args.output_dir, args.retries, args.timeout): assembly
            for assembly in assemblies
        }
        for future in as_completed(futures):
            assembly = futures[future]
            try:
                row = future.result()
                print(f"METADATA_OK group={assembly.group} assembly={assembly.accession} "
                      f"feature={row['feature_status']} report={row['report_status']}")
            except Exception as error:  # preserve all failures in the audit before exiting
                row = {
                    "cellular_group": assembly.group,
                    "assembly_accession": assembly.accession,
                    "feature_table": "",
                    "assembly_report": "",
                    "feature_status": "failed",
                    "report_status": "failed",
                    "error": str(error),
                }
                print(f"METADATA_FAILED group={assembly.group} assembly={assembly.accession}: {error}")
            rows.append(row)

    summary = args.output_dir / "download_summary.tsv"
    write_summary(summary, rows)
    failures = sum(bool(row["error"]) for row in rows)
    print(f"Assemblies: {len(rows)}; failures: {failures}; summary: {summary}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
