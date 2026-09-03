#!/usr/bin/env python3
"""Join one viCAT alignment table to a compact run-wide reference subset."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import duckdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diamond", required=True, type=Path)
    parser.add_argument("--reference-subset", required=True, type=Path)
    parser.add_argument("--reference-subset-metadata", required=True, type=Path)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--memory-limit", required=True)
    parser.add_argument("--temp-directory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-metadata", required=True, type=Path)
    return parser.parse_args()


def sql_literal(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def configure_connection(args: argparse.Namespace) -> duckdb.DuckDBPyConnection:
    if args.threads < 1:
        raise SystemExit("--threads must be a positive integer")
    args.temp_directory.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect()
    connection.execute(f"SET threads = {args.threads}")
    connection.execute(f"SET memory_limit = {sql_literal(args.memory_limit)}")
    connection.execute(
        f"SET temp_directory = {sql_literal(args.temp_directory.resolve())}"
    )
    connection.execute("SET preserve_insertion_order = false")
    return connection


def read_competitive_mode(path: Path) -> bool:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise SystemExit(
            "viCAT reference-subset metadata must contain exactly one data row"
        )
    value = rows[0].get("competitive_mode", "").strip().lower()
    if value not in {"true", "false"}:
        raise SystemExit(
            "viCAT reference-subset metadata has an invalid competitive_mode"
        )
    return value == "true"


def main() -> None:
    args = parse_args()
    competitive_mode = read_competitive_mode(args.reference_subset_metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    connection = configure_connection(args)
    try:
        connection.execute(
            f"""
            COPY (
                SELECT
                    d.qseqid, d.sseqid, d.pident, d.length, d.qlen, d.slen,
                    d.qstart, d.qend, d.sstart, d.send, d.evalue, d.bitscore,
                    d.qcovhsp, d.scovhsp,
                    r.* EXCLUDE (reference_id),
                    r.reference_id IS NOT NULL AS manifest_matched,
                    r.reference_id IS NOT NULL AS taxonomy_matched
                FROM read_csv(
                    {sql_literal(args.diamond)}, delim='\t', header=true,
                    auto_detect=true
                ) d
                LEFT JOIN read_parquet({sql_literal(args.reference_subset)}) r
                  ON d.sseqid = r.reference_id
            ) TO {sql_literal(args.output)} (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
        missing_references, alignment_count, hit_orf_count, competitive_modes = (
            connection.execute(
                f"""
                SELECT
                    count(*) FILTER (WHERE NOT manifest_matched),
                    count(*),
                    count(DISTINCT qseqid),
                    count(DISTINCT competitive_mode)
                FROM read_parquet({sql_literal(args.output)})
                """
            ).fetchone()
        )
        if missing_references:
            args.output.unlink(missing_ok=True)
            raise SystemExit(
                f"{missing_references} DIAMOND alignment(s) absent from the "
                "run-wide viCAT reference subset"
            )
        if competitive_modes > 1:
            args.output.unlink(missing_ok=True)
            raise SystemExit("Mixed competitive modes in viCAT reference subset")
    finally:
        connection.close()

    with args.output_metadata.open("w", encoding="utf-8", newline="") as handle:
        handle.write(
            "alignment_count\thit_orf_count\tcompetitive_mode\tthreads\t"
            "memory_limit\tstorage_format\n"
        )
        handle.write(
            f"{alignment_count}\t{hit_orf_count}\t"
            f"{'true' if competitive_mode else 'false'}\t"
            f"{args.threads}\t{args.memory_limit}\tparquet_zstd\n"
        )
    print(
        f"viCAT hits prepared: alignments={alignment_count} "
        f"hit_orfs={hit_orf_count} output={args.output}"
    )


if __name__ == "__main__":
    main()
