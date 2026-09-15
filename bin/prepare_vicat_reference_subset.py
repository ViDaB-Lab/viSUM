#!/usr/bin/env python3
"""Subset viCAT reference metadata once for all alignments in a pipeline run."""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb

from evidence_schema import TAXONOMY_COLUMNS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diamond", required=True, type=Path, nargs="+")
    parser.add_argument("--taxonomy-lookup", required=True, type=Path)
    parser.add_argument("--reference-manifest", type=Path)
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


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    diamond_paths = ", ".join(sql_literal(path) for path in args.diamond)
    lookup = sql_literal(args.taxonomy_lookup)
    rank_sql = ", ".join(
        f't."{column}" AS "{column}"' for column in TAXONOMY_COLUMNS
    )
    connection = configure_connection(args)
    try:
        connection.execute(
            f"""
            CREATE TEMP TABLE relevant_reference_ids AS
            SELECT DISTINCT sseqid AS reference_id
            FROM read_csv(
                [{diamond_paths}], delim='\t', header=true,
                auto_detect=true, union_by_name=true
            )
            """
        )
        if args.reference_manifest is not None:
            manifest = sql_literal(args.reference_manifest)
            connection.execute(
                f"""
                CREATE TEMP TABLE selected_references AS
                SELECT
                    r.reference_id,
                    m.reference_class, m.source_protein_id, m.cellular_group,
                    m.source_accession, m.organism_name, m.taxid,
                    coalesce(t.member_votu_count, 0) AS member_votu_count,
                    coalesce(t.classified_votu_count, 0) AS classified_votu_count,
                    coalesce(t.taxonomy_coverage, 0) AS taxonomy_coverage,
                    coalesce(t.taxonomy_conflict, false) AS reference_taxonomy_conflict,
                    coalesce(t.taxonomy_conflict_rank, '') AS reference_taxonomy_conflict_rank,
                    coalesce(t.taxonomy_methods, '') AS taxonomy_methods,
                    {rank_sql},
                    m.reference_id IS NOT NULL AS manifest_matched,
                    m.reference_class != 'viral'
                        OR t.representative_protein_id IS NOT NULL AS taxonomy_matched,
                    true AS competitive_mode
                FROM relevant_reference_ids r
                LEFT JOIN read_parquet({manifest}) m
                  ON r.reference_id = m.reference_id
                LEFT JOIN read_parquet({lookup}) t
                  ON m.reference_class = 'viral'
                 AND m.source_protein_id = t.representative_protein_id
                """
            )
        else:
            connection.execute(
                f"""
                CREATE TEMP TABLE selected_references AS
                SELECT
                    r.reference_id,
                    'viral' AS reference_class,
                    r.reference_id AS source_protein_id,
                    '' AS cellular_group, '' AS source_accession,
                    '' AS organism_name, '' AS taxid,
                    coalesce(t.member_votu_count, 0) AS member_votu_count,
                    coalesce(t.classified_votu_count, 0) AS classified_votu_count,
                    coalesce(t.taxonomy_coverage, 0) AS taxonomy_coverage,
                    coalesce(t.taxonomy_conflict, false) AS reference_taxonomy_conflict,
                    coalesce(t.taxonomy_conflict_rank, '') AS reference_taxonomy_conflict_rank,
                    coalesce(t.taxonomy_methods, '') AS taxonomy_methods,
                    {rank_sql},
                    true AS manifest_matched,
                    t.representative_protein_id IS NOT NULL AS taxonomy_matched,
                    false AS competitive_mode
                FROM relevant_reference_ids r
                LEFT JOIN read_parquet({lookup}) t
                  ON r.reference_id = t.representative_protein_id
                """
            )

        reference_count, missing_manifest, missing_taxonomy = connection.execute(
            """
            SELECT
                count(*),
                count(*) FILTER (WHERE NOT manifest_matched),
                count(*) FILTER (WHERE NOT taxonomy_matched)
            FROM selected_references
            """
        ).fetchone()
        if missing_manifest:
            raise SystemExit(
                f"{missing_manifest} DIAMOND reference ID(s) absent from the "
                "competitive reference manifest"
            )
        if missing_taxonomy:
            raise SystemExit(
                f"{missing_taxonomy} viral reference protein(s) absent from "
                "the viCAT taxonomy lookup"
            )
        connection.execute(
            f"""
            COPY (
                SELECT * EXCLUDE (manifest_matched, taxonomy_matched)
                FROM selected_references
            ) TO {sql_literal(args.output)} (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
    finally:
        connection.close()

    with args.output_metadata.open("w", encoding="utf-8", newline="") as handle:
        handle.write(
            "sample_count\treference_count\tcompetitive_mode\tthreads\t"
            "memory_limit\tstorage_format\n"
        )
        handle.write(
            f"{len(args.diamond)}\t{reference_count}\t"
            f"{'true' if args.reference_manifest is not None else 'false'}\t"
            f"{args.threads}\t{args.memory_limit}\tparquet_zstd\n"
        )
    print(
        f"viCAT reference subset prepared: samples={len(args.diamond)} "
        f"references={reference_count} output={args.output}"
    )


if __name__ == "__main__":
    main()
