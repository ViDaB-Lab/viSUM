#!/usr/bin/env python3
"""Join viCAT alignments to reference metadata once and cache them as Parquet."""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb

from evidence_schema import TAXONOMY_COLUMNS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diamond", required=True, type=Path)
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


def joined_query(args: argparse.Namespace) -> str:
    diamond = sql_literal(args.diamond)
    lookup = sql_literal(args.taxonomy_lookup)
    rank_sql = ", ".join(
        f't."{column}" AS "{column}"' for column in TAXONOMY_COLUMNS
    )
    alignment_columns = """
        d.qseqid, d.sseqid, d.pident, d.length, d.qlen, d.slen,
        d.qstart, d.qend, d.sstart, d.send, d.evalue, d.bitscore,
        d.qcovhsp, d.scovhsp
    """
    if args.reference_manifest is not None:
        manifest = sql_literal(args.reference_manifest)
        return f"""
            SELECT
                {alignment_columns},
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
            FROM read_csv(
                {diamond}, delim='\t', header=true, auto_detect=true
            ) d
            LEFT JOIN read_parquet({manifest}) m
              ON d.sseqid = m.reference_id
            LEFT JOIN read_parquet({lookup}) t
              ON m.reference_class = 'viral'
             AND m.source_protein_id = t.representative_protein_id
        """
    return f"""
        SELECT
            {alignment_columns},
            'viral' AS reference_class, d.sseqid AS source_protein_id,
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
        FROM read_csv(
            {diamond}, delim='\t', header=true, auto_detect=true
        ) d
        LEFT JOIN read_parquet({lookup}) t
          ON d.sseqid = t.representative_protein_id
    """


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    connection = configure_connection(args)
    try:
        connection.execute(
            f"COPY ({joined_query(args)}) TO {sql_literal(args.output)} "
            "(FORMAT PARQUET, COMPRESSION ZSTD)"
        )
        missing_manifest, missing_taxonomy, alignment_count, hit_orf_count = (
            connection.execute(
                f"""
                SELECT
                    count(*) FILTER (WHERE NOT manifest_matched),
                    count(*) FILTER (WHERE NOT taxonomy_matched),
                    count(*),
                    count(DISTINCT qseqid)
                FROM read_parquet({sql_literal(args.output)})
                """
            ).fetchone()
        )
        if missing_manifest:
            args.output.unlink(missing_ok=True)
            raise SystemExit(
                f"{missing_manifest} DIAMOND alignment(s) absent from the "
                "competitive reference manifest"
            )
        if missing_taxonomy:
            args.output.unlink(missing_ok=True)
            raise SystemExit(
                f"{missing_taxonomy} viral DIAMOND alignment(s) reference "
                "proteins absent from the viCAT taxonomy lookup"
            )
    finally:
        connection.close()

    with args.output_metadata.open("w", encoding="utf-8", newline="") as handle:
        handle.write(
            "alignment_count\thit_orf_count\tcompetitive_mode\tthreads\t"
            "memory_limit\tstorage_format\n"
        )
        handle.write(
            f"{alignment_count}\t{hit_orf_count}\t"
            f"{'true' if args.reference_manifest is not None else 'false'}\t"
            f"{args.threads}\t{args.memory_limit}\tparquet_zstd\n"
        )
    print(
        f"viCAT hits prepared: alignments={alignment_count} "
        f"hit_orfs={hit_orf_count} output={args.output}"
    )


if __name__ == "__main__":
    main()
