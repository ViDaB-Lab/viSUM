#!/usr/bin/env python3
"""Join one viCAT alignment table to a compact run-wide reference subset."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import duckdb

from evidence_schema import TAXONOMY_COLUMNS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--diamond", type=Path)
    parser.add_argument("--reference-subset", type=Path)
    parser.add_argument("--reference-subset-metadata", type=Path)
    parser.add_argument("--viral-diamond", type=Path)
    parser.add_argument("--viral-reference-subset", type=Path)
    parser.add_argument("--viral-reference-subset-metadata", type=Path)
    parser.add_argument("--nonviral-diamond", type=Path)
    parser.add_argument("--nonviral-metadata", type=Path)
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


def validate_input_mode(args: argparse.Namespace) -> bool:
    legacy = (args.diamond, args.reference_subset, args.reference_subset_metadata)
    dual = (
        args.viral_diamond,
        args.viral_reference_subset,
        args.viral_reference_subset_metadata,
        args.nonviral_diamond,
        args.nonviral_metadata,
    )
    if any(legacy) and any(dual):
        raise SystemExit("Legacy and dual-database viCAT inputs cannot be combined")
    if all(legacy):
        return False
    if all(dual):
        return True
    if any(legacy) or any(dual):
        raise SystemExit("Incomplete viCAT hit-preparation input set")
    raise SystemExit("No viCAT alignment inputs were provided")


def main() -> None:
    args = parse_args()
    dual_mode = validate_input_mode(args)
    competitive_mode = (
        True
        if dual_mode
        else read_competitive_mode(args.reference_subset_metadata)
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    connection = configure_connection(args)
    try:
        if dual_mode:
            taxonomy_sql = ", ".join(
                f'r."{column}" AS "{column}"' for column in TAXONOMY_COLUMNS
            )
            blank_taxonomy_sql = ", ".join(
                f"'' AS \"{column}\"" for column in TAXONOMY_COLUMNS
            )
            connection.execute(
                f"""
                COPY (
                    WITH viral_hits AS (
                        SELECT
                            d.qseqid, d.sseqid, d.pident, d.length, d.qlen,
                            d.slen, d.qstart, d.qend, d.sstart, d.send,
                            d.evalue, d.bitscore, d.qcovhsp, d.scovhsp,
                            'viral' AS reference_class,
                            '' AS nonviral_reference_class,
                            r.source_protein_id, r.cellular_group,
                            r.source_accession, r.organism_name, r.taxid,
                            '' AS source_classes, '' AS replicon_accessions,
                            '' AS replicon_types,
                            false AS provirus_flank_eligible,
                            '' AS classification_note,
                            0::BIGINT AS cluster_member_count,
                            r.member_votu_count, r.classified_votu_count,
                            r.taxonomy_coverage,
                            r.reference_taxonomy_conflict,
                            r.reference_taxonomy_conflict_rank,
                            r.taxonomy_methods,
                            {taxonomy_sql},
                            r.reference_id IS NOT NULL AS manifest_matched,
                            r.reference_id IS NOT NULL AS taxonomy_matched,
                            true AS competitive_mode
                        FROM read_csv(
                            {sql_literal(args.viral_diamond)}, delim='\t',
                            header=true, auto_detect=true
                        ) d
                        LEFT JOIN read_parquet(
                            {sql_literal(args.viral_reference_subset)}
                        ) r ON d.sseqid = r.reference_id
                    ),
                    nonviral_hits AS (
                        SELECT
                            d.qseqid, d.sseqid, d.pident, d.length, d.qlen,
                            d.slen, d.qstart, d.qend, d.sstart, d.send,
                            d.evalue, d.bitscore, d.qcovhsp, d.scovhsp,
                            'cellular' AS reference_class,
                            n.reference_class AS nonviral_reference_class,
                            n.source_protein_id, n.cellular_group,
                            n.source_accession, n.organism_name, n.taxid,
                            n.source_classes, n.replicon_accessions,
                            n.replicon_types, n.provirus_flank_eligible,
                            n.classification_note, n.cluster_member_count,
                            0::BIGINT AS member_votu_count,
                            0::BIGINT AS classified_votu_count,
                            0.0::DOUBLE AS taxonomy_coverage,
                            false AS reference_taxonomy_conflict,
                            '' AS reference_taxonomy_conflict_rank,
                            '' AS taxonomy_methods,
                            {blank_taxonomy_sql},
                            n.reference_id IS NOT NULL AS manifest_matched,
                            true AS taxonomy_matched,
                            true AS competitive_mode
                        FROM read_csv(
                            {sql_literal(args.nonviral_diamond)}, delim='\t',
                            header=true, auto_detect=true
                        ) d
                        LEFT JOIN read_parquet(
                            {sql_literal(args.nonviral_metadata)}
                        ) n ON d.sseqid = n.reference_id
                    )
                    SELECT * FROM viral_hits
                    UNION ALL
                    SELECT * FROM nonviral_hits
                ) TO {sql_literal(args.output)}
                (FORMAT PARQUET, COMPRESSION ZSTD)
                """
            )
        else:
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
                    LEFT JOIN read_parquet(
                        {sql_literal(args.reference_subset)}
                    ) r ON d.sseqid = r.reference_id
                ) TO {sql_literal(args.output)}
                (FORMAT PARQUET, COMPRESSION ZSTD)
                """
            )
        (
            missing_references,
            alignment_count,
            hit_orf_count,
            competitive_modes,
            viral_alignment_count,
            nonviral_alignment_count,
        ) = connection.execute(
            f"""
            SELECT
                count(*) FILTER (WHERE NOT manifest_matched),
                count(*),
                count(DISTINCT qseqid),
                count(DISTINCT competitive_mode),
                count(*) FILTER (WHERE reference_class = 'viral'),
                count(*) FILTER (WHERE reference_class != 'viral')
            FROM read_parquet({sql_literal(args.output)})
            """
        ).fetchone()
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
            "alignment_count\tviral_alignment_count\tnonviral_alignment_count\t"
            "hit_orf_count\tcompetitive_mode\tdatabase_architecture\tthreads\t"
            "memory_limit\tstorage_format\n"
        )
        handle.write(
            f"{alignment_count}\t{viral_alignment_count}\t"
            f"{nonviral_alignment_count}\t{hit_orf_count}\t"
            f"{'true' if competitive_mode else 'false'}\t"
            f"{'separate_viral_nonviral' if dual_mode else 'legacy'}\t"
            f"{args.threads}\t{args.memory_limit}\tparquet_zstd\n"
        )
    print(
        f"viCAT hits prepared: alignments={alignment_count} "
        f"hit_orfs={hit_orf_count} output={args.output}"
    )


if __name__ == "__main__":
    main()
