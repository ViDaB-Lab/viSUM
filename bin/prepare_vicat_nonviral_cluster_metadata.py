#!/usr/bin/env python3
"""Validate viCAT nonviral clusters and package representative metadata."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import duckdb


REQUIRED_METADATA_COLUMNS = {
    "reference_id",
    "source_protein_id",
    "reference_class",
    "source_classes",
    "cellular_group",
    "source_accession",
    "replicon_accessions",
    "replicon_types",
    "provirus_flank_eligible",
    "classification_note",
    "organism_name",
    "taxid",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--output-parquet", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--work-database", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    return parser.parse_args()


def sql_path(path: Path) -> str:
    return "'" + str(path.resolve()).replace("'", "''") + "'"


def table_columns(connection: duckdb.DuckDBPyConnection, expression: str) -> set[str]:
    return {
        row[0]
        for row in connection.execute(f"DESCRIBE SELECT * FROM {expression}").fetchall()
    }


def scalar(connection: duckdb.DuckDBPyConnection, query: str) -> int:
    return int(connection.execute(query).fetchone()[0])


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise SystemExit("ERROR: --threads must be a positive integer")
    for source in (args.metadata, args.membership):
        if not source.is_file():
            raise SystemExit(f"ERROR: Input not found: {source}")

    for destination in (
        args.output_tsv,
        args.output_parquet,
        args.output_summary,
        args.work_database,
    ):
        destination.parent.mkdir(parents=True, exist_ok=True)

    connection = duckdb.connect(str(args.work_database))
    try:
        connection.execute(f"SET threads = {args.threads}")
        metadata_source = (
            f"read_csv({sql_path(args.metadata)}, delim='\\t', header=true, "
            "all_varchar=true)"
        )
        membership_source = (
            f"read_csv({sql_path(args.membership)}, delim='\\t', header=true, "
            "all_varchar=true)"
        )

        missing_columns = REQUIRED_METADATA_COLUMNS - table_columns(
            connection, metadata_source
        )
        if missing_columns:
            raise ValueError(
                f"Nonviral metadata is missing columns: {sorted(missing_columns)}"
            )
        membership_columns = table_columns(connection, membership_source)
        expected_membership_columns = {
            "reference_class",
            "representative_id",
            "member_id",
        }
        if membership_columns != expected_membership_columns:
            raise ValueError(
                "Cluster membership columns differ from the required schema: "
                f"observed {sorted(membership_columns)}"
            )

        connection.execute(
            f"CREATE OR REPLACE TABLE source_metadata AS SELECT * FROM {metadata_source}"
        )
        connection.execute(
            "CREATE OR REPLACE TABLE cluster_membership AS "
            f"SELECT * FROM {membership_source}"
        )
        metadata = "source_metadata"
        membership = "cluster_membership"

        duplicate_metadata = scalar(
            connection,
            f"""
            SELECT count(*) FROM (
                SELECT reference_id FROM {metadata}
                GROUP BY reference_id HAVING count(*) != 1
            )
            """,
        )
        if duplicate_metadata:
            raise ValueError(
                f"Nonviral metadata has {duplicate_metadata} duplicated reference IDs"
            )

        invalid_members = scalar(
            connection,
            f"""
            SELECT count(*) FROM (
                SELECT member_id FROM {membership}
                GROUP BY member_id HAVING count(*) != 1
            )
            """,
        )
        if invalid_members:
            raise ValueError(
                f"Cluster membership has {invalid_members} non-unique member IDs"
            )

        missing_metadata = scalar(
            connection,
            f"""
            SELECT count(*) FROM {membership} c
            LEFT JOIN {metadata} m ON c.member_id = m.reference_id
            WHERE m.reference_id IS NULL
            """,
        )
        unused_metadata = scalar(
            connection,
            f"""
            SELECT count(*) FROM {metadata} m
            LEFT JOIN {membership} c ON c.member_id = m.reference_id
            WHERE c.member_id IS NULL
            """,
        )
        class_mismatches = scalar(
            connection,
            f"""
            SELECT count(*) FROM {membership} c
            JOIN {metadata} m ON c.member_id = m.reference_id
            WHERE c.reference_class != m.reference_class
            """,
        )
        missing_representatives = scalar(
            connection,
            f"""
            SELECT count(*) FROM (
                SELECT DISTINCT representative_id FROM {membership}
            ) r
            LEFT JOIN {metadata} m ON r.representative_id = m.reference_id
            WHERE m.reference_id IS NULL
            """,
        )
        if missing_metadata or unused_metadata or class_mismatches or missing_representatives:
            raise ValueError(
                "Cluster/metadata validation failed: "
                f"missing_metadata={missing_metadata}, "
                f"unused_metadata={unused_metadata}, "
                f"class_mismatches={class_mismatches}, "
                f"missing_representatives={missing_representatives}"
            )

        connection.execute(
            f"""
            CREATE OR REPLACE TEMP VIEW representative_metadata AS
            WITH cluster_counts AS (
                SELECT
                    reference_class,
                    representative_id,
                    count(*)::BIGINT AS cluster_member_count
                FROM {membership}
                GROUP BY reference_class, representative_id
            )
            SELECT
                m.reference_id,
                m.source_protein_id,
                m.reference_class,
                c.cluster_member_count,
                m.source_classes,
                m.cellular_group,
                m.source_accession,
                m.replicon_accessions,
                m.replicon_types,
                CAST(m.provirus_flank_eligible AS BOOLEAN) AS provirus_flank_eligible,
                m.classification_note,
                m.organism_name,
                m.taxid
            FROM cluster_counts c
            JOIN {metadata} m ON c.representative_id = m.reference_id
            ORDER BY m.reference_class, m.reference_id
            """
        )

        for destination in (args.output_tsv, args.output_parquet):
            if destination.exists():
                destination.unlink()
        connection.execute(
            f"""
            COPY representative_metadata TO {sql_path(args.output_tsv)}
            (FORMAT CSV, HEADER true, DELIMITER '\\t', COMPRESSION gzip)
            """
        )
        connection.execute(
            f"""
            COPY representative_metadata TO {sql_path(args.output_parquet)}
            (FORMAT PARQUET, COMPRESSION zstd)
            """
        )

        summary_rows = connection.execute(
            f"""
            WITH inputs AS (
                SELECT reference_class, count(*)::BIGINT AS input_proteins
                FROM {membership} GROUP BY reference_class
            ), representatives AS (
                SELECT reference_class, count(*)::BIGINT AS representatives
                FROM representative_metadata GROUP BY reference_class
            )
            SELECT
                i.reference_class,
                i.input_proteins,
                r.representatives,
                i.input_proteins - r.representatives AS proteins_removed,
                round(r.representatives * 100.0 / i.input_proteins, 4) AS retained_percent
            FROM inputs i JOIN representatives r USING (reference_class)
            ORDER BY i.reference_class
            """
        ).fetchall()

        with args.output_summary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "reference_class",
                    "input_proteins",
                    "representatives",
                    "proteins_removed",
                    "retained_percent",
                ]
            )
            writer.writerows(summary_rows)

        total_inputs = sum(int(row[1]) for row in summary_rows)
        total_representatives = sum(int(row[2]) for row in summary_rows)
        print(f"Cluster members: {total_inputs}")
        print(f"Representatives: {total_representatives}")
        print(f"Summary: {args.output_summary}")
    except (duckdb.Error, OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    finally:
        connection.close()


if __name__ == "__main__":
    main()
