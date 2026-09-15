#!/usr/bin/env python3
"""Build a vOTU-balanced taxonomy lookup for exact-deduplicated viCAT proteins."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

try:
    import duckdb
except ImportError as exc:  # pragma: no cover - exercised by the CLI environment
    raise SystemExit(
        "DuckDB is required. Install it with: "
        "mamba install -c conda-forge python-duckdb"
    ) from exc


RANKS = [
    ("realm", "r__Realm", "r__"),
    ("kingdom", "k__Kingdom", "k__"),
    ("phylum", "p__Phylum", "p__"),
    ("tax_class", "c__Class", "c__"),
    ("tax_order", "o__Order", "o__"),
    ("family", "f__Family", "f__"),
    ("genus", "g__Genus", "g__"),
    ("species", "s__Species", "s__"),
]
REQUIRED_METADATA_COLUMNS = {
    "uvig",
    "votu",
    "ictv_taxonomy",
    "ictv_taxonomy_method",
}
MEMORY_LIMIT_PATTERN = re.compile(r"^[1-9][0-9]*(?:\.[0-9]+)?(?:KB|MB|GB|TB)$", re.I)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Join the MMseqs2 representative/member table to MetaVR metadata, "
            "collapse repeated UViGs within each vOTU, and compute a conservative "
            "taxonomy LCA across distinct vOTUs."
        )
    )
    parser.add_argument("--cluster-members", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-limit", default="64GB")
    parser.add_argument("--prefix", default="IMGVR5_UViG")
    parser.add_argument("--expected-member-count", type=int)
    parser.add_argument("--expected-representative-count", type=int)
    parser.add_argument(
        "--allow-missing-metadata",
        action="store_true",
        help="Retain unmatched UViGs as unclassified instead of failing.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse completed intermediate tables in the work database.",
    )
    return parser.parse_args()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def clean_text(expression: str) -> str:
    return f"""
        CASE
            WHEN {expression} IS NULL THEN NULL
            WHEN trim({expression}) IN ('', '\\N') THEN NULL
            WHEN lower(trim({expression})) IN ('na', 'n/a', 'none', 'null') THEN NULL
            ELSE trim({expression})
        END
    """


def taxonomy_value(prefix: str) -> str:
    extracted = (
        "trim(regexp_extract(ictv_taxonomy, "
        f"'(^|;){prefix}([^;]*)', 2))"
    )
    return f"""
        CASE
            WHEN {extracted} = '' THEN NULL
            WHEN lower({extracted}) IN
                ('unclassified', 'unclassified virus', 'viruses', 'virus')
                THEN NULL
            ELSE {extracted}
        END
    """


def table_exists(connection: duckdb.DuckDBPyConnection, table: str) -> bool:
    return bool(
        connection.execute(
            "SELECT count(*) FROM information_schema.tables WHERE table_name = ?",
            [table],
        ).fetchone()[0]
    )


def create_stage(
    connection: duckdb.DuckDBPyConnection,
    table: str,
    query: str,
    resume: bool,
) -> None:
    if table_exists(connection, table):
        if resume:
            print(f"RESUME: reusing completed table {table}", flush=True)
            return
        raise RuntimeError(
            f"Build table already exists: {table}. Use --resume or a new work directory."
        )

    print(f"BUILD: {table}", flush=True)
    connection.execute(f"CREATE TABLE {table} AS {query}")
    connection.execute("CHECKPOINT")


def validate_args(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise SystemExit("--threads must be a positive integer")
    if not MEMORY_LIMIT_PATTERN.fullmatch(args.memory_limit):
        raise SystemExit("--memory-limit must look like 64GB, 400GB, or 1TB")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.prefix):
        raise SystemExit("--prefix may contain only letters, numbers, periods, underscores, and hyphens")
    for path, label in (
        (args.cluster_members, "cluster membership table"),
        (args.metadata, "metadata table"),
    ):
        if not path.is_file():
            raise SystemExit(f"Missing {label}: {path}")


def metadata_columns(connection: duckdb.DuckDBPyConnection, path: Path) -> set[str]:
    relation = connection.execute(
        f"""
        SELECT *
        FROM read_csv(
            {sql_string(path)},
            delim = '\t',
            header = true,
            all_varchar = true,
            sample_size = 10000
        )
        LIMIT 0
        """
    )
    return {column[0] for column in relation.description}


def build_metadata_sql(path: Path) -> str:
    rank_columns = []
    prefix_by_rank = {
        "realm": "r__",
        "kingdom": "k__",
        "phylum": "p__",
        "tax_class": "c__",
        "tax_order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }
    for rank, _, _ in RANKS:
        rank_columns.append(f"{taxonomy_value(prefix_by_rank[rank])} AS {rank}")

    uvig = clean_text("uvig")
    votu = clean_text("votu")
    method = clean_text("ictv_taxonomy_method")
    return f"""
        SELECT
            {uvig} AS uvig,
            coalesce({votu}, {uvig}) AS votu_key,
            {method} AS taxonomy_method,
            {', '.join(rank_columns)}
        FROM read_csv(
            {sql_string(path)},
            delim = '\t',
            header = true,
            all_varchar = true,
            sample_size = 10000
        )
        WHERE {uvig} IS NOT NULL
    """


def member_reader(path: Path) -> str:
    return f"""
        read_csv(
            {sql_string(path)},
            delim = '\t',
            header = false,
            columns = {{
                'representative_protein_id': 'VARCHAR',
                'member_protein_id': 'VARCHAR'
            }},
            strict_mode = true
        )
    """


def rep_votu_counts_sql(member_path: Path) -> str:
    rank_aggregates = []
    for rank, _, _ in RANKS:
        rank_aggregates.extend(
            [
                f"min({rank}) FILTER (WHERE {rank} IS NOT NULL) AS {rank}_candidate",
                f"count(DISTINCT {rank}) FILTER (WHERE {rank} IS NOT NULL) "
                f"AS {rank}_distinct_count",
            ]
        )

    return f"""
        WITH members AS (
            SELECT
                trim(representative_protein_id) AS representative_protein_id,
                trim(member_protein_id) AS member_protein_id,
                split_part(trim(member_protein_id), '|', 1) AS member_uvig
            FROM {member_reader(member_path)}
        ),
        joined AS (
            SELECT
                members.*,
                metadata_norm.uvig AS metadata_uvig,
                metadata_norm.votu_key,
                metadata_norm.taxonomy_method,
                {', '.join(rank for rank, _, _ in RANKS)}
            FROM members
            LEFT JOIN metadata_norm
                ON metadata_norm.uvig = members.member_uvig
        )
        SELECT
            representative_protein_id,
            coalesce(votu_key, '__MISSING__:' || member_uvig) AS votu_key,
            count(*)::UBIGINT AS member_protein_count,
            count(DISTINCT member_uvig)::UBIGINT AS member_uvig_count,
            count(*) FILTER (WHERE metadata_uvig IS NULL)::UBIGINT
                AS missing_metadata_member_count,
            list_distinct(
                list(taxonomy_method) FILTER (WHERE taxonomy_method IS NOT NULL)
            ) AS taxonomy_methods,
            {', '.join(rank_aggregates)}
        FROM joined
        GROUP BY representative_protein_id, coalesce(votu_key, '__MISSING__:' || member_uvig)
    """


def rep_votu_sql() -> str:
    rank_values = []
    preceding_conditions = []
    conflict_cases = []
    for index, (rank, _, _) in enumerate(RANKS, start=1):
        preceding_conditions.append(f"{rank}_distinct_count <= 1")
        rank_values.append(
            "CASE WHEN "
            + " AND ".join(preceding_conditions)
            + f" THEN {rank}_candidate ELSE NULL END AS {rank}"
        )
        conflict_cases.append(
            f"WHEN {rank}_distinct_count > 1 THEN {index}"
        )

    return f"""
        SELECT
            representative_protein_id,
            votu_key,
            member_protein_count,
            member_uvig_count,
            missing_metadata_member_count,
            taxonomy_methods,
            CASE {' '.join(conflict_cases)} ELSE 0 END AS internal_conflict_rank,
            {', '.join(rank_values)}
        FROM rep_votu_counts
    """


def representative_counts_sql() -> str:
    rank_aggregates = []
    for rank, _, _ in RANKS:
        rank_aggregates.extend(
            [
                f"min({rank}) FILTER (WHERE {rank} IS NOT NULL) AS {rank}_candidate",
                f"count(DISTINCT {rank}) FILTER (WHERE {rank} IS NOT NULL) "
                f"AS {rank}_distinct_count",
            ]
        )
    any_rank = " OR ".join(f"{rank} IS NOT NULL" for rank, _, _ in RANKS)

    return f"""
        SELECT
            representative_protein_id,
            split_part(representative_protein_id, '|', 1) AS representative_uvig,
            sum(member_protein_count)::UBIGINT AS member_protein_count,
            sum(member_uvig_count)::UBIGINT AS member_uvig_count,
            count(*)::UBIGINT AS member_votu_count,
            count(*) FILTER (WHERE {any_rank})::UBIGINT AS classified_votu_count,
            sum(missing_metadata_member_count)::UBIGINT AS missing_metadata_member_count,
            min(internal_conflict_rank) FILTER (WHERE internal_conflict_rank > 0)
                AS earliest_internal_conflict_rank,
            list_sort(list_distinct(flatten(list(taxonomy_methods)))) AS taxonomy_methods,
            {', '.join(rank_aggregates)}
        FROM rep_votu
        GROUP BY representative_protein_id
    """


def representative_final_sql() -> str:
    external_cases = []
    for index, (rank, _, _) in enumerate(RANKS, start=1):
        external_cases.append(f"WHEN {rank}_distinct_count > 1 THEN {index}")
    external_conflict = f"CASE {' '.join(external_cases)} ELSE 0 END"

    final_rank_values = []
    for index, (rank, _, _) in enumerate(RANKS, start=1):
        final_rank_values.append(
            f"CASE WHEN conflict_rank = 0 OR conflict_rank > {index} "
            f"THEN {rank}_candidate ELSE NULL END AS {rank}"
        )

    return f"""
        WITH conflict_resolution AS (
            SELECT
                *,
                CASE
                    WHEN earliest_internal_conflict_rank IS NULL
                        THEN {external_conflict}
                    WHEN {external_conflict} = 0
                        THEN earliest_internal_conflict_rank
                    ELSE least(earliest_internal_conflict_rank, {external_conflict})
                END AS conflict_rank
            FROM representative_counts
        )
        SELECT
            representative_protein_id,
            representative_uvig,
            member_protein_count,
            member_uvig_count,
            member_votu_count,
            classified_votu_count,
            missing_metadata_member_count,
            classified_votu_count::DOUBLE / member_votu_count::DOUBLE
                AS taxonomy_coverage,
            conflict_rank > 0 AS taxonomy_conflict,
            conflict_rank,
            array_to_string(taxonomy_methods, '|') AS taxonomy_methods,
            {', '.join(final_rank_values)}
        FROM conflict_resolution
    """


def final_projection(runtime: bool) -> str:
    conflict_rank_name = (
        "CASE conflict_rank "
        + " ".join(
            f"WHEN {index} THEN '{rank.replace('tax_', '')}'"
            for index, (rank, _, _) in enumerate(RANKS, start=1)
        )
        + " ELSE '' END"
    )
    prefixed_columns = []
    taxonomy_tokens = ["'d__Viruses'"]
    for rank, output_column, prefix in RANKS:
        token = f"coalesce('{prefix}' || {rank}, '{prefix}unclassified')"
        prefixed_columns.append(f"{token} AS \"{output_column}\"")
        taxonomy_tokens.append(token)

    deepest_rank = (
        "CASE "
        + " ".join(
            f"WHEN {rank} IS NOT NULL THEN '{rank.replace('tax_', '')}'"
            for rank, _, _ in reversed(RANKS)
        )
        + " ELSE 'domain' END"
    )
    common = f"""
        representative_protein_id,
        'virus' AS classification,
        concat_ws(';', {', '.join(taxonomy_tokens)}) AS taxonomy,
        {deepest_rank} AS classification_rank,
        member_votu_count,
        classified_votu_count,
        taxonomy_coverage,
        taxonomy_conflict,
        {conflict_rank_name} AS taxonomy_conflict_rank,
        coalesce(taxonomy_methods, '') AS taxonomy_methods,
        'd__Viruses' AS "d__Domain",
        {', '.join(prefixed_columns)}
    """
    if runtime:
        return common
    return f"""
        {common},
        representative_uvig,
        member_protein_count,
        member_uvig_count,
        missing_metadata_member_count
    """


def write_summary(
    connection: duckdb.DuckDBPyConnection,
    path: Path,
    metadata_count: int,
    member_count: int,
    representative_count: int,
    rep_votu_count: int,
    missing_count: int,
) -> None:
    conflict_count, unclassified_count = connection.execute(
        """
        SELECT
            count(*) FILTER (WHERE taxonomy_conflict),
            count(*) FILTER (WHERE classified_votu_count = 0)
        FROM representative_final
        """
    ).fetchone()
    rows = [
        ("metadata_uvigs", metadata_count),
        ("cluster_member_rows", member_count),
        ("representative_proteins", representative_count),
        ("representative_votu_pairs", rep_votu_count),
        ("missing_metadata_members", missing_count),
        ("representatives_with_taxonomy_conflicts", conflict_count),
        ("representatives_without_classified_votus", unclassified_count),
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    validate_args(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)

    runtime_output = args.output_dir / f"{args.prefix}.vicat_taxonomy_lookup.parquet"
    audit_output = args.output_dir / f"{args.prefix}.vicat_taxonomy_audit.parquet"
    summary_output = args.output_dir / f"{args.prefix}.vicat_taxonomy_build_summary.tsv"
    build_database = args.work_dir / f"{args.prefix}.vicat_taxonomy_build.duckdb"

    for output in (runtime_output, audit_output, summary_output):
        if output.exists():
            raise SystemExit(f"Refusing to overwrite existing output: {output}")
    if build_database.exists() and not args.resume:
        raise SystemExit(
            f"Build database already exists: {build_database}. Use --resume to reuse it."
        )

    connection = duckdb.connect(str(build_database))
    connection.execute(f"SET threads = {args.threads}")
    connection.execute(f"SET memory_limit = {sql_string(args.memory_limit)}")
    connection.execute(f"SET temp_directory = {sql_string(args.work_dir / 'spill')}")
    connection.execute("SET preserve_insertion_order = false")

    columns = metadata_columns(connection, args.metadata)
    missing_columns = sorted(REQUIRED_METADATA_COLUMNS - columns)
    if missing_columns:
        raise SystemExit(
            "Metadata table is missing required columns: " + ", ".join(missing_columns)
        )

    create_stage(
        connection,
        "metadata_norm",
        build_metadata_sql(args.metadata),
        args.resume,
    )
    metadata_count, distinct_uvigs = connection.execute(
        "SELECT count(*), count(DISTINCT uvig) FROM metadata_norm"
    ).fetchone()
    if metadata_count != distinct_uvigs:
        raise SystemExit(
            f"Metadata UViG IDs are not unique: {metadata_count} rows, "
            f"{distinct_uvigs} distinct IDs"
        )

    create_stage(
        connection,
        "rep_votu_counts",
        rep_votu_counts_sql(args.cluster_members),
        args.resume,
    )
    create_stage(connection, "rep_votu", rep_votu_sql(), args.resume)
    create_stage(
        connection,
        "representative_counts",
        representative_counts_sql(),
        args.resume,
    )
    create_stage(
        connection,
        "representative_final",
        representative_final_sql(),
        args.resume,
    )

    member_count, missing_count = connection.execute(
        """
        SELECT
            sum(member_protein_count)::UBIGINT,
            sum(missing_metadata_member_count)::UBIGINT
        FROM rep_votu_counts
        """
    ).fetchone()
    representative_count = connection.execute(
        "SELECT count(*) FROM representative_final"
    ).fetchone()[0]
    rep_votu_count = connection.execute("SELECT count(*) FROM rep_votu").fetchone()[0]

    if args.expected_member_count is not None and member_count != args.expected_member_count:
        raise SystemExit(
            f"Cluster member count mismatch: observed {member_count}, "
            f"expected {args.expected_member_count}"
        )
    if (
        args.expected_representative_count is not None
        and representative_count != args.expected_representative_count
    ):
        raise SystemExit(
            f"Representative count mismatch: observed {representative_count}, "
            f"expected {args.expected_representative_count}"
        )
    if missing_count and not args.allow_missing_metadata:
        raise SystemExit(
            f"{missing_count} cluster members could not be joined to MetaVR metadata. "
            "Investigate the IDs or rerun with --allow-missing-metadata."
        )

    print("WRITE: runtime taxonomy lookup", flush=True)
    connection.execute(
        f"""
        COPY (
            SELECT {final_projection(runtime=True)}
            FROM representative_final
            ORDER BY representative_protein_id
        )
        TO {sql_string(runtime_output)}
        (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880)
        """
    )
    print("WRITE: taxonomy audit table", flush=True)
    connection.execute(
        f"""
        COPY (
            SELECT {final_projection(runtime=False)}
            FROM representative_final
        )
        TO {sql_string(audit_output)}
        (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880)
        """
    )

    write_summary(
        connection,
        summary_output,
        metadata_count,
        member_count,
        representative_count,
        rep_votu_count,
        missing_count,
    )
    connection.close()

    print(f"Runtime lookup: {runtime_output}")
    print(f"Audit table: {audit_output}")
    print(f"Build summary: {summary_output}")
    print(f"Resumable build database: {build_database}")


if __name__ == "__main__":
    main()
