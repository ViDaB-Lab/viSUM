#!/usr/bin/env python3
"""Label viral/cellular viCAT references and build an auditable manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import TextIO

import duckdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--viral-representatives", required=True, type=Path)
    parser.add_argument("--viral-taxonomy-lookup", required=True, type=Path)
    parser.add_argument("--cellular-representatives", required=True, type=Path)
    parser.add_argument("--cellular-metadata", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def sql_path(path: Path) -> str:
    return "'" + str(path.resolve()).replace("'", "''") + "'"


def normalize_cellular_metadata(source: Path, destination: Path) -> set[str]:
    """Write a canonical TSV, accepting real tabs or shell-safe literal ``\\t``."""
    required = {"protein_id", "cellular_group", "source_accession"}
    row_count = 0

    with source.open(encoding="utf-8-sig", newline="") as source_handle:
        first_line = source_handle.readline()
        if not first_line:
            raise ValueError(f"Cellular metadata is empty: {source}")

        if "\t" in first_line:
            replace_literal_tabs = False
        elif "\\t" in first_line:
            replace_literal_tabs = True
            print(
                "WARNING: Cellular metadata uses literal \\t separators; "
                "normalizing them to tab characters."
            )
        else:
            raise ValueError(
                "Cellular metadata is not tab-delimited. "
                f"Observed header: {first_line.rstrip()}"
            )

        def normalized_lines():
            line = first_line
            while line:
                yield line.replace("\\t", "\t") if replace_literal_tabs else line
                line = source_handle.readline()

        reader = csv.DictReader(normalized_lines(), delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Cellular metadata has no header: {source}")
        fieldnames = [name.strip() for name in reader.fieldnames]
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError(
                f"Cellular metadata has duplicate columns: {fieldnames}"
            )
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(
                f"Cellular metadata is missing columns: {sorted(missing)}. "
                f"Observed columns: {fieldnames}"
            )

        with destination.open("w", encoding="utf-8", newline="") as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row_number, row in enumerate(reader, start=2):
                if None in row:
                    raise ValueError(
                        "Cellular metadata row has more values than its header "
                        f"at line {row_number}"
                    )
                normalized_row = {
                    key.strip(): (value if value is not None else "")
                    for key, value in row.items()
                }
                writer.writerow(normalized_row)
                row_count += 1

    if row_count == 0:
        raise ValueError(f"Cellular metadata has a header but no data rows: {source}")
    return set(fieldnames)


def label_fasta(source: Path, output: TextIO, label: str, ids: TextIO | None) -> int:
    count = 0
    sequence_seen = False
    with open_text(source) as handle:
        for raw_line in handle:
            if raw_line.startswith(">"):
                if count and not sequence_seen:
                    raise ValueError(
                        f"FASTA record before header {raw_line.strip()} has no sequence in {source}"
                    )
                source_id = raw_line[1:].strip().split(maxsplit=1)[0]
                if not source_id:
                    raise ValueError(f"Empty FASTA identifier in {source}")
                count += 1
                sequence_seen = False
                output.write(f">{label}|{source_id}\n")
                if ids is not None:
                    ids.write(source_id + "\n")
            else:
                sequence = raw_line.strip()
                if sequence:
                    if count == 0:
                        raise ValueError(f"Sequence data precedes the first FASTA header in {source}")
                    sequence_seen = True
                    output.write(sequence + "\n")
        if count and not sequence_seen:
            raise ValueError(f"Final FASTA record has no sequence in {source}")
    if count == 0:
        raise ValueError(f"No FASTA records found in {source}")
    return count


def table_columns(connection: duckdb.DuckDBPyConnection, expression: str) -> set[str]:
    rows = connection.execute(f"DESCRIBE SELECT * FROM {expression}").fetchall()
    return {row[0] for row in rows}


def quoted_or_default(columns: set[str], name: str, default: str = "''") -> str:
    return f'CAST("{name}" AS VARCHAR)' if name in columns else default


def main() -> None:
    args = parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    cellular_ids = args.work_dir / "cellular_representative_ids.tsv"
    normalized_cellular_metadata = args.work_dir / "cellular_metadata.normalized.tsv"
    normalized_cellular_columns = normalize_cellular_metadata(
        args.cellular_metadata, normalized_cellular_metadata
    )

    with args.output_fasta.open("w", encoding="utf-8", newline="\n") as output:
        viral_count = label_fasta(args.viral_representatives, output, "VIRAL", None)
        with cellular_ids.open("w", encoding="utf-8", newline="\n") as ids:
            ids.write("protein_id\n")
            cellular_count = label_fasta(
                args.cellular_representatives, output, "CELLULAR", ids
            )

    connection = duckdb.connect(str(args.work_dir / "reference_manifest.duckdb"))
    viral_expression = f"read_parquet({sql_path(args.viral_taxonomy_lookup)})"
    cellular_expression = (
        f"read_csv({sql_path(normalized_cellular_metadata)}, delim='\\t', header=true, "
        "auto_detect=true, all_varchar=true)"
    )
    viral_columns = table_columns(connection, viral_expression)
    cellular_columns = table_columns(connection, cellular_expression)

    required_viral = {
        "representative_protein_id", "taxonomy", "classification_rank",
        "taxonomy_conflict", "taxonomy_conflict_rank",
    }
    missing_viral = required_viral - viral_columns
    if missing_viral:
        raise ValueError(
            f"Viral taxonomy lookup is missing columns: {sorted(missing_viral)}"
        )
    if cellular_columns != normalized_cellular_columns:
        raise ValueError(
            "Cellular metadata columns changed while loading the normalized TSV: "
            f"expected {sorted(normalized_cellular_columns)}, "
            f"observed {sorted(cellular_columns)}"
        )

    duplicate_metadata = connection.execute(
        f"""
        SELECT count(*) FROM (
            SELECT protein_id
            FROM {cellular_expression}
            GROUP BY protein_id
            HAVING count(*) > 1
        )
        """
    ).fetchone()[0]
    if duplicate_metadata:
        raise ValueError(
            f"Cellular metadata contains {duplicate_metadata} duplicate protein IDs"
        )

    duplicate_viral = connection.execute(
        f"""
        SELECT count(*) FROM (
            SELECT representative_protein_id
            FROM {viral_expression}
            GROUP BY representative_protein_id
            HAVING count(*) > 1
        )
        """
    ).fetchone()[0]
    if duplicate_viral:
        raise ValueError(
            f"Viral taxonomy lookup contains {duplicate_viral} duplicate protein IDs"
        )

    ids_expression = (
        f"read_csv({sql_path(cellular_ids)}, delim='\\t', header=true, "
        "all_varchar=true)"
    )
    missing_metadata = connection.execute(
        f"""
        SELECT count(*)
        FROM {ids_expression} ids
        LEFT JOIN {cellular_expression} metadata USING (protein_id)
        WHERE metadata.protein_id IS NULL
        """
    ).fetchone()[0]
    if missing_metadata:
        raise ValueError(
            f"Cellular metadata is missing {missing_metadata} representative protein IDs"
        )

    unused_metadata = connection.execute(
        f"""
        SELECT count(*)
        FROM {cellular_expression} metadata
        LEFT JOIN {ids_expression} ids USING (protein_id)
        WHERE ids.protein_id IS NULL
        """
    ).fetchone()[0]

    viral_lookup_count = connection.execute(
        f"SELECT count(*) FROM {viral_expression}"
    ).fetchone()[0]
    if viral_lookup_count != viral_count:
        raise ValueError(
            "Viral representative FASTA/taxonomy lookup counts disagree: "
            f"{viral_count} FASTA records versus {viral_lookup_count} lookup rows"
        )

    manifest_sql = f"""
        COPY (
            SELECT
                'VIRAL|' || CAST(v.representative_protein_id AS VARCHAR) AS reference_id,
                'viral' AS reference_class,
                CAST(v.representative_protein_id AS VARCHAR) AS source_protein_id,
                '' AS cellular_group,
                '' AS source_accession,
                '' AS organism_name,
                '' AS taxid,
                CAST(v.taxonomy AS VARCHAR) AS taxonomy,
                CAST(v.classification_rank AS VARCHAR) AS classification_rank,
                CAST(v.taxonomy_conflict AS BOOLEAN) AS taxonomy_conflict,
                CAST(v.taxonomy_conflict_rank AS VARCHAR) AS taxonomy_conflict_rank,
                true AS viral_taxonomy_eligible
            FROM {viral_expression} v
            UNION ALL
            SELECT
                'CELLULAR|' || CAST(ids.protein_id AS VARCHAR) AS reference_id,
                'cellular' AS reference_class,
                CAST(ids.protein_id AS VARCHAR) AS source_protein_id,
                CAST(metadata.cellular_group AS VARCHAR) AS cellular_group,
                CAST(metadata.source_accession AS VARCHAR) AS source_accession,
                {quoted_or_default(cellular_columns, 'organism_name')} AS organism_name,
                {quoted_or_default(cellular_columns, 'taxid')} AS taxid,
                {quoted_or_default(cellular_columns, 'taxonomy')} AS taxonomy,
                {quoted_or_default(cellular_columns, 'classification_rank')} AS classification_rank,
                false AS taxonomy_conflict,
                '' AS taxonomy_conflict_rank,
                false AS viral_taxonomy_eligible
            FROM {ids_expression} ids
            INNER JOIN {cellular_expression} metadata USING (protein_id)
        ) TO {sql_path(args.output_manifest)} (FORMAT PARQUET, COMPRESSION ZSTD)
    """
    connection.execute(manifest_sql)

    manifest_counts = dict(
        connection.execute(
            f"""
            SELECT reference_class, count(*)
            FROM read_parquet({sql_path(args.output_manifest)})
            GROUP BY reference_class
            """
        ).fetchall()
    )
    connection.close()
    if manifest_counts != {"viral": viral_count, "cellular": cellular_count}:
        raise ValueError(
            f"Unexpected reference manifest counts: {manifest_counts}"
        )

    with args.output_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerow(("viral_representatives", viral_count))
        writer.writerow(("cellular_representatives", cellular_count))
        writer.writerow(("combined_references", viral_count + cellular_count))
        writer.writerow(("unused_cellular_metadata_rows", unused_metadata))
        writer.writerow(("reference_labels", "VIRAL|,CELLULAR|"))

    print(f"Viral representatives: {viral_count}")
    print(f"Cellular representatives: {cellular_count}")
    print(f"Combined references: {viral_count + cellular_count}")


if __name__ == "__main__":
    main()
