process PREPARE_VICAT_NONVIRAL_DATABASE {

    tag 'vicat_nonviral_database'

    conda "${projectDir}/envs/vicat.yml"
    cpus 1
    memory '4 GB'
    time '2h'

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'vicat_nonviral_database_setup_metadata.tsv'

    input:
    tuple val(database_path), val(database_source)

    output:
    tuple path('vicat_nonviral_database'),
          path('vicat_nonviral_database_setup_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    DB="${database_path}"
    [[ -d "\$DB" && -r "\$DB" ]] || {
        echo "ERROR: viCAT nonviral database directory is unavailable: \$DB" >&2
        exit 1
    }

    REQUIRED=(
        vicat_nonviral.dmnd
        vicat_nonviral_representatives.faa.gz
        vicat_nonviral_representative_metadata.parquet
        vicat_nonviral_database_metadata.tsv
        vicat_nonviral_cluster_summary.tsv
        SHA256SUMS
        .visum_db_complete
    )
    for name in "\${REQUIRED[@]}"; do
        [[ -s "\$DB/\$name" ]] || {
            echo "ERROR: viCAT nonviral database is missing: \$DB/\$name" >&2
            exit 1
        }
    done

    diamond dbinfo --db "\$DB/vicat_nonviral.dmnd" >/dev/null
    read -r DIAMOND_COUNT METADATA_COUNT INVALID_CLASSES INVALID_FLANKS DUPLICATES <<< "\$(
        python - "\$DB/vicat_nonviral.dmnd" \
            "\$DB/vicat_nonviral_representative_metadata.parquet" <<'PY'
import duckdb
import subprocess
import sys

database, metadata = sys.argv[1:]
dbinfo = subprocess.run(
    ["diamond", "dbinfo", "--db", database],
    check=True,
    capture_output=True,
    text=True,
).stdout.splitlines()
diamond_count = next(
    int(line.split()[-1]) for line in dbinfo if line.strip().startswith("Sequences")
)
connection = duckdb.connect()
metadata_count = connection.execute(
    "SELECT count(*) FROM read_parquet(?)", [metadata]
).fetchone()[0]
allowed = (
    "CELLULAR_CHROMOSOME", "CELLULAR_UNPLACED", "PLASMID", "PLASTID",
    "MITOCHONDRIAL", "SHARED_NONVIRAL",
)
invalid_classes = connection.execute(
    "SELECT count(*) FROM read_parquet(?) WHERE reference_class NOT IN (?, ?, ?, ?, ?, ?)",
    [metadata, *allowed],
).fetchone()[0]
invalid_flanks = connection.execute(
    "SELECT count(*) FROM read_parquet(?) "
    "WHERE provirus_flank_eligible != "
    "(reference_class = 'CELLULAR_CHROMOSOME')",
    [metadata],
).fetchone()[0]
duplicates = connection.execute(
    "SELECT count(*) FROM ("
    "SELECT reference_id FROM read_parquet(?) "
    "GROUP BY reference_id HAVING count(*) != 1)",
    [metadata],
).fetchone()[0]
print(diamond_count, metadata_count, invalid_classes, invalid_flanks, duplicates)
PY
    )"

    [[ "\$DIAMOND_COUNT" == "\$METADATA_COUNT" ]] || {
        echo "ERROR: Nonviral DIAMOND and metadata counts disagree" >&2
        exit 1
    }
    [[ "\$INVALID_CLASSES" == 0 && "\$INVALID_FLANKS" == 0 && "\$DUPLICATES" == 0 ]] || {
        echo "ERROR: Nonviral reference metadata invariants failed" >&2
        exit 1
    }
    (cd "\$DB" && sha256sum --check --quiet SHA256SUMS)

    ln -s "\$DB" vicat_nonviral_database
    printf 'database_path\tdatabase_source\tvalidation\tdiamond_sequences\tmetadata_rows\n' \
        > vicat_nonviral_database_setup_metadata.tsv
    printf '%s\t%s\tpassed\t%s\t%s\n' \
        "\$DB" "${database_source}" "\$DIAMOND_COUNT" "\$METADATA_COUNT" \
        >> vicat_nonviral_database_setup_metadata.tsv

    echo "VICAT_NONVIRAL_DB database=\$DB sequences=\$DIAMOND_COUNT"
    """
}
