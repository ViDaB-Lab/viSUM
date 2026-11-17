process PREPARE_VICAT_NONVIRAL_DATABASE {

    tag 'vicat_nonviral_database'

    conda "${projectDir}/envs/vicat_db_build.yml"
    // Persistent assets must be checked even when Nextflow resumes old tasks.
    cache false
    cpus { (classified_dir || source_manifest) ? Math.min(params.vicat_nonviral_build_cpus as int, params.max_cpus as int) : 1 }
    memory { (classified_dir || source_manifest) ? params.vicat_nonviral_build_memory : '4 GB' }
    time { params.vicat_nonviral_build_time }

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        overwrite: true,
        failOnError: true,
        pattern: 'vicat_nonviral_database_setup_metadata.tsv'

    input:
    tuple val(database_path), val(database_source), val(classified_dir),
          val(source_manifest), val(package_root), val(metadata_root), path(build_helpers)

    output:
    tuple path('vicat_nonviral_database'),
          path('vicat_nonviral_database_setup_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    DB="${database_path}"
    ACTION='reused-existing-database'
    if [[ ! -e "\$DB" && ! -L "\$DB" ]]; then
        [[ -n "${classified_dir}" || -n "${source_manifest}" ]] || {
            echo 'ERROR: No viCAT nonviral database exists. Supply --vicat_nonviral_classified_dir' >&2
            echo 'or enable --vicat_nonviral_auto_download, or supply all three local NCBI source inputs.' >&2
            exit 1
        }
        mkdir -p "\$(dirname "\$DB")"
        command -v flock >/dev/null || { echo 'ERROR: flock is required for nonviral setup.' >&2; exit 1; }
        exec 9>"\${DB}.install.lock"
        flock 9
        if [[ ! -e "\$DB" && ! -L "\$DB" ]]; then
            args=()
            if [[ -n "${classified_dir}" ]]; then
                args=(--classified-dir "${classified_dir}")
            elif [[ -z "${package_root}" && -z "${metadata_root}" ]]; then
                args=(--manifest "${source_manifest}" --download-sources)
            else
                args=(--manifest "${source_manifest}" --package-root "${package_root}" --metadata-root "${metadata_root}")
            fi
            python prepare_vicat_nonviral_build.py "\${args[@]}" \
                --destination "\$DB" --work-dir "\${DB}.build" --threads "${task.cpus}"
            ACTION='built-from-user-supplied-nonviral-files'
            if [[ -n "${source_manifest}" && -z "${package_root}" && -z "${metadata_root}" ]]; then
                ACTION='downloaded-and-built-from-assembly-manifest'
            fi
        fi
    fi
    [[ -d "\$DB" && -r "\$DB" ]] || {
        echo "ERROR: viCAT nonviral database directory is unavailable: \$DB" >&2
        exit 1
    }

    python "${projectDir}/bin/validate_reference_database.py" vicat_nonviral "\$DB"

    diamond dbinfo --db "\$DB/vicat_nonviral.dmnd" >/dev/null
    DIAMOND_COUNT=\$(diamond dbinfo --db "\$DB/vicat_nonviral.dmnd" \
        | python3 "${projectDir}/bin/parse_diamond_dbinfo.py")
    read -r METADATA_COUNT INVALID_CLASSES INVALID_FLANKS DUPLICATES <<< "\$(
        python - "\$DB/vicat_nonviral_representative_metadata.parquet" <<'PY'
import duckdb
import sys

metadata = sys.argv[1]
connection = duckdb.connect()
metadata_count = connection.execute(
    "SELECT count(*) FROM read_parquet(?)", [metadata]
).fetchone()[0]
allowed = (
    "CELLULAR_CHROMOSOME", "CELLULAR_UNPLACED", "PLASMID", "PLASTID",
    "MITOCHONDRIAL", "SHARED_NONVIRAL",
)
invalid_classes = connection.execute(
    "SELECT count(*) FROM read_parquet(?) WHERE reference_class IS NULL OR reference_class NOT IN (?, ?, ?, ?, ?, ?)",
    [metadata, *allowed],
).fetchone()[0]
invalid_flanks = connection.execute(
    "SELECT count(*) FROM read_parquet(?) "
    "WHERE provirus_flank_eligible IS DISTINCT FROM "
    "(reference_class = 'CELLULAR_CHROMOSOME')",
    [metadata],
).fetchone()[0]
duplicates = connection.execute(
    "SELECT count(*) FROM ("
    "SELECT reference_id FROM read_parquet(?) "
    "GROUP BY reference_id HAVING count(*) != 1 OR reference_id IS NULL OR reference_id = '')",
    [metadata],
).fetchone()[0]
print(metadata_count, invalid_classes, invalid_flanks, duplicates)
PY
    )"

    [[ "\$DIAMOND_COUNT" == "\$METADATA_COUNT" && "\$METADATA_COUNT" -gt 0 ]] || {
        echo "ERROR: Nonviral DIAMOND and metadata counts disagree" >&2
        exit 1
    }
    [[ "\$INVALID_CLASSES" == 0 && "\$INVALID_FLANKS" == 0 && "\$DUPLICATES" == 0 ]] || {
        echo "ERROR: Nonviral reference metadata invariants failed" >&2
        exit 1
    }
    (cd "\$DB" && sha256sum --check --quiet SHA256SUMS)

    ln -s "\$DB" vicat_nonviral_database
    printf 'database_path\tdatabase_source\tvalidation\tdiamond_sequences\tmetadata_rows\tinstallation_action\n' \
        > vicat_nonviral_database_setup_metadata.tsv
    printf '%s\t%s\tpassed\t%s\t%s\t%s\n' \
        "\$DB" "${database_source}" "\$DIAMOND_COUNT" "\$METADATA_COUNT" "\$ACTION" \
        >> vicat_nonviral_database_setup_metadata.tsv

    echo "VICAT_NONVIRAL_DB database=\$DB sequences=\$DIAMOND_COUNT"
    """
}
