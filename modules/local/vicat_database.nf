process PREPARE_VICAT_DATABASE {

    tag 'vicat_database'

    conda "${projectDir}/envs/vicat_db_build.yml"

    // Validation of an existing database is lightweight. Reserve the large
    // build resources only when raw MetaVR inputs were supplied.
    cpus {
        source_proteins
            ? Math.min(params.vicat_build_cpus as int, params.max_cpus as int)
            : 1
    }
    memory { source_proteins ? params.vicat_build_memory : '4 GB' }
    time { source_proteins ? params.vicat_build_time : '1h' }

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'vicat_database_setup_metadata.tsv'

    input:
    tuple val(database_path),
          val(database_source),
          val(source_proteins),
          val(source_metadata),
          val(protein_sha256),
          val(metadata_sha256),
          val(expected_uvigs),
          val(expected_proteins),
          val(expected_representatives)

    output:
    tuple path('vicat_database'),
          path('vicat_database_setup_metadata.tsv'),
          emit: database

    script:
    def taxonomyMemory = params.vicat_build_memory.toString().replaceAll('\\s+', '')
    """
    set -euo pipefail

    DB_DEST="${database_path}"
    DB_SOURCE="${database_source}"
    SOURCE_PROTEINS="${source_proteins}"
    SOURCE_METADATA="${source_metadata}"

    locate_database_files() {
        local candidate="\$1"

        if [[ -s "\$candidate/IMGVR5_UViG_representatives.dmnd" && \
              -s "\$candidate/IMGVR5_UViG.vicat_taxonomy_lookup.parquet" ]]; then
            VICAT_DIAMOND="\$candidate/IMGVR5_UViG_representatives.dmnd"
            VICAT_TAXONOMY="\$candidate/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
            return 0
        fi

        # Accept the development build layout so an already completed local
        # database can be tested without copying its large files.
        if [[ -s "\$candidate/diamond/IMGVR5_UViG_representatives.dmnd" && \
              -s "\$candidate/taxonomy/IMGVR5_UViG.vicat_taxonomy_lookup.parquet" ]]; then
            VICAT_DIAMOND="\$candidate/diamond/IMGVR5_UViG_representatives.dmnd"
            VICAT_TAXONOMY="\$candidate/taxonomy/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
            return 0
        fi

        return 1
    }

    validate_database() {
        local candidate="\$1"
        local diamond_count lookup_count

        [[ -d "\$candidate" && -r "\$candidate" ]] || return 1
        locate_database_files "\$candidate" || return 1
        diamond dbinfo --db "\$VICAT_DIAMOND" >/dev/null 2>&1 || return 1
        diamond_count=\$(diamond dbinfo --db "\$VICAT_DIAMOND" \
            | awk '\$1 == "Sequences" {print \$2}')
        lookup_count=\$(python - "\$VICAT_TAXONOMY" <<'PY'
import duckdb
import sys

connection = duckdb.connect()
print(connection.execute(
    "SELECT count(*) FROM read_parquet(?)", [sys.argv[1]]
).fetchone()[0])
PY
        )

        [[ "\$diamond_count" == "${expected_representatives}" ]] || return 1
        [[ "\$lookup_count" == "${expected_representatives}" ]] || return 1

        if [[ -f "\$candidate/SHA256SUMS" ]]; then
            (cd "\$candidate" && sha256sum --check --quiet SHA256SUMS) || return 1
        fi
    }

    write_metadata() {
        local action="\$1"
        locate_database_files "\$DB_DEST"
        printf 'database_path\tdatabase_source\tvalidation\tinstallation_action\tdiamond_database\ttaxonomy_lookup\trepresentative_proteins\n' \
            > vicat_database_setup_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\t%s\t%s\n' \
            "\$DB_DEST" "\$DB_SOURCE" "\$action" "\$VICAT_DIAMOND" \
            "\$VICAT_TAXONOMY" "${expected_representatives}" \
            >> vicat_database_setup_metadata.tsv
    }

    emit_database() {
        local action="\$1"
        ln -s "\$DB_DEST" vicat_database
        write_metadata "\$action"
        echo "VICAT_DB database=\$DB_DEST source=\$DB_SOURCE action=\$action"
    }

    if [[ "\$DB_SOURCE" == 'user-supplied' ]]; then
        if ! validate_database "\$DB_DEST"; then
            echo "ERROR: The user-supplied viCAT database failed validation:" >&2
            echo "       \$DB_DEST" >&2
            echo "Expected a DIAMOND database and vOTU-balanced taxonomy lookup" >&2
            echo "containing ${expected_representatives} representative proteins." >&2
            exit 1
        fi
        emit_database 'skipped-user-supplied-database'
        exit 0
    fi

    DB_PARENT=\$(dirname "\$DB_DEST")
    BUILD_DIR="\${DB_DEST}.build"
    mkdir -p "\$DB_PARENT"

    command -v flock >/dev/null 2>&1 || {
        echo "ERROR: flock is required for safe viCAT database setup." >&2
        exit 1
    }
    exec 9>"\${DB_DEST}.install.lock"
    flock 9

    if validate_database "\$DB_DEST"; then
        emit_database 'skipped-existing-managed-database'
        exit 0
    fi

    if [[ -e "\$DB_DEST" ]]; then
        echo "ERROR: The managed viCAT database exists but failed validation:" >&2
        echo "       \$DB_DEST" >&2
        echo "viSUM will not overwrite it automatically." >&2
        exit 1
    fi

    if [[ -z "\$SOURCE_PROTEINS" || -z "\$SOURCE_METADATA" ]]; then
        echo "ERROR: viCAT is enabled, but no valid database is available." >&2
        echo "Supply --vicat_db PATH, or provide both --vicat_metavr_proteins" >&2
        echo "and --vicat_metavr_metadata to build MetaVR5 locally." >&2
        exit 1
    fi
    [[ -f "\$SOURCE_PROTEINS" ]] || {
        echo "ERROR: MetaVR protein file not found: \$SOURCE_PROTEINS" >&2
        exit 1
    }
    [[ -f "\$SOURCE_METADATA" ]] || {
        echo "ERROR: MetaVR metadata file not found: \$SOURCE_METADATA" >&2
        exit 1
    }

    echo "Building the viCAT MetaVR5 database. Temporary files are retained"
    echo "after interruption for restart and removed after successful validation."
    bash "${projectDir}/bin/build_vicat_database.sh" \
        --proteins "\$SOURCE_PROTEINS" \
        --metadata "\$SOURCE_METADATA" \
        --destination "\$DB_DEST" \
        --work-dir "\$BUILD_DIR" \
        --threads "${task.cpus}" \
        --memory-limit "${taxonomyMemory}" \
        --protein-sha256 "${protein_sha256}" \
        --metadata-sha256 "${metadata_sha256}" \
        --expected-uvigs "${expected_uvigs}" \
        --expected-proteins "${expected_proteins}" \
        --expected-representatives "${expected_representatives}"

    validate_database "\$DB_DEST" || {
        echo "ERROR: The completed viCAT database failed final validation." >&2
        exit 1
    }
    emit_database 'built-from-user-supplied-metavr-files'
    """
}
