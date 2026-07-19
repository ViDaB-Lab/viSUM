process PREPARE_VIRSORTER2_DATABASE {

    tag 'virsorter2_database'

    conda "bioconda::virsorter=${params.virsorter2_version}"

    cpus params.threads
    memory params.memory
    time params.time

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'virsorter2_database_metadata.tsv'

    input:
    tuple val(database_path), val(database_source), val(auto_download)

    output:
    tuple path('virsorter2_db'),
          path('virsorter2_database_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    DB_DEST="${database_path}"
    DB_SOURCE="${database_source}"
    AUTO_DOWNLOAD="${auto_download}"

    validate_database() {
        local candidate="\$1"
        local group

        [[ -d "\$candidate" ]] || return 1
        [[ -r "\$candidate" ]] || return 1
        [[ -s "\$candidate/hmm/viral/combined.hmm" ]] || return 1

        for group in dsDNAphage ssDNA NCLDV RNA lavidaviridae; do
            [[ -s "\$candidate/group/\$group/model" ]] || return 1
        done

        # VirSorter2 uses Snakemake environments installed under the database.
        find "\$candidate/conda_envs" -type f -path '*/conda-meta/history' \
            -print -quit 2>/dev/null | grep -q . || return 1
    }

    write_metadata() {
        local action="\$1"
        local virsorter2_version
        virsorter2_version=\$(virsorter --version 2>&1 | head -n 1)

        printf 'database_path\tdatabase_source\tvalidation\tinstallation_action\tvirsorter2_version\n' \
            > virsorter2_database_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\n' \
            "\$DB_DEST" "\$DB_SOURCE" "\$action" "\$virsorter2_version" \
            >> virsorter2_database_metadata.tsv
    }

    emit_database() {
        local action="\$1"
        ln -s "\$DB_DEST" virsorter2_db
        write_metadata "\$action"

        echo "VirSorter2 database: \$DB_DEST"
        echo "Database source: \$DB_SOURCE"
        echo "Validation: passed"
        echo "Installation action: \$action"
    }

    if [[ "\$DB_SOURCE" == 'user-supplied' ]]; then
        if ! validate_database "\$DB_DEST"; then
            echo "ERROR: The user-supplied VirSorter2 database is missing, unreadable, or incomplete:" >&2
            echo "       \$DB_DEST" >&2
            echo "Expected all five classifier models, the combined viral HMM, and installed dependency environments." >&2
            echo "Automatic setup was not attempted because --virsorter2_db was supplied." >&2
            exit 1
        fi

        emit_database 'skipped-user-supplied-database'
        exit 0
    fi

    DB_PARENT=\$(dirname "\$DB_DEST")
    INSTALLING_MARKER="\$DB_DEST/.visum_installing"
    mkdir -p "\$DB_PARENT"

    if ! command -v flock >/dev/null 2>&1; then
        echo "ERROR: The 'flock' command is required for safe automatic database setup." >&2
        exit 1
    fi

    # Only one pipeline run may inspect or install this managed database at a time.
    exec 9>"\${DB_DEST}.install.lock"
    flock 9

    if validate_database "\$DB_DEST"; then
        rm -f "\$INSTALLING_MARKER"
        emit_database 'skipped-existing-managed-database'
        exit 0
    fi

    if [[ -e "\$DB_DEST" && ! -f "\$INSTALLING_MARKER" ]]; then
        echo "ERROR: A managed VirSorter2 database path exists but failed validation:" >&2
        echo "       \$DB_DEST" >&2
        echo "viSUM will not overwrite a directory it cannot identify as its own interrupted installation." >&2
        exit 1
    fi

    if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
        echo "ERROR: No valid managed VirSorter2 database was found at:" >&2
        echo "       \$DB_DEST" >&2
        echo "Enable --virsorter2_auto_download true or supply --virsorter2_db PATH." >&2
        exit 1
    fi

    if [[ -f "\$INSTALLING_MARKER" ]]; then
        echo "Removing an interrupted viSUM-managed VirSorter2 installation before retrying."
        rm -rf "\$DB_DEST"
    fi

    mkdir -p "\$DB_DEST"
    date -Iseconds > "\$INSTALLING_MARKER"

    echo "Installing the VirSorter2 database and internal dependencies directly into the managed path..."
    virsorter setup -d "\$DB_DEST" -j ${task.cpus}

    if ! validate_database "\$DB_DEST"; then
        echo "ERROR: VirSorter2 setup completed, but the database failed validation." >&2
        echo "The installation marker was retained so viSUM can safely retry this managed setup." >&2
        exit 1
    fi

    rm -f "\$INSTALLING_MARKER"
    date -Iseconds > "\$DB_DEST/.visum_db_complete"
    emit_database 'downloaded-and-installed'
    """
}
