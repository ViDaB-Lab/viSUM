process PREPARE_GENOMAD_DATABASE {

    tag 'genomad_database'

    conda "bioconda::genomad=${params.genomad_version}"
    cache false

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'genomad_database_metadata.tsv'

    input:
    tuple val(database_path), val(database_source), val(auto_download)

    output:
    tuple path('genomad_db'),
          path('genomad_database_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    DB_DEST="${database_path}"
    DB_SOURCE="${database_source}"
    AUTO_DOWNLOAD="${auto_download}"

    validate_database() {
        local candidate="\$1"

        [[ -d "\$candidate" ]] || return 1
        [[ -r "\$candidate" ]] || return 1

        python "${projectDir}/bin/validate_reference_database.py" genomad \
            "\$candidate" --check-mmseqs
    }

    write_metadata() {
        local action="\$1"
        local genomad_version
        genomad_version=\$(genomad --version 2>&1 | head -n 1)

        printf 'database_path\tdatabase_source\tvalidation\tinstallation_action\tgenomad_version\n' \
            > genomad_database_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\n' \
            "\$DB_DEST" "\$DB_SOURCE" "\$action" "\$genomad_version" \
            >> genomad_database_metadata.tsv
    }

    if [[ "\$DB_SOURCE" == 'user-supplied' ]]; then
        if ! validate_database "\$DB_DEST"; then
            echo "ERROR: The user-supplied geNomad database is missing, unreadable, or incomplete:" >&2
            echo "       \$DB_DEST" >&2
            echo "Automatic download was not attempted because --genomad_db was supplied." >&2
            exit 1
        fi

        ln -s "\$DB_DEST" genomad_db
        write_metadata 'skipped-user-supplied-database'

        echo "geNomad database: \$DB_DEST"
        echo "Database source: user-supplied"
        echo "Validation: passed"
        echo "Installation: skipped"
        exit 0
    fi

    DB_PARENT=\$(dirname "\$DB_DEST")
    mkdir -p "\$DB_PARENT"

    if ! command -v flock >/dev/null 2>&1; then
        echo "ERROR: The 'flock' command is required for safe automatic database setup." >&2
        exit 1
    fi

    # Only one pipeline run may inspect/install the managed database at a time.
    exec 9>"\${DB_DEST}.install.lock"
    flock 9

    if validate_database "\$DB_DEST"; then
        ln -s "\$DB_DEST" genomad_db
        write_metadata 'skipped-existing-managed-database'

        echo "geNomad database: \$DB_DEST"
        echo "Database source: viSUM-managed"
        echo "Validation: passed"
        echo "Installation: skipped"
        exit 0
    fi

    if [[ -e "\$DB_DEST" ]]; then
        echo "ERROR: A geNomad database path exists but failed validation:" >&2
        echo "       \$DB_DEST" >&2
        echo "viSUM will not overwrite or repair it automatically." >&2
        exit 1
    fi

    if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
        echo "ERROR: No valid managed geNomad database was found at:" >&2
        echo "       \$DB_DEST" >&2
        echo "Enable --genomad_auto_download true or supply --genomad_db PATH." >&2
        exit 1
    fi

    STAGE_ROOT="\${DB_PARENT}/.genomad-download.\$\$"
    cleanup_stage() {
        rm -rf "\$STAGE_ROOT"
    }
    trap cleanup_stage EXIT

    mkdir -p "\$STAGE_ROOT"
    echo "Downloading the geNomad database into a staged directory..."
    genomad download-database "\$STAGE_ROOT"

    STAGED_DB="\$STAGE_ROOT/genomad_db"
    if ! validate_database "\$STAGED_DB"; then
        echo "ERROR: The downloaded geNomad database failed validation." >&2
        exit 1
    fi

    mv "\$STAGED_DB" "\$DB_DEST"
    date -Iseconds > "\$DB_DEST/.visum_db_complete"

    ln -s "\$DB_DEST" genomad_db
    write_metadata 'downloaded-and-installed'

    echo "geNomad database: \$DB_DEST"
    echo "Database source: viSUM-managed"
    echo "Validation: passed"
    echo "Installation: completed"
    """
}
