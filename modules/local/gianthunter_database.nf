process PREPARE_GIANTHUNTER_DATABASE {

    tag 'gianthunter_database'

    // Setup and analysis share one pinned environment so Nextflow creates only
    // one persistent GiantHunter environment.
    conda "${projectDir}/envs/gianthunter.yml"

    cpus params.threads
    memory params.gianthunter_memory
    time params.gianthunter_time

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'gianthunter_database_metadata.tsv'

    input:
    tuple val(database_path),
          val(database_source),
          val(auto_download),
          val(database_release),
          val(database_url),
          val(database_sha256)

    output:
    tuple path('gianthunter_database'),
          path('gianthunter_database_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    DB_DEST="${database_path}"
    DB_SOURCE="${database_source}"
    AUTO_DOWNLOAD="${auto_download}"
    DATABASE_RELEASE="${database_release}"
    DATABASE_URL="${database_url}"
    DATABASE_SHA256="${database_sha256}"

    validate_database() {
        local candidate="\$1"
        local required_file

        [[ -d "\$candidate" && -r "\$candidate" ]] || return 1

        for required_file in \
            database.dmnd \
            names.csv \
            nodes.csv \
            ProkaryoticGroup.pkl \
            proteins.csv \
            RefVirus.dmnd \
            RefVirus.faa \
            RefVirus_anno.pkl \
            taxid.csv \
            transformer.pth; do
            [[ -s "\$candidate/\$required_file" ]] || return 1
        done

        diamond dbinfo --db "\$candidate/database.dmnd" >/dev/null 2>&1 || return 1
        diamond dbinfo --db "\$candidate/RefVirus.dmnd" >/dev/null 2>&1 || return 1
    }

    gianthunter_version() {
        python -c "from importlib.metadata import version; print(version('gianthunter'))" \
            2>/dev/null || true
    }

    write_metadata() {
        local action="\$1"
        local version
        version="\$(gianthunter_version)"
        [[ -n "\$version" ]] || version='unknown'

        printf 'database_path\tdatabase_source\tvalidation\tinstallation_action\tdatabase_release\tdatabase_url\tdatabase_sha256\tgianthunter_version\n' \
            > gianthunter_database_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\t%s\t%s\t%s\n' \
            "\$DB_DEST" "\$DB_SOURCE" "\$action" "\$DATABASE_RELEASE" \
            "\$DATABASE_URL" "\$DATABASE_SHA256" "\$version" \
            >> gianthunter_database_metadata.tsv
    }

    emit_database() {
        local action="\$1"

        ln -s "\$DB_DEST" gianthunter_database
        write_metadata "\$action"

        echo "GiantHunter database: \$DB_DEST"
        echo "Database source: \$DB_SOURCE"
        echo "Validation: passed"
        echo "Installation action: \$action"
    }

    if ! gianthunter --help >/dev/null 2>&1; then
        echo "ERROR: The pinned GiantHunter environment failed its CLI self-check." >&2
        exit 1
    fi

    if [[ "\$DB_SOURCE" == 'user-supplied' ]]; then
        if ! validate_database "\$DB_DEST"; then
            echo "ERROR: The user-supplied GiantHunter database is missing, unreadable, or incomplete:" >&2
            echo "       \$DB_DEST" >&2
            echo "Expected the official GiantHunter v1 database, including both DIAMOND databases," >&2
            echo "taxonomy tables, protein-cluster tables, and transformer.pth." >&2
            echo "Automatic setup was not attempted because --gianthunter_db was supplied." >&2
            exit 1
        fi

        emit_database 'skipped-user-supplied-database'
        exit 0
    fi

    DB_PARENT=\$(dirname "\$DB_DEST")
    INSTALLING_MARKER="\$DB_DEST/.visum_installing"
    mkdir -p "\$DB_PARENT"

    if ! command -v flock >/dev/null 2>&1; then
        echo "ERROR: The 'flock' command is required for safe automatic GiantHunter setup." >&2
        exit 1
    fi

    exec 9>"\${DB_DEST}.install.lock"
    flock 9

    if validate_database "\$DB_DEST"; then
        rm -f "\$INSTALLING_MARKER"
        emit_database 'skipped-existing-managed-database'
        exit 0
    fi

    if [[ -e "\$DB_DEST" && ! -f "\$INSTALLING_MARKER" ]]; then
        echo "ERROR: A managed GiantHunter database path exists but failed validation:" >&2
        echo "       \$DB_DEST" >&2
        echo "viSUM will not overwrite a directory it cannot identify as an interrupted managed installation." >&2
        exit 1
    fi

    if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
        echo "ERROR: No valid managed GiantHunter database was found at:" >&2
        echo "       \$DB_DEST" >&2
        echo "Enable --gianthunter_auto_download true or supply --gianthunter_db PATH." >&2
        exit 1
    fi

    if [[ -f "\$INSTALLING_MARKER" ]]; then
        echo "Removing an interrupted viSUM-managed GiantHunter installation before retrying."
        rm -rf "\$DB_DEST"
    fi

    STAGE_ROOT="\${DB_PARENT}/.gianthunter-download.\$\$"
    ARCHIVE="\$STAGE_ROOT/gianthunter_db_v1.zip"
    STAGED_DB="\$STAGE_ROOT/gianthunter_db_v1"
    cleanup_stage() {
        rm -rf "\$STAGE_ROOT"
    }
    trap cleanup_stage EXIT

    mkdir -p "\$STAGE_ROOT"
    echo "Downloading the official GiantHunter \$DATABASE_RELEASE database..."
    curl --fail --location --retry 3 --output "\$ARCHIVE" "\$DATABASE_URL"
    printf '%s  %s\n' "\$DATABASE_SHA256" "\$ARCHIVE" | sha256sum --check -
    unzip -q "\$ARCHIVE" -d "\$STAGE_ROOT"

    if ! validate_database "\$STAGED_DB"; then
        echo "ERROR: The downloaded GiantHunter database failed validation." >&2
        exit 1
    fi

    date -Iseconds > "\$STAGED_DB/.visum_installing"
    mv "\$STAGED_DB" "\$DB_DEST"
    rm -f "\$INSTALLING_MARKER"
    date -Iseconds > "\$DB_DEST/.visum_db_complete"

    emit_database 'downloaded-and-installed'
    """
}
