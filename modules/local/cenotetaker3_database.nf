process PREPARE_CENOTETAKER3_DATABASE {

    tag 'cenotetaker3_database'

    conda "bioconda::cenote-taker3=${params.ct3_version} conda-forge::wget"

    cpus params.threads
    memory params.memory
    time params.time

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'cenotetaker3_database_metadata.tsv'

    input:
    tuple val(database_path),
          val(database_source),
          val(auto_download),
          val(hmm_database_version)

    output:
    tuple path('cenotetaker3_db'),
          path('cenotetaker3_database_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    DB_DEST="${database_path}"
    DB_SOURCE="${database_source}"
    AUTO_DOWNLOAD="${auto_download}"
    HMM_DB_VERSION="${hmm_database_version}"

    validate_mmseqs_database() {
        local database_prefix="\$1"

        [[ -s "\$database_prefix" ]] || return 1
        [[ -s "\${database_prefix}.dbtype" ]] || return 1
    }

    validate_taxonomy_database() {
        local database_prefix="\$1"

        validate_mmseqs_database "\$database_prefix" || return 1
        [[ -s "\${database_prefix}_mapping" ]] || return 1
        [[ -s "\${database_prefix}_taxonomy" ]] || return 1
    }

    validate_database() {
        local candidate="\$1"
        local hmm_name
        local hmm_extension

        [[ -d "\$candidate" ]] || return 1
        [[ -r "\$candidate" ]] || return 1

        for hmm_name in \
            Virion_HMMs \
            DNA_rep_HMMs \
            RDRP_HMMs \
            Useful_Annotation_HMMs \
            phrogs_for_ct; do
            for hmm_extension in h3f h3i h3m h3p; do
                [[ -s "\$candidate/hmmscan_DBs/\$HMM_DB_VERSION/\${hmm_name}.\${hmm_extension}" ]] \
                    || return 1
            done
        done

        validate_taxonomy_database \
            "\$candidate/mmseqs_DBs/ct3_hallmark.taxDB" || return 1
        validate_taxonomy_database \
            "\$candidate/mmseqs_DBs/refseq_virus_prot_taxDB" || return 1
        validate_mmseqs_database "\$candidate/mmseqs_DBs/CDD" || return 1
        [[ -s "\$candidate/viral_cdds_and_pfams_191028.txt" ]] || return 1
    }

    write_metadata() {
        local action="\$1"
        local cenotetaker3_version
        cenotetaker3_version=\$(cenotetaker3 --version 2>&1 | head -n 1 || true)
        if [[ -z "\$cenotetaker3_version" ]]; then
            cenotetaker3_version='unknown'
        fi

        printf 'database_path\tdatabase_source\tvalidation\tinstallation_action\tcenotetaker3_version\thmm_database_version\thhsuite_databases\n' \
            > cenotetaker3_database_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\t%s\tnot-installed\n' \
            "\$DB_DEST" "\$DB_SOURCE" "\$action" "\$cenotetaker3_version" \
            "\$HMM_DB_VERSION" >> cenotetaker3_database_metadata.tsv
    }

    emit_database() {
        local action="\$1"

        ln -s "\$DB_DEST" cenotetaker3_db
        write_metadata "\$action"

        echo "Cenote-Taker 3 database: \$DB_DEST"
        echo "Database source: \$DB_SOURCE"
        echo "Validation: passed"
        echo "Installation action: \$action"
        echo "HHsuite databases: not installed"
    }

    if [[ "\$DB_SOURCE" == 'user-supplied' ]]; then
        if ! validate_database "\$DB_DEST"; then
            echo "ERROR: The user-supplied Cenote-Taker 3 database is missing, unreadable, or incomplete:" >&2
            echo "       \$DB_DEST" >&2
            echo "Expected CT3 HMM database \$HMM_DB_VERSION, hallmark and RefSeq taxonomy databases," >&2
            echo "the MMseqs2 CDD database, and the viral domain list." >&2
            echo "Automatic setup was not attempted because --ct3_db was supplied." >&2
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
        echo "ERROR: A managed Cenote-Taker 3 database path exists but failed validation:" >&2
        echo "       \$DB_DEST" >&2
        echo "viSUM will not overwrite a directory it cannot identify as its own interrupted installation." >&2
        exit 1
    fi

    if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
        echo "ERROR: No valid managed Cenote-Taker 3 database was found at:" >&2
        echo "       \$DB_DEST" >&2
        echo "Enable --ct3_auto_download true or supply --ct3_db PATH." >&2
        exit 1
    fi

    if [[ -f "\$INSTALLING_MARKER" ]]; then
        echo "Removing an interrupted viSUM-managed Cenote-Taker 3 installation before retrying."
        rm -rf "\$DB_DEST"
    fi

    mkdir -p "\$DB_DEST"
    date -Iseconds > "\$INSTALLING_MARKER"

    echo "Installing the required Cenote-Taker 3 databases into the managed path..."
    if ! get_ct3_dbs \
        -o "\$DB_DEST" \
        --hmm T \
        --hallmark_tax T \
        --refseq_tax T \
        --mmseqs_cdd T \
        --domain_list T; then
        echo "ERROR: The Cenote-Taker 3 database downloader returned an error." >&2
        echo "The installation marker was retained so viSUM can safely retry this managed setup." >&2
        exit 1
    fi

    # get_ct3_dbs does not consistently propagate failed download subprocesses,
    # so the installed files must be checked even when the command exits zero.
    if ! validate_database "\$DB_DEST"; then
        echo "ERROR: Cenote-Taker 3 database setup completed, but the database failed validation." >&2
        echo "The installation marker was retained so viSUM can safely retry this managed setup." >&2
        exit 1
    fi

    rm -f "\$INSTALLING_MARKER"
    date -Iseconds > "\$DB_DEST/.visum_db_complete"
    emit_database 'downloaded-and-installed'
    """
}
