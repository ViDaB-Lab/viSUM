process PREPARE_CHECKV_DATABASE {

    tag 'checkv_database'

    // Database setup and the future analysis process will share this pinned
    // environment, so Nextflow only needs to create it once.
    conda "${projectDir}/envs/checkv.yml"

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'checkv_database_metadata.tsv'

    input:
    tuple val(database_path), val(database_source), val(auto_download)

    output:
    tuple path('checkv_database'),
          path('checkv_database_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    DB_DEST="${database_path}"
    DB_SOURCE="${database_source}"
    AUTO_DOWNLOAD="${auto_download}"

    validate_database() {
        local candidate="\$1"
        local hmm_file_count

        [[ -d "\$candidate" && -r "\$candidate" ]] || return 1
        [[ -s "\$candidate/genome_db/checkv_reps.faa" ]] || return 1
        [[ -s "\$candidate/genome_db/checkv_reps.fna" ]] || return 1
        [[ -s "\$candidate/genome_db/checkv_reps.tsv" ]] || return 1
        [[ -s "\$candidate/genome_db/checkv_reps.dmnd" ]] || return 1
        [[ -s "\$candidate/hmm_db/checkv_hmms.tsv" ]] || return 1
        [[ -s "\$candidate/hmm_db/genome_lengths.tsv" ]] || return 1
        [[ -d "\$candidate/hmm_db/checkv_hmms" ]] || return 1

        hmm_file_count=\$(find "\$candidate/hmm_db/checkv_hmms" \
            -maxdepth 1 -type f -size +0c | wc -l)
        [[ "\$hmm_file_count" -gt 0 ]] || return 1

        diamond dbinfo --db "\$candidate/genome_db/checkv_reps.dmnd" \
            >/dev/null 2>&1 || return 1
    }

    checkv_version() {
        checkv --version 2>&1 | head -n 1 || \
            python -c "from importlib.metadata import version; print(version('checkv'))"
    }

    database_release() {
        local candidate="\$1"
        local release

        if [[ -s "\$candidate/database_version" ]]; then
            head -n 1 "\$candidate/database_version"
            return
        fi
        if [[ -s "\$candidate/.visum_database_release" ]]; then
            head -n 1 "\$candidate/.visum_database_release"
            return
        fi

        release=\$(basename "\$candidate")
        printf '%s\n' "\$release"
    }

    write_metadata() {
        local action="\$1"
        local version release
        version=\$(checkv_version)
        release=\$(database_release "\$DB_DEST")

        printf 'database_path\tdatabase_source\tvalidation\tinstallation_action\tdatabase_release\tcheckv_version\n' \
            > checkv_database_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\t%s\n' \
            "\$DB_DEST" "\$DB_SOURCE" "\$action" "\$release" "\$version" \
            >> checkv_database_metadata.tsv
    }

    emit_database() {
        local action="\$1"

        ln -s "\$DB_DEST" checkv_database
        write_metadata "\$action"

        echo "CheckV database: \$DB_DEST"
        echo "Database source: \$DB_SOURCE"
        echo "Validation: passed"
        echo "Installation action: \$action"
    }

    if ! checkv --help >/dev/null 2>&1; then
        echo "ERROR: The pinned CheckV environment failed its CLI self-check." >&2
        exit 1
    fi

    if [[ "\$DB_SOURCE" == 'user-supplied' ]]; then
        if ! validate_database "\$DB_DEST"; then
            echo "ERROR: The user-supplied CheckV database is missing, unreadable, or incomplete:" >&2
            echo "       \$DB_DEST" >&2
            echo "Expected genome_db with the CheckV FASTA/TSV/DIAMOND files and" >&2
            echo "hmm_db with checkv_hmms.tsv, genome_lengths.tsv, and populated checkv_hmms/." >&2
            echo "Automatic download was not attempted because --checkv_db was supplied." >&2
            exit 1
        fi

        emit_database 'skipped-user-supplied-database'
        exit 0
    fi

    DB_PARENT=\$(dirname "\$DB_DEST")
    mkdir -p "\$DB_PARENT"

    if ! command -v flock >/dev/null 2>&1; then
        echo "ERROR: The 'flock' command is required for safe automatic CheckV setup." >&2
        exit 1
    fi

    exec 9>"\${DB_DEST}.install.lock"
    flock 9

    if validate_database "\$DB_DEST"; then
        emit_database 'skipped-existing-managed-database'
        exit 0
    fi

    if [[ -e "\$DB_DEST" ]]; then
        echo "ERROR: A managed CheckV database path exists but failed validation:" >&2
        echo "       \$DB_DEST" >&2
        echo "viSUM will not overwrite or repair it automatically." >&2
        exit 1
    fi

    if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
        echo "ERROR: No valid managed CheckV database was found at:" >&2
        echo "       \$DB_DEST" >&2
        echo "Enable --checkv_auto_download true or supply --checkv_db PATH." >&2
        exit 1
    fi

    STAGE_ROOT="\${DB_PARENT}/.checkv-download.\$\$"
    cleanup_stage() {
        rm -rf "\$STAGE_ROOT"
    }
    trap cleanup_stage EXIT

    mkdir -p "\$STAGE_ROOT"
    echo "Downloading the current official CheckV database into a staged directory..."
    checkv download_database "\$STAGE_ROOT"

    mapfile -t CANDIDATES < <(find "\$STAGE_ROOT" -mindepth 1 -maxdepth 2 \
        -type d -name genome_db -printf '%h\n')
    if [[ "\${#CANDIDATES[@]}" -ne 1 ]]; then
        echo "ERROR: CheckV download produced \${#CANDIDATES[@]} candidate database directories; expected one." >&2
        exit 1
    fi

    STAGED_DB="\${CANDIDATES[0]}"
    if ! validate_database "\$STAGED_DB"; then
        echo "ERROR: The downloaded CheckV database failed validation." >&2
        exit 1
    fi

    RELEASE=\$(basename "\$STAGED_DB")
    printf '%s\n' "\$RELEASE" > "\$STAGED_DB/.visum_database_release"
    date -Iseconds > "\$STAGED_DB/.visum_db_complete"
    mv "\$STAGED_DB" "\$DB_DEST"

    emit_database 'downloaded-and-installed'
    """
}
