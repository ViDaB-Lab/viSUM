process PREPARE_VIRBOT_DATABASE {

    tag 'virbot_database'

    // Setup and analysis intentionally share one environment so Nextflow can
    // reuse a single persistent Conda environment for the whole program.
    conda "${projectDir}/envs/virbot.yml"

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'virbot_installation_metadata.tsv'

    input:
    tuple val(installation_path),
          val(installation_source),
          val(auto_download),
          val(repository),
          val(requested_revision)

    output:
    tuple path('virbot_installation'),
          path('virbot_database'),
          path('virbot_installation_metadata.tsv'),
          emit: installation

    script:
    """
    set -euo pipefail

    INSTALL_DEST="${installation_path}"
    INSTALL_SOURCE="${installation_source}"
    AUTO_DOWNLOAD="${auto_download}"
    REPOSITORY="${repository}"
    REQUESTED_REVISION="${requested_revision}"

    report_error() {
        echo "Validation failed: \$1" >&2
    }

    is_lfs_pointer() {
        local candidate="\$1"
        head -n 1 "\$candidate" 2>/dev/null \
            | grep -q '^version https://git-lfs.github.com/spec/'
    }

    validate_installation() {
        local candidate="\$1"
        local report_errors="\${2:-false}"
        local required_file
        local database="\$candidate/virbot/data/ref"

        if [[ ! -d "\$candidate" || ! -r "\$candidate" ]]; then
            [[ "\$report_errors" == 'true' ]] && \
                report_error "VirBot installation is missing or unreadable: \$candidate"
            return 1
        fi

        for required_file in \
            "\$candidate/virbot/VirBot.py" \
            "\$candidate/virbot/__init__.py" \
            "\$database/VirBot.hmm" \
            "\$database/VirBot.dmnd" \
            "\$database/VirBot_hmm_threshold.txt" \
            "\$database/VirBot_hmm_taxa_full.txt" \
            "\$database/VirBot_RNAvirus_acc.txt"; do
            if [[ ! -s "\$required_file" ]]; then
                [[ "\$report_errors" == 'true' ]] && \
                    report_error "missing or empty required file: \$required_file"
                return 1
            fi
            if is_lfs_pointer "\$required_file"; then
                [[ "\$report_errors" == 'true' ]] && \
                    report_error "Git-LFS pointer was not replaced by database content: \$required_file"
                return 1
            fi
        done
    }

    validate_runtime() {
        local candidate="\$1"
        local report_errors="\${2:-false}"
        local database="\$candidate/virbot/data/ref"

        if ! python "\$candidate/virbot/VirBot.py" --help >/dev/null; then
            [[ "\$report_errors" == 'true' ]] && \
                report_error "VirBot could not load its Python code and reference tables"
            return 1
        fi

        if ! hmmstat "\$database/VirBot.hmm" >/dev/null; then
            [[ "\$report_errors" == 'true' ]] && \
                report_error "VirBot.hmm failed HMMER validation"
            return 1
        fi

        if ! diamond dbinfo --db "\$database/VirBot.dmnd" >/dev/null; then
            [[ "\$report_errors" == 'true' ]] && \
                report_error "VirBot.dmnd failed DIAMOND validation"
            return 1
        fi
    }

    validate_bundle() {
        local candidate="\$1"
        local report_errors="\${2:-false}"

        validate_installation "\$candidate" "\$report_errors" || return 1
        validate_runtime "\$candidate" "\$report_errors" || return 1
    }

    installed_revision() {
        git -C "\$1" rev-parse HEAD 2>/dev/null || true
    }

    managed_revision_matches() {
        [[ "\$(installed_revision "\$INSTALL_DEST")" == "\$REQUESTED_REVISION" ]]
    }

    virbot_version() {
        awk -F '"' '/^__version__/ { print \$2; exit }' \
            "\$1/virbot/__init__.py"
    }

    write_metadata() {
        local action="\$1"
        local revision
        local version

        revision="\$(installed_revision "\$INSTALL_DEST")"
        [[ -n "\$revision" ]] || revision='unknown'
        version="\$(virbot_version "\$INSTALL_DEST")"

        printf 'installation_path\tinstallation_source\tdatabase_path\tvalidation\tinstallation_action\trepository\trequested_revision\tinstalled_revision\tvirbot_version\n' \
            > virbot_installation_metadata.tsv
        printf '%s\t%s\t%s\tpassed\t%s\t%s\t%s\t%s\t%s\n' \
            "\$INSTALL_DEST" "\$INSTALL_SOURCE" \
            "\$INSTALL_DEST/virbot/data/ref" "\$action" "\$REPOSITORY" \
            "\$REQUESTED_REVISION" "\$revision" "\$version" \
            >> virbot_installation_metadata.tsv
    }

    emit_installation() {
        local action="\$1"

        ln -s "\$INSTALL_DEST" virbot_installation
        ln -s "\$INSTALL_DEST/virbot/data/ref" virbot_database
        write_metadata "\$action"

        echo "VirBot installation: \$INSTALL_DEST"
        echo "VirBot database: \$INSTALL_DEST/virbot/data/ref"
        echo "Installation source: \$INSTALL_SOURCE"
        echo "Validation: passed"
        echo "Installation action: \$action"
    }

    if [[ "\$INSTALL_SOURCE" == 'user-supplied' ]]; then
        if ! validate_bundle "\$INSTALL_DEST"; then
            echo "ERROR: The user-supplied VirBot installation is incomplete:" >&2
            echo "       \$INSTALL_DEST" >&2
            validate_bundle "\$INSTALL_DEST" true || true
            echo "The path must contain VirBot's code and extracted virbot/data/ref database." >&2
            echo "Automatic setup was not attempted because --virbot_dir was supplied." >&2
            exit 1
        fi

        emit_installation 'skipped-user-supplied-installation'
        exit 0
    fi

    INSTALL_PARENT=\$(dirname "\$INSTALL_DEST")
    mkdir -p "\$INSTALL_PARENT"

    if ! command -v flock >/dev/null 2>&1; then
        echo "ERROR: The 'flock' command is required for safe automatic VirBot setup." >&2
        exit 1
    fi

    exec 9>"\${INSTALL_DEST}.install.lock"
    flock 9

    if validate_bundle "\$INSTALL_DEST"; then
        if managed_revision_matches; then
            emit_installation 'skipped-existing-managed-installation'
            exit 0
        fi

        current_revision="\$(installed_revision "\$INSTALL_DEST")"
        echo "ERROR: The managed VirBot installation has a different revision:" >&2
        echo "       installed: \${current_revision:-unknown}" >&2
        echo "       requested: \$REQUESTED_REVISION" >&2
        echo "viSUM will not replace an existing managed installation automatically." >&2
        exit 1
    fi

    if [[ -e "\$INSTALL_DEST" ]]; then
        echo "ERROR: A managed VirBot path exists but failed validation:" >&2
        echo "       \$INSTALL_DEST" >&2
        validate_bundle "\$INSTALL_DEST" true || true
        echo "viSUM will not overwrite or repair it automatically." >&2
        exit 1
    fi

    if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
        echo "ERROR: No valid managed VirBot installation was found at:" >&2
        echo "       \$INSTALL_DEST" >&2
        echo "Enable --virbot_auto_download true or supply --virbot_dir PATH." >&2
        exit 1
    fi

    STAGE_ROOT="\${INSTALL_PARENT}/.virbot-download.\$\$"
    STAGED_INSTALLATION="\$STAGE_ROOT/VirBot"
    cleanup_stage() {
        rm -rf "\$STAGE_ROOT"
    }
    trap cleanup_stage EXIT

    mkdir -p "\$STAGE_ROOT"
    echo "Fetching the pinned VirBot code and official Git-LFS reference bundle..."
    GIT_LFS_SKIP_SMUDGE=1 git init "\$STAGED_INSTALLATION"
    git -C "\$STAGED_INSTALLATION" remote add origin "\$REPOSITORY"
    GIT_LFS_SKIP_SMUDGE=1 git -C "\$STAGED_INSTALLATION" \
        fetch --depth 1 origin "\$REQUESTED_REVISION"
    GIT_LFS_SKIP_SMUDGE=1 git -C "\$STAGED_INSTALLATION" \
        checkout --detach FETCH_HEAD
    git -C "\$STAGED_INSTALLATION" lfs install --local
    git -C "\$STAGED_INSTALLATION" lfs pull \
        --include='virbot/data/ref.zip' --exclude=''

    (
        cd "\$STAGED_INSTALLATION/virbot/data"
        md5sum --check zip.md5
        unzip -q ref.zip
        # The extracted references are the runtime database. Discard the
        # archive afterward so a managed install does not retain a duplicate.
        rm -f ref.zip
    )
    rm -rf "\$STAGED_INSTALLATION/.git/lfs/objects"

    if ! validate_bundle "\$STAGED_INSTALLATION"; then
        echo "ERROR: The downloaded VirBot installation failed validation." >&2
        validate_bundle "\$STAGED_INSTALLATION" true || true
        exit 1
    fi

    mv "\$STAGED_INSTALLATION" "\$INSTALL_DEST"
    date -Iseconds > "\$INSTALL_DEST/.visum_install_complete"

    emit_installation 'downloaded-and-installed'
    """
}
