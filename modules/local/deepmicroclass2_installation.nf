process PREPARE_DEEPMICROCLASS2 {

    tag 'deepmicroclass2_installation'

    conda 'conda-forge::git'

    publishDir "${params.outdir}/software_setup",
        mode: 'copy',
        pattern: 'deepmicroclass2_installation_metadata.tsv'

    input:
    tuple val(installation_path),
          val(installation_source),
          val(auto_download),
          val(repository),
          val(requested_revision)

    output:
    tuple path('deepmicroclass2_installation'),
          path('deepmicroclass2_installation_metadata.tsv'),
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

    validate_installation() {
        local candidate="\$1"
        local report_errors="\${2:-false}"
        local required_file

        if [[ ! -d "\$candidate" || ! -r "\$candidate" ]]; then
            [[ "\$report_errors" == 'true' ]] && \
                report_error "DeepMicroClass2 installation is missing or unreadable: \$candidate"
            return 1
        fi

        for required_file in \
            "\$candidate/predict.py" \
            "\$candidate/models.py" \
            "\$candidate/utils.py" \
            "\$candidate/model/8-class/DeepMicroClass-best-500-8class.ckpt" \
            "\$candidate/model/8-class/DeepMicroClass-best-1000-8class.ckpt" \
            "\$candidate/model/8-class/DeepMicroClass-best-3000-8class.ckpt"; do
            if [[ ! -s "\$required_file" ]]; then
                [[ "\$report_errors" == 'true' ]] && \
                    report_error "missing or empty required file: \$required_file"
                return 1
            fi
        done
    }

    installed_revision() {
        git -C "\$1" rev-parse HEAD 2>/dev/null || true
    }

    managed_revision_matches() {
        [[ "\$(installed_revision "\$INSTALL_DEST")" == "\$REQUESTED_REVISION" ]]
    }

    write_metadata() {
        local action="\$1"
        local revision

        revision="\$(installed_revision "\$INSTALL_DEST")"
        [[ -n "\$revision" ]] || revision='unknown'

        printf 'installation_path\tinstallation_source\tvalidation\tinstallation_action\trepository\trequested_revision\tinstalled_revision\tmodel_mode\n' \
            > deepmicroclass2_installation_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\t%s\t%s\t%s\n' \
            "\$INSTALL_DEST" "\$INSTALL_SOURCE" "\$action" "\$REPOSITORY" \
            "\$REQUESTED_REVISION" "\$revision" '8class' \
            >> deepmicroclass2_installation_metadata.tsv
    }

    emit_installation() {
        local action="\$1"

        ln -s "\$INSTALL_DEST" deepmicroclass2_installation
        write_metadata "\$action"

        echo "DeepMicroClass2 installation: \$INSTALL_DEST"
        echo "Installation source: \$INSTALL_SOURCE"
        echo "Validation: passed"
        echo "Installation action: \$action"
    }

    if [[ "\$INSTALL_SOURCE" == 'user-supplied' ]]; then
        if ! validate_installation "\$INSTALL_DEST"; then
            echo "ERROR: The user-supplied DeepMicroClass2 installation is incomplete:" >&2
            echo "       \$INSTALL_DEST" >&2
            validate_installation "\$INSTALL_DEST" true || true
            echo "Automatic setup was not attempted because --deepmicroclass2_dir was supplied." >&2
            exit 1
        fi

        emit_installation 'skipped-user-supplied-installation'
        exit 0
    fi

    INSTALL_PARENT=\$(dirname "\$INSTALL_DEST")
    mkdir -p "\$INSTALL_PARENT"

    if ! command -v flock >/dev/null 2>&1; then
        echo "ERROR: The 'flock' command is required for safe automatic DeepMicroClass2 setup." >&2
        exit 1
    fi

    exec 9>"\${INSTALL_DEST}.install.lock"
    flock 9

    if validate_installation "\$INSTALL_DEST"; then
        if managed_revision_matches; then
            emit_installation 'skipped-existing-managed-installation'
            exit 0
        fi

        current_revision="\$(installed_revision "\$INSTALL_DEST")"
        echo "ERROR: The managed DeepMicroClass2 installation has a different revision:" >&2
        echo "       installed: \${current_revision:-unknown}" >&2
        echo "       requested: \$REQUESTED_REVISION" >&2
        echo "viSUM will not replace an existing managed installation automatically." >&2
        exit 1
    fi

    if [[ -e "\$INSTALL_DEST" ]]; then
        echo "ERROR: A managed DeepMicroClass2 path exists but failed validation:" >&2
        echo "       \$INSTALL_DEST" >&2
        validate_installation "\$INSTALL_DEST" true || true
        echo "viSUM will not overwrite or repair it automatically." >&2
        exit 1
    fi

    if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
        echo "ERROR: No valid managed DeepMicroClass2 installation was found at:" >&2
        echo "       \$INSTALL_DEST" >&2
        echo "Enable --deepmicroclass2_auto_download true or supply --deepmicroclass2_dir PATH." >&2
        exit 1
    fi

    STAGE_ROOT="\${INSTALL_PARENT}/.deepmicroclass2-download.\$\$"
    STAGED_INSTALLATION="\$STAGE_ROOT/DeepMicroClass2"
    cleanup_stage() {
        rm -rf "\$STAGE_ROOT"
    }
    trap cleanup_stage EXIT

    mkdir -p "\$STAGE_ROOT"
    echo "Cloning the pinned DeepMicroClass2 code and bundled checkpoints..."
    git clone "\$REPOSITORY" "\$STAGED_INSTALLATION"
    git -C "\$STAGED_INSTALLATION" checkout --detach "\$REQUESTED_REVISION"

    if ! validate_installation "\$STAGED_INSTALLATION"; then
        echo "ERROR: The downloaded DeepMicroClass2 installation failed validation." >&2
        validate_installation "\$STAGED_INSTALLATION" true || true
        exit 1
    fi

    mv "\$STAGED_INSTALLATION" "\$INSTALL_DEST"
    date -Iseconds > "\$INSTALL_DEST/.visum_install_complete"

    emit_installation 'downloaded-and-installed'
    """
}
