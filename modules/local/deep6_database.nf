process PREPARE_DEEP6_DATABASE {

    tag 'deep6_database'

    conda 'conda-forge::git'

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'deep6_database_metadata.tsv'

    input:
    tuple val(installation_path),
          val(installation_source),
          val(model_path),
          val(model_source),
          val(auto_download),
          val(repository),
          val(requested_revision)

    output:
    tuple path('deep6_installation'),
          path('deep6_models'),
          path('deep6_database_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    INSTALL_DEST="${installation_path}"
    INSTALL_SOURCE="${installation_source}"
    MODEL_DEST="${model_path}"
    MODEL_SOURCE="${model_source}"
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
                report_error "Deep6 installation is missing or unreadable: \$candidate"
            return 1
        fi

        for required_file in \
            "\$candidate/Master/deep6.py" \
            "\$candidate/Master/deep6_functions.py"; do
            if [[ ! -s "\$required_file" ]]; then
                [[ "\$report_errors" == 'true' ]] && \
                    report_error "missing or empty required file: \$required_file"
                return 1
            fi
        done
    }

    model_for_length() {
        local candidate="\$1"
        local model_length="\$2"

        find "\$candidate" -maxdepth 1 -type f \
            -name "model_\${model_length}_*.h5" -size +0c -print \
            | sort
    }

    validate_hdf5_signature() {
        local model_file="\$1"
        local signature

        signature=\$(head -c 8 "\$model_file" | od -An -tx1 | tr -d '[:space:]')
        [[ "\$signature" == '894844460d0a1a0a' ]]
    }

    validate_models() {
        local candidate="\$1"
        local report_errors="\${2:-false}"
        local model_length
        local model_count
        local model_file

        if [[ ! -d "\$candidate" || ! -r "\$candidate" ]]; then
            [[ "\$report_errors" == 'true' ]] && \
                report_error "Deep6 model directory is missing or unreadable: \$candidate"
            return 1
        fi

        for model_length in 250 500 1000 1500; do
            model_count=\$(model_for_length "\$candidate" "\$model_length" | wc -l)
            if [[ "\$model_count" -ne 1 ]]; then
                if [[ "\$report_errors" == 'true' ]]; then
                    report_error "expected exactly one non-empty model_\${model_length}_*.h5 file in \$candidate; found \$model_count"
                fi
                return 1
            fi

            model_file=\$(model_for_length "\$candidate" "\$model_length")
            if ! validate_hdf5_signature "\$model_file"; then
                [[ "\$report_errors" == 'true' ]] && \
                    report_error "model does not have a valid HDF5 signature: \$model_file"
                return 1
            fi
        done
    }

    validate_bundle() {
        local installation="\$1"
        local models="\$2"
        local report_errors="\${3:-false}"

        validate_installation "\$installation" "\$report_errors" || return 1
        validate_models "\$models" "\$report_errors" || return 1
    }

    managed_revision_matches() {
        local installed_revision

        installed_revision=\$(
            git -C "\$INSTALL_DEST" rev-parse HEAD 2>/dev/null || true
        )
        [[ "\$installed_revision" == "\$REQUESTED_REVISION" ]]
    }

    write_metadata() {
        local action="\$1"
        local installed_revision
        local model_250
        local model_500
        local model_1000
        local model_1500

        installed_revision=\$(
            git -C "\$INSTALL_DEST" rev-parse HEAD 2>/dev/null || true
        )
        [[ -n "\$installed_revision" ]] || installed_revision='unknown'

        model_250=\$(basename "\$(model_for_length "\$MODEL_DEST" 250)")
        model_500=\$(basename "\$(model_for_length "\$MODEL_DEST" 500)")
        model_1000=\$(basename "\$(model_for_length "\$MODEL_DEST" 1000)")
        model_1500=\$(basename "\$(model_for_length "\$MODEL_DEST" 1500)")

        printf 'installation_path\tinstallation_source\tmodel_path\tmodel_source\tvalidation\tinstallation_action\trepository\trequested_revision\tinstalled_revision\tmodel_250\tmodel_500\tmodel_1000\tmodel_1500\n' \
            > deep6_database_metadata.tsv
        printf '%s\t%s\t%s\t%s\tpassed\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "\$INSTALL_DEST" "\$INSTALL_SOURCE" "\$MODEL_DEST" "\$MODEL_SOURCE" \
            "\$action" "\$REPOSITORY" "\$REQUESTED_REVISION" "\$installed_revision" \
            "\$model_250" "\$model_500" "\$model_1000" "\$model_1500" \
            >> deep6_database_metadata.tsv
    }

    emit_database() {
        local action="\$1"

        ln -s "\$INSTALL_DEST" deep6_installation
        ln -s "\$MODEL_DEST" deep6_models
        write_metadata "\$action"

        echo "Deep6 installation: \$INSTALL_DEST"
        echo "Deep6 models: \$MODEL_DEST"
        echo "Installation source: \$INSTALL_SOURCE"
        echo "Model source: \$MODEL_SOURCE"
        echo "Validation: passed"
        echo "Installation action: \$action"
    }

    if [[ "\$MODEL_SOURCE" == 'user-supplied' ]] && ! validate_models "\$MODEL_DEST"; then
        echo "ERROR: The user-supplied Deep6 model directory is incomplete:" >&2
        echo "       \$MODEL_DEST" >&2
        validate_models "\$MODEL_DEST" true || true
        echo "Automatic model download was not attempted because --deep6_model was supplied." >&2
        exit 1
    fi

    if [[ "\$INSTALL_SOURCE" == 'user-supplied' ]]; then
        if ! validate_bundle "\$INSTALL_DEST" "\$MODEL_DEST"; then
            echo "ERROR: The user-supplied Deep6 installation or model directory is incomplete:" >&2
            echo "       installation: \$INSTALL_DEST" >&2
            echo "       models:       \$MODEL_DEST" >&2
            validate_bundle "\$INSTALL_DEST" "\$MODEL_DEST" true || true
            echo "Automatic setup was not attempted because --deep6_dir was supplied." >&2
            exit 1
        fi

        emit_database 'skipped-user-supplied-installation'
        exit 0
    fi

    INSTALL_PARENT=\$(dirname "\$INSTALL_DEST")
    mkdir -p "\$INSTALL_PARENT"

    if ! command -v flock >/dev/null 2>&1; then
        echo "ERROR: The 'flock' command is required for safe automatic Deep6 setup." >&2
        exit 1
    fi

    exec 9>"\${INSTALL_DEST}.install.lock"
    flock 9

    if validate_bundle "\$INSTALL_DEST" "\$MODEL_DEST"; then
        if managed_revision_matches; then
            emit_database 'skipped-existing-managed-installation'
            exit 0
        fi

        installed_revision=\$(
            git -C "\$INSTALL_DEST" rev-parse HEAD 2>/dev/null || true
        )
        echo "ERROR: The managed Deep6 installation is valid but has a different revision:" >&2
        echo "       installed: \${installed_revision:-unknown}" >&2
        echo "       requested: \$REQUESTED_REVISION" >&2
        echo "viSUM will not replace an existing managed installation automatically." >&2
        exit 1
    fi

    if [[ -e "\$INSTALL_DEST" ]]; then
        echo "ERROR: A managed Deep6 path exists but failed validation:" >&2
        echo "       \$INSTALL_DEST" >&2
        validate_bundle "\$INSTALL_DEST" "\$MODEL_DEST" true || true
        echo "viSUM will not overwrite or repair it automatically." >&2
        exit 1
    fi

    if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
        echo "ERROR: No valid managed Deep6 installation was found at:" >&2
        echo "       \$INSTALL_DEST" >&2
        echo "Enable --deep6_auto_download true or supply --deep6_dir PATH." >&2
        exit 1
    fi

    STAGE_ROOT="\${INSTALL_PARENT}/.deep6-download.\$\$"
    STAGED_INSTALLATION="\$STAGE_ROOT/Deep6"
    cleanup_stage() {
        rm -rf "\$STAGE_ROOT"
    }
    trap cleanup_stage EXIT

    mkdir -p "\$STAGE_ROOT"
    echo "Cloning the pinned Deep6 code and pretrained models into a staged directory..."
    git clone "\$REPOSITORY" "\$STAGED_INSTALLATION"
    git -C "\$STAGED_INSTALLATION" checkout --detach "\$REQUESTED_REVISION"

    STAGED_MODELS="\$STAGED_INSTALLATION/Models"
    if [[ "\$MODEL_SOURCE" == 'user-supplied' ]]; then
        STAGED_MODELS="\$MODEL_DEST"
    fi

    if ! validate_bundle "\$STAGED_INSTALLATION" "\$STAGED_MODELS"; then
        echo "ERROR: The downloaded Deep6 bundle failed validation." >&2
        validate_bundle "\$STAGED_INSTALLATION" "\$STAGED_MODELS" true || true
        exit 1
    fi

    mv "\$STAGED_INSTALLATION" "\$INSTALL_DEST"
    date -Iseconds > "\$INSTALL_DEST/.visum_db_complete"

    if [[ "\$MODEL_SOURCE" != 'user-supplied' ]]; then
        MODEL_DEST="\$INSTALL_DEST/Models"
    fi

    emit_database 'downloaded-and-installed'
    """
}
