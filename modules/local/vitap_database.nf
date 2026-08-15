process PREPARE_VITAP_DATABASE {

    tag 'vitap_database'

    conda "${projectDir}/envs/vitap.yml"

    cpus { Math.min(params.vitap_cpus as int, params.max_cpus as int) }
    memory { params.vitap_memory }
    time { params.vitap_time }

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'vitap_database_metadata.tsv'

    input:
    tuple val(database_path),
          val(database_source),
          val(vmr_source),
          val(database_label),
          val(auto_download),
          val(update_database),
          val(cleanup_source),
          path(vmr_helper),
          path(update_helper)

    output:
    tuple path('vitap_database'),
          path('vitap_database_metadata.tsv'),
          emit: database

    script:
    """
    set -euo pipefail

    DATABASE_PATH="${database_path}"
    DATABASE_SOURCE="${database_source}"
    VMR_SOURCE="${vmr_source}"
    REQUESTED_LABEL="${database_label}"
    AUTO_DOWNLOAD="${auto_download}"
    UPDATE_DATABASE="${update_database}"
    CLEANUP_SOURCE="${cleanup_source}"
    ICTV_CURRENT_URL="${params.vitap_current_vmr_url}"

    validate_database() {
        local candidate="\$1"
        local gff_file db_base rank

        [[ -d "\$candidate" && -r "\$candidate" ]] || return 1
        mapfile -t GFF_FILES < <(find "\$candidate" -maxdepth 1 -type f \
            -name 'VMR_genome_*.gff' -size +0c | sort)
        [[ "\${#GFF_FILES[@]}" -eq 1 ]] || return 1
        gff_file="\${GFF_FILES[0]}"
        db_base="\${gff_file%.gff}"

        [[ -s "\${db_base}.fasta" ]] || return 1
        [[ -s "\${db_base}.faa" ]] || return 1
        [[ -s "\${db_base}.dmnd" ]] || return 1
        [[ -s "\$candidate/uniref90.dmnd" ]] || return 1
        [[ -s "\$candidate/uniref90.accession2taxid" ]] || return 1
        [[ -s "\$candidate/taxdmp/nodes.dmp" ]] || return 1
        [[ -s "\$candidate/taxdmp/names.dmp" ]] || return 1

        for rank in Realm Kingdom Phylum Class Order Family Genus Species; do
            [[ -s "\$candidate/\${rank}_genome.threshold" ]] || return 1
        done

        find "\$candidate" -maxdepth 1 -type f -name '*.csv' -size +0c \
            | grep -q . || return 1
        diamond dbinfo --db "\${db_base}.dmnd" >/dev/null 2>&1 || return 1
        diamond dbinfo --db "\$candidate/uniref90.dmnd" >/dev/null 2>&1 || return 1
    }

    vitap_version() {
        python -c "from importlib.metadata import version; print(version('VITAP'))"
    }

    write_metadata() {
        local resolved_db="\$1"
        local action="\$2"
        local vmr_file="\$3"
        local vmr_url="\$4"
        local vmr_sha="\$5"
        local release="\$6"
        local version
        version=\$(vitap_version)

        printf 'database_path\tdatabase_source\tvalidation\tinstallation_action\tdatabase_release\tvmr_file\tvmr_source\tvmr_sha256\tvitap_version\trequested_cpus\tupstream_update_cpu_control\n' \
            > vitap_database_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "\$resolved_db" "\$DATABASE_SOURCE" "\$action" "\$release" \
            "\$vmr_file" "\$vmr_url" "\$vmr_sha" "\$version" \
            "${task.cpus}" 'explicit-wrapper-injection' \
            >> vitap_database_metadata.tsv
    }

    emit_database() {
        local resolved_db="\$1"
        local action="\$2"
        local vmr_file="\${3:-unknown}"
        local vmr_url="\${4:-unknown}"
        local vmr_sha="\${5:-unknown}"
        local release="\${6:-unknown}"

        ln -s "\$resolved_db" vitap_database
        write_metadata "\$resolved_db" "\$action" "\$vmr_file" \
            "\$vmr_url" "\$vmr_sha" "\$release"
        echo "VITAP_DB database=\$resolved_db source=\$DATABASE_SOURCE action=\$action"
    }

    emit_existing_managed_database() {
        local resolved_db="\$1"
        local action="\$2"
        local release vmr_file vmr_url vmr_sha recorded
        release=\$(basename "\$resolved_db" | sed 's/^DB_//')
        recorded="\$resolved_db/.visum_vmr_metadata.tsv"
        if [[ -s "\$recorded" ]]; then
            vmr_file=\$(awk -F '\t' 'NR == 2 {print \$2}' "\$recorded")
            vmr_url=\$(awk -F '\t' 'NR == 2 {print \$3}' "\$recorded")
            vmr_sha=\$(awk -F '\t' 'NR == 2 {print \$4}' "\$recorded")
        else
            vmr_file='unknown-legacy-managed-database'
            vmr_url='unknown-legacy-managed-database'
            vmr_sha='unknown-legacy-managed-database'
        fi
        emit_database "\$resolved_db" "\$action" "\$vmr_file" \
            "\$vmr_url" "\$vmr_sha" "\$release"
    }

    if ! VITAP --help >/dev/null 2>&1; then
        echo 'ERROR: The pinned VITAP 1.12 environment failed its CLI self-check.' >&2
        exit 1
    fi

    if [[ "\$DATABASE_SOURCE" == 'user-supplied' ]]; then
        RESOLVED_DATABASE=\$(readlink -f "\$DATABASE_PATH")
        if ! validate_database "\$RESOLVED_DATABASE"; then
            echo 'ERROR: The user-supplied VITAP database failed validation:' >&2
            echo "       \$DATABASE_PATH" >&2
            echo 'Expected the VMR FASTA/FAA/GFF/DIAMOND files, all eight rank' >&2
            echo 'thresholds, uniref90.dmnd, uniref90.accession2taxid, and taxdmp/.' >&2
            exit 1
        fi
        emit_existing_managed_database "\$RESOLVED_DATABASE" \
            'skipped-user-supplied-database'
        exit 0
    fi

    DB_ROOT="\$DATABASE_PATH"
    mkdir -p "\$DB_ROOT/sources"
    command -v flock >/dev/null 2>&1 || {
        echo 'ERROR: flock is required for safe automatic VITAP setup.' >&2
        exit 1
    }
    exec 9>"\$DB_ROOT/vitap_db.install.lock"
    flock 9

    if [[ "\$UPDATE_DATABASE" != 'true' && -L "\$DB_ROOT/current" ]]; then
        RESOLVED_DATABASE=\$(readlink -f "\$DB_ROOT/current")
        if validate_database "\$RESOLVED_DATABASE"; then
            emit_existing_managed_database "\$RESOLVED_DATABASE" \
                'skipped-existing-managed-database'
            exit 0
        fi
        echo 'ERROR: The managed VITAP current link points to an invalid database:' >&2
        echo "       \$RESOLVED_DATABASE" >&2
        exit 1
    fi

    if [[ -n "\$VMR_SOURCE" ]]; then
        [[ -f "\$VMR_SOURCE" ]] || {
            echo "ERROR: User-supplied VMR not found: \$VMR_SOURCE" >&2
            exit 1
        }
        SOURCE_FILE="\$DB_ROOT/sources/\$(basename "\$VMR_SOURCE")"
        if [[ "\$(readlink -f "\$VMR_SOURCE")" != "\$(readlink -m "\$SOURCE_FILE")" ]]; then
            cp -f "\$VMR_SOURCE" "\$SOURCE_FILE"
        fi
        SOURCE_URL='user-supplied'
    else
        if [[ "\$AUTO_DOWNLOAD" != 'true' ]]; then
            echo 'ERROR: No valid managed VITAP database or VMR source is available.' >&2
            echo 'Enable --vitap_auto_download true, supply --vitap_vmr PATH,' >&2
            echo 'or supply a completed database with --vitap_db PATH.' >&2
            exit 1
        fi
        DOWNLOAD_TMP="\$DB_ROOT/sources/.current-vmr.\$\$.xlsx"
        SOURCE_URL=\$(curl --fail --location --retry 5 --continue-at - \
            --output "\$DOWNLOAD_TMP" --write-out '%{url_effective}' \
            "\$ICTV_CURRENT_URL")
        SOURCE_NAME=\$(basename "\${SOURCE_URL%%\?*}")
        [[ "\$SOURCE_NAME" == *.xlsx ]] || {
            echo "ERROR: ICTV current-VMR URL resolved to an unexpected file: \$SOURCE_URL" >&2
            exit 1
        }
        SOURCE_FILE="\$DB_ROOT/sources/\$SOURCE_NAME"
        if [[ -s "\$SOURCE_FILE" ]]; then
            rm -f "\$DOWNLOAD_TMP"
        else
            mv "\$DOWNLOAD_TMP" "\$SOURCE_FILE"
        fi
    fi

    SOURCE_STEM=\$(basename "\$SOURCE_FILE")
    SOURCE_STEM="\${SOURCE_STEM%.*}"
    PREPARED_CSV="\$DB_ROOT/sources/\${SOURCE_STEM}.vitap.csv"
    PREPARED_METADATA="\$DB_ROOT/sources/\${SOURCE_STEM}.vitap_vmr_metadata.tsv"

    PREPARE_ARGS=(
        --input "\$SOURCE_FILE"
        --output-csv "\$PREPARED_CSV"
        --metadata "\$PREPARED_METADATA"
        --source-url "\$SOURCE_URL"
    )
    if [[ -n "\$REQUESTED_LABEL" ]]; then
        PREPARE_ARGS+=(--label "\$REQUESTED_LABEL")
    fi
    [[ -s prepare_vitap_vmr.py ]] || {
        echo 'ERROR: Nextflow did not stage prepare_vitap_vmr.py.' >&2
        exit 1
    }
    [[ -s run_vitap_update.py ]] || {
        echo 'ERROR: Nextflow did not stage run_vitap_update.py.' >&2
        exit 1
    }
    python prepare_vitap_vmr.py "\${PREPARE_ARGS[@]}"

    RELEASE=\$(awk -F '\t' 'NR == 2 {print \$1}' "\$PREPARED_METADATA")
    VMR_SHA=\$(awk -F '\t' 'NR == 2 {print \$4}' "\$PREPARED_METADATA")
    TARGET_DATABASE="\$DB_ROOT/DB_\$RELEASE"

    if validate_database "\$TARGET_DATABASE"; then
        ln -sfn "DB_\$RELEASE" "\$DB_ROOT/current"
        emit_database "\$TARGET_DATABASE" 'skipped-current-release-already-installed' \
            "\$(basename "\$SOURCE_FILE")" "\$SOURCE_URL" "\$VMR_SHA" "\$RELEASE"
        exit 0
    fi
    if [[ -e "\$TARGET_DATABASE" ]]; then
        echo 'ERROR: A final VITAP database directory exists but failed validation:' >&2
        echo "       \$TARGET_DATABASE" >&2
        echo 'Move it aside or remove it after inspecting the incomplete contents.' >&2
        exit 1
    fi

    BUILD_ROOT="\$DB_ROOT/.build-\$RELEASE"
    mkdir -p "\$BUILD_ROOT"
    if [[ ! -e "\$BUILD_ROOT/VMR_Genome" ]]; then
        mkdir -p "\$DB_ROOT/VMR_Genome"
        ln -s "\$DB_ROOT/VMR_Genome" "\$BUILD_ROOT/VMR_Genome"
    fi

    echo 'Building the VITAP 1.12 database. This includes UniRef90 and can use'
    echo 'roughly 200-250 GB of temporary plus final storage. Interrupted build'
    echo 'files are retained under the managed database root for diagnosis/retry.'
    export OMP_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export POLARS_MAX_THREADS="\$OMP_NUM_THREADS"
    PROCESS_DIR="\$PWD"
    (
        cd "\$BUILD_ROOT"
        printf 'Y\n' | python "\$PROCESS_DIR/run_vitap_update.py" \
            --vmr "\$PREPARED_CSV" \
            --out "\$BUILD_ROOT/\${SOURCE_STEM}_reformat.csv" \
            --db "\$RELEASE" \
            --threads "${task.cpus}"
    )

    STAGED_DATABASE="\$BUILD_ROOT/DB_\$RELEASE"
    if ! validate_database "\$STAGED_DATABASE"; then
        echo 'ERROR: VITAP reported completion, but the database failed validation.' >&2
        exit 1
    fi
    cp "\$PREPARED_METADATA" "\$STAGED_DATABASE/.visum_vmr_metadata.tsv"
    date -Iseconds > "\$STAGED_DATABASE/.visum_db_complete"
    mv "\$STAGED_DATABASE" "\$TARGET_DATABASE"
    ln -sfn "DB_\$RELEASE" "\$DB_ROOT/current"

    if [[ "\$CLEANUP_SOURCE" == 'true' ]]; then
        rm -f "\$TARGET_DATABASE/uniref90.fasta" \
              "\$TARGET_DATABASE/uniref90.fasta.gz" \
              "\$TARGET_DATABASE/taxdmp.zip"
    fi
    rm -f "\$BUILD_ROOT/VMR_Genome"
    rmdir "\$BUILD_ROOT" 2>/dev/null || true

    emit_database "\$TARGET_DATABASE" 'downloaded-and-built' \
        "\$(basename "\$SOURCE_FILE")" "\$SOURCE_URL" "\$VMR_SHA" "\$RELEASE"
    """
}
