process PREPARE_VCONTACT3_DATABASE {

    tag 'vcontact3_database'

    conda "${projectDir}/envs/vcontact3.yml"

    cpus { Math.min(params.vcontact3_cpus as int, params.max_cpus as int) }
    memory { params.vcontact3_memory }
    time { params.vcontact3_time }

    publishDir "${params.outdir}/database_setup",
        mode: 'copy',
        pattern: 'vcontact3_database_metadata.tsv'

    input:
    tuple val(database_path),
          val(database_source),
          val(auto_download),
          val(update_database),
          val(cleanup_archive),
          path(database_validator)

    output:
    tuple path('vcontact3_database'),
          path('vcontact3_database_metadata.tsv'),
          emit: database

    // A shell block keeps Bash variable expansion separate from Nextflow's
    // !{...} interpolation and avoids the escaping failures seen in older
    // embedded setup processes.
    shell:
    '''
    set -euo pipefail

    DATABASE_PATH="!{database_path}"
    DATABASE_SOURCE="!{database_source}"
    AUTO_DOWNLOAD="!{auto_download}"
    UPDATE_DATABASE="!{update_database}"
    CLEANUP_ARCHIVE="!{cleanup_archive}"
    REQUESTED_VERSION='latest'
    VALIDATOR="$PWD/validate_vcontact3_database.py"

    [[ -s "$VALIDATOR" ]] || {
        echo 'ERROR: Nextflow did not stage validate_vcontact3_database.py.' >&2
        exit 1
    }

    vcontact3_version() {
        python -c "from importlib.metadata import version; print(version('vcontact3'))"
    }

    require_compatible_runtime() {
        python -c "from importlib.metadata import version; from packaging.version import Version; current=Version(version('vcontact3')); assert current >= Version('3.2.0'), f'vConTACT3 {current} is incompatible with the current v232+ database; version 3.2.0 or newer is required'"
        command -v mmseqs >/dev/null 2>&1 || {
            echo 'ERROR: The vConTACT3 environment does not provide MMseqs2.' >&2
            exit 1
        }
    }

    validate_database() {
        python "$VALIDATOR" "$1" --field summary >/dev/null 2>&1
    }

    installed_version() {
        python "$VALIDATOR" "$1" --field version
    }

    manifest_path() {
        python "$VALIDATOR" "$1" --field manifest
    }

    database_domains() {
        python "$VALIDATOR" "$1" --field domains
    }

    checked_file_count() {
        python "$VALIDATOR" "$1" --field checked_files
    }

    archive_status() {
        local root="$1"
        local version="$2"
        if find "$root" -maxdepth 1 -type f -name "v${version}.tar.*" -size +0c | grep -q .; then
            printf 'retained\n'
        else
            printf 'not-retained\n'
        fi
    }

    write_metadata() {
        local resolved_path="$1"
        local action="$2"
        local version manifest domains checked archive_state software_version
        version=$(installed_version "$resolved_path")
        manifest=$(manifest_path "$resolved_path")
        domains=$(database_domains "$resolved_path")
        checked=$(checked_file_count "$resolved_path")
        archive_state=$(archive_status "$(dirname "$manifest")" "$version")
        software_version=$(vcontact3_version)

        printf 'database_path\tdatabase_source\tvalidation\tinstallation_action\trequested_version\tinstalled_version\tmanifest_path\tdomains\tvalidated_file_count\tarchive_status\tvcontact3_version\tupdate_requested\n' > vcontact3_database_metadata.tsv
        printf '%s\t%s\tpassed\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$resolved_path" "$DATABASE_SOURCE" "$action" \
            "$REQUESTED_VERSION" "$version" "$manifest" "$domains" \
            "$checked" "$archive_state" "$software_version" \
            "$UPDATE_DATABASE" >> vcontact3_database_metadata.tsv
    }

    emit_database() {
        local resolved_path="$1"
        local action="$2"
        local version
        version=$(installed_version "$resolved_path")
        ln -s "$resolved_path" vcontact3_database
        write_metadata "$resolved_path" "$action"
        echo "VCONTACT3_DB database=$resolved_path version=$version action=$action"
    }

    require_compatible_runtime

    if [[ "$DATABASE_SOURCE" == 'user-supplied' ]]; then
        RESOLVED_DATABASE=$(readlink -f "$DATABASE_PATH")
        if ! validate_database "$RESOLVED_DATABASE"; then
            echo 'ERROR: The user-supplied vConTACT3 database failed validation:' >&2
            echo "       $DATABASE_PATH" >&2
            echo 'Expected an official three-digit JSON manifest, both prokaryotic' >&2
            echo 'and eukaryotic RefSeq data, identity tables, and the v232+ VOGDB files.' >&2
            exit 1
        fi
        emit_database "$RESOLVED_DATABASE" 'skipped-user-supplied-database'
        exit 0
    fi

    DB_ROOT="$DATABASE_PATH"
    mkdir -p "$DB_ROOT/releases"
    command -v flock >/dev/null 2>&1 || {
        echo 'ERROR: flock is required for safe automatic vConTACT3 setup.' >&2
        exit 1
    }
    exec 9>"$DB_ROOT/vcontact3_db.install.lock"
    flock 9

    CURRENT_DATABASE=''
    if [[ -L "$DB_ROOT/current" ]]; then
        CURRENT_DATABASE=$(readlink -f "$DB_ROOT/current")
    fi

    if [[ -n "$CURRENT_DATABASE" ]] && \
       validate_database "$CURRENT_DATABASE" && \
       [[ "$UPDATE_DATABASE" != 'true' ]]; then
        emit_database "$CURRENT_DATABASE" 'skipped-existing-managed-database'
        exit 0
    fi

    if [[ -n "$CURRENT_DATABASE" ]] && ! validate_database "$CURRENT_DATABASE"; then
        echo 'ERROR: The managed vConTACT3 current link is invalid:' >&2
        echo "       $CURRENT_DATABASE" >&2
        echo 'Inspect or remove the link after reviewing the incomplete release.' >&2
        exit 1
    fi

    if [[ "$AUTO_DOWNLOAD" != 'true' ]]; then
        if [[ -n "$CURRENT_DATABASE" ]] && validate_database "$CURRENT_DATABASE"; then
            emit_database "$CURRENT_DATABASE" 'skipped-existing-managed-database-update-disabled'
            exit 0
        fi
        echo 'ERROR: No valid managed vConTACT3 database was found at:' >&2
        echo "       $DB_ROOT" >&2
        echo 'Enable --vcontact3_auto_download true or supply --vcontact3_db PATH.' >&2
        exit 1
    fi

    # Stage beside the managed root so the validated directory can be renamed
    # atomically into releases/<version> on the same filesystem.
    STAGE_ROOT="${DB_ROOT}.download-$$"
    cleanup_stage() {
        rm -rf "$STAGE_ROOT"
    }
    trap cleanup_stage EXIT
    mkdir -p "$STAGE_ROOT"

    echo 'Downloading the latest compatible official vConTACT3 database...'
    vcontact3 prepare_databases \
        --get-version latest \
        --set-location "$STAGE_ROOT"

    if ! validate_database "$STAGE_ROOT"; then
        echo 'ERROR: The downloaded vConTACT3 database failed validation.' >&2
        exit 1
    fi

    VERSION=$(installed_version "$STAGE_ROOT")
    MANIFEST=$(manifest_path "$STAGE_ROOT")
    VERSION_DIRECTORY="$STAGE_ROOT/v$VERSION"
    [[ -d "$VERSION_DIRECTORY" ]] || {
        echo "ERROR: The downloaded database is missing v$VERSION/." >&2
        exit 1
    }

    FINAL_RELEASE="$DB_ROOT/releases/$VERSION"
    if [[ -e "$FINAL_RELEASE" ]]; then
        if validate_database "$FINAL_RELEASE" && \
           [[ "$(installed_version "$FINAL_RELEASE")" == "$VERSION" ]]; then
            ln -sfn "releases/$VERSION" "$DB_ROOT/current"
            emit_database "$FINAL_RELEASE" 'skipped-latest-managed-database-already-installed'
            exit 0
        fi
        echo "ERROR: Managed vConTACT3 release $VERSION exists but failed validation." >&2
        echo 'Inspect or move the incomplete release before retrying.' >&2
        exit 1
    fi

    if [[ "$CLEANUP_ARCHIVE" == 'true' ]]; then
        find "$STAGE_ROOT" -maxdepth 1 -type f -name "v${VERSION}.tar.*" -delete
    fi

    date -Iseconds > "$STAGE_ROOT/.visum_db_complete"
    printf '%s\n' "$VERSION" > "$STAGE_ROOT/.visum_database_release"
    mv "$STAGE_ROOT" "$FINAL_RELEASE"
    ln -sfn "releases/$VERSION" "$DB_ROOT/current"

    if ! validate_database "$FINAL_RELEASE"; then
        echo 'ERROR: The installed managed vConTACT3 database failed final validation.' >&2
        exit 1
    fi
    emit_database "$FINAL_RELEASE" 'downloaded-latest-and-installed'
    '''
}
