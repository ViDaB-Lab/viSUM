process RUN_VITAP {

    tag "${prefix}"

    conda "${projectDir}/envs/vitap.yml"

    cpus { Math.min(params.vitap_cpus as int, params.max_cpus as int) }
    memory params.vitap_memory
    time params.vitap_time

    // Publish the stable review files. The full upstream directory remains a
    // task output for provenance, but is not copied because its alignments can
    // be very large and are reproducible from the refined FASTA and database.
    publishDir { "${params.outdir}/${prefix}_results/vitap" },
        mode: 'copy',
        pattern: '*.vitap_*.*'

    input:
    tuple val(prefix),
          val(type),
          path(refined_fasta),
          path(provirus_region_map),
          path(provirus_boundary_audit),
          path(provirus_refinement_summary)
    tuple path(vitap_database), path(vitap_database_metadata)
    val include_low_confidence

    output:
    tuple val(prefix),
          val(type),
          path(provirus_region_map),
          path("${prefix}.vitap_best_determined_lineages.tsv"),
          path("${prefix}.vitap_all_lineages.tsv"),
          path("${prefix}.vitap_uniref90_fallback.tsv"),
          path("${prefix}.vitap.log"),
          path("${prefix}.vitap_run_metadata.tsv"),
          emit: results
    path("${prefix}.vitap_raw"), emit: raw_output

    script:
    def lowConfidenceArgument = include_low_confidence ? '--low_conf' : ''
    """
    set -euo pipefail

    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export NUMEXPR_NUM_THREADS="${task.cpus}"
    export POLARS_MAX_THREADS="${task.cpus}"

    RAW_DIRECTORY="${prefix}.vitap_raw"
    BEST_FILE="${prefix}.vitap_best_determined_lineages.tsv"
    ALL_FILE="${prefix}.vitap_all_lineages.tsv"
    FALLBACK_FILE="${prefix}.vitap_uniref90_fallback.tsv"
    LOG_FILE="${prefix}.vitap.log"
    METADATA_FILE="${prefix}.vitap_run_metadata.tsv"
    mkdir -p "\$RAW_DIRECTORY"

    INPUT_COUNT=\$(awk '/^>/ { count++ } END { print count + 0 }' "${refined_fasta}")

    if [[ "\$INPUT_COUNT" -eq 0 ]]; then
        printf 'Genome_ID\tlineage\tlineage_score/participation_index\tConfidence_level\n' \
            > "\$RAW_DIRECTORY/best_determined_lineages.tsv"
        printf 'Genome_ID\tlineage\tlineage_score/participation_index\n' \
            > "\$RAW_DIRECTORY/all_lineages.tsv"
        printf 'genome_id\ttaxa_name\tparticipation_index\ttaxon_level\n' \
            > "\$RAW_DIRECTORY/target_uniref90_taxa_fallback.out"
        printf 'VITAP skipped: provirus refinement produced no candidate sequences.\n' \
            > "\$LOG_FILE"
        RUN_STATUS='skipped_no_refined_candidates'
    else
        if ! VITAP assignment \
            -i "${refined_fasta}" \
            -d "${vitap_database}" \
            -p "${task.cpus}" \
            -o "\$RAW_DIRECTORY" \
            ${lowConfidenceArgument} \
            > "\$LOG_FILE" 2>&1; then
            echo "ERROR: VITAP assignment failed for sample '${prefix}'." >&2
            echo "Last 50 lines of \$LOG_FILE:" >&2
            tail -n 50 "\$LOG_FILE" >&2 || true
            exit 1
        fi
        RUN_STATUS='completed'
    fi

    for required_file in best_determined_lineages.tsv all_lineages.tsv; do
        if [[ ! -s "\$RAW_DIRECTORY/\$required_file" ]]; then
            echo "ERROR: VITAP completed without required output: \$required_file" >&2
            exit 1
        fi
    done

    if [[ ! -f "\$RAW_DIRECTORY/target_uniref90_taxa_fallback.out" ]]; then
        printf 'genome_id\ttaxa_name\tparticipation_index\ttaxon_level\n' \
            > "\$RAW_DIRECTORY/target_uniref90_taxa_fallback.out"
    fi

    # VITAP writes these tables through Python's csv module, which may use
    # CRLF even on Linux. Remove only the trailing carriage return before
    # comparing the schema; preserve the original upstream files unchanged.
    BEST_HEADER=\$(head -n 1 "\$RAW_DIRECTORY/best_determined_lineages.tsv" | tr -d '\r')
    ALL_HEADER=\$(head -n 1 "\$RAW_DIRECTORY/all_lineages.tsv" | tr -d '\r')
    if [[ "\$BEST_HEADER" != \$'Genome_ID\tlineage\tlineage_score/participation_index\tConfidence_level' ]]; then
        echo 'ERROR: VITAP best-lineage output has an unexpected schema.' >&2
        printf 'Observed header: %s\n' "\$BEST_HEADER" >&2
        exit 1
    fi
    if [[ "\$ALL_HEADER" != \$'Genome_ID\tlineage\tlineage_score/participation_index' ]]; then
        echo 'ERROR: VITAP all-lineages output has an unexpected schema.' >&2
        printf 'Observed header: %s\n' "\$ALL_HEADER" >&2
        exit 1
    fi

    cp "\$RAW_DIRECTORY/best_determined_lineages.tsv" "\$BEST_FILE"
    cp "\$RAW_DIRECTORY/all_lineages.tsv" "\$ALL_FILE"
    cp "\$RAW_DIRECTORY/target_uniref90_taxa_fallback.out" "\$FALLBACK_FILE"

    BEST_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' "\$BEST_FILE")
    ALL_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' "\$ALL_FILE")
    FALLBACK_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' "\$FALLBACK_FILE")
    VITAP_VERSION=\$(python -c "from importlib.metadata import version; print(version('VITAP'))")
    DATABASE_RELEASE=\$(awk -F '\t' 'NR == 2 { print \$5 }' "${vitap_database_metadata}")
    RESOLVED_DATABASE_PATH=\$(readlink -f "${vitap_database}")
    if [[ -z "\$RESOLVED_DATABASE_PATH" || ! -d "\$RESOLVED_DATABASE_PATH" ]]; then
        echo "ERROR: Could not resolve VITAP database path from '${vitap_database}'." >&2
        exit 1
    fi

    printf 'sample_id\tinput_type\trefined_sequence_count\tbest_lineage_row_count\tall_lineage_row_count\tuniref90_fallback_row_count\tinclude_low_confidence\tthreads\tvitap_version\tdatabase_path\tdatabase_release\trun_status\tinput_fasta\n' \
        > "\$METADATA_FILE"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$INPUT_COUNT" "\$BEST_COUNT" "\$ALL_COUNT" \
        "\$FALLBACK_COUNT" "${include_low_confidence}" "${task.cpus}" "\$VITAP_VERSION" \
        "\$RESOLVED_DATABASE_PATH" "\$DATABASE_RELEASE" "\$RUN_STATUS" \
        "${refined_fasta.name}" >> "\$METADATA_FILE"

    echo "VITAP sample=${prefix} type=${type} refined=\$INPUT_COUNT best=\$BEST_COUNT all=\$ALL_COUNT fallback=\$FALLBACK_COUNT"
    """
}
