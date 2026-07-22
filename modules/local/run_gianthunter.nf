process RUN_GIANTHUNTER {

    tag "${prefix}"

    conda "${projectDir}/envs/gianthunter.yml"

    cpus params.threads
    memory params.gianthunter_memory
    time params.gianthunter_time

    publishDir { "${params.outdir}/${prefix}_results/gianthunter" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)
    tuple path(gianthunter_database), path(gianthunter_database_metadata)

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.gianthunter_prediction.tsv"),
          path("${prefix}.gianthunter_gene_annotations.tsv"),
          path("${prefix}.gianthunter.log"),
          path("${prefix}.gianthunter_run_metadata.tsv"),
          emit: results

    script:
    """
    set -euo pipefail

    if [[ "${type}" != 'dna' ]]; then
        echo "ERROR: GiantHunter received non-DNA sample '${prefix}' with type '${type}'." >&2
        exit 1
    fi

    RAW_DIRECTORY="${prefix}.gianthunter_raw"
    RAW_PREDICTION="\$RAW_DIRECTORY/final_prediction/gianthunter_prediction.tsv"
    RAW_ANNOTATIONS="\$RAW_DIRECTORY/final_prediction/supplementary/gene_annotation.tsv"
    FINAL_PREDICTION="${prefix}.gianthunter_prediction.tsv"
    FINAL_ANNOTATIONS="${prefix}.gianthunter_gene_annotations.tsv"
    LOG_FILE="${prefix}.gianthunter.log"
    METADATA_FILE="${prefix}.gianthunter_run_metadata.tsv"
    EXPECTED_HEADER='Accession\tLength\tGiantVirus\tPotentialLineage\tScore'
    NO_HIT_HEADER='Accession\tLength\tPotentialLineage\tScore\tGenus\tGenusCluster'
    NO_REFERENCE_HITS='false'

    INPUT_COUNT=\$(awk 'NR > 1 && NF { count++ } END { print count + 0 }' "${header_map}")
    ELIGIBLE_COUNT=\$(awk -F '\t' -v minimum="${params.gianthunter_min_length}" \
        'NR > 1 && (\$7 + 0) >= minimum { count++ } END { print count + 0 }' \
        "${header_map}")

    if [[ "\$ELIGIBLE_COUNT" -eq 0 ]]; then
        printf '%b\n' "\$EXPECTED_HEADER" > "\$FINAL_PREDICTION"
        : > "\$FINAL_ANNOTATIONS"
        printf 'GiantHunter skipped: no input sequences met the minimum length of %s nt.\n' \
            "${params.gianthunter_min_length}" > "\$LOG_FILE"
    else
        if ! gianthunter \
            --contigs "${normalized_fasta}" \
            --len "${params.gianthunter_min_length}" \
            --threads "${task.cpus}" \
            --dbdir "${gianthunter_database}" \
            --outpth "\$RAW_DIRECTORY" \
            --midfolder midfolder \
            --reject "${params.gianthunter_reject}" \
            --query_cover "${params.gianthunter_query_cover}" \
            > "\$LOG_FILE" 2>&1; then
            echo "ERROR: GiantHunter prediction failed for sample '${prefix}'." >&2
            echo "Last 50 lines of \$LOG_FILE:" >&2
            tail -n 50 "\$LOG_FILE" >&2 || true
            exit 1
        fi

        if [[ ! -s "\$RAW_PREDICTION" ]]; then
            echo "ERROR: GiantHunter completed without gianthunter_prediction.tsv." >&2
            exit 1
        fi

        ACTUAL_HEADER=\$(head -n 1 "\$RAW_PREDICTION" | tr -d '\r')
        if [[ "\$ACTUAL_HEADER" == "\$(printf '%b' "\$NO_HIT_HEADER")" ]]; then
            # GiantHunter intentionally emits this alternate table and exits
            # successfully when no predicted protein hits RefVirus.dmnd.
            NO_REFERENCE_HITS='true'
        elif [[ "\$ACTUAL_HEADER" != "\$(printf '%b' "\$EXPECTED_HEADER")" ]]; then
            echo "ERROR: GiantHunter produced an unexpected prediction-table header:" >&2
            echo "       \$ACTUAL_HEADER" >&2
            exit 1
        fi

        cp "\$RAW_PREDICTION" "\$FINAL_PREDICTION"
    fi

    GIANT_VIRUS_CALL_COUNT=\$(awk -F '\t' \
        'NR > 1 && \$3 == "GiantVirus" { count++ } END { print count + 0 }' \
        "\$FINAL_PREDICTION")

    if [[ "\$ELIGIBLE_COUNT" -gt 0 ]]; then
        if [[ -f "\$RAW_ANNOTATIONS" ]]; then
            cp "\$RAW_ANNOTATIONS" "\$FINAL_ANNOTATIONS"
        else
            # Gene annotations are supplemental context. Their absence must
            # not invalidate an otherwise valid GiantHunter prediction table.
            : > "\$FINAL_ANNOTATIONS"
            echo "WARNING: GiantHunter produced no gene-annotation table; n_genes will be blank." \
                >> "\$LOG_FILE"
        fi
    fi

    GIANTHUNTER_VERSION=\$(
        python -c "from importlib.metadata import version; print(version('gianthunter'))" \
            2>/dev/null || true
    )
    [[ -n "\$GIANTHUNTER_VERSION" ]] || GIANTHUNTER_VERSION='unknown'

    if [[ "\$ELIGIBLE_COUNT" -eq 0 ]]; then
        RUN_STATUS='skipped_no_sequences_meeting_minimum_length'
    elif [[ "\$NO_REFERENCE_HITS" == 'true' ]]; then
        RUN_STATUS='completed_no_reference_protein_hits'
    elif [[ "\$GIANT_VIRUS_CALL_COUNT" -eq 0 ]]; then
        RUN_STATUS='completed_no_giant_virus_calls'
    else
        RUN_STATUS='completed_with_giant_virus_calls'
    fi

    printf 'sample_id\tinput_type\tinput_sequence_count\teligible_sequence_count\tgiant_virus_call_count\tgianthunter_version\tdatabase_path\tminimum_length\treject_threshold\tquery_cover\tthreads\trun_status\tprediction_file\n' \
        > "\$METADATA_FILE"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$INPUT_COUNT" "\$ELIGIBLE_COUNT" \
        "\$GIANT_VIRUS_CALL_COUNT" "\$GIANTHUNTER_VERSION" \
        "${gianthunter_database}" "${params.gianthunter_min_length}" \
        "${params.gianthunter_reject}" "${params.gianthunter_query_cover}" \
        "${task.cpus}" "\$RUN_STATUS" "\$FINAL_PREDICTION" \
        >> "\$METADATA_FILE"

    echo "GIANTHUNTER sample=${prefix} type=${type} calls=\$GIANT_VIRUS_CALL_COUNT prediction=\$FINAL_PREDICTION"
    """
}
