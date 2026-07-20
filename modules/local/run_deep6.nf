process RUN_DEEP6 {

    tag "${prefix}"

    conda "${projectDir}/envs/deep6.yml"

    cpus 1
    memory params.deep6_memory
    time params.deep6_time

    publishDir { "${params.outdir}/${prefix}_results/deep6" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)
    tuple path(deep6_installation),
          path(deep6_models),
          path(deep6_database_metadata)

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.deep6_scores.tsv"),
          path("${prefix}.deep6.log"),
          path("${prefix}.deep6_run_metadata.tsv"),
          emit: results

    script:
    """
    set -euo pipefail

    if [[ "${type}" != 'rna' ]]; then
        echo "ERROR: Deep6 received non-RNA sample '${prefix}' with type '${type}'." >&2
        exit 1
    fi

    DEEP6_INPUT="${prefix}.deep6_input.fasta"
    FINAL_SCORE_FILE="${prefix}.deep6_scores.tsv"
    LOG_FILE="${prefix}.deep6.log"
    METADATA_FILE="${prefix}.deep6_run_metadata.tsv"

    # Deep6 explicitly encodes A, C, G, and T. Convert uracil to its DNA-alphabet
    # equivalent for RNA FASTA files while preserving identifiers and lengths.
    URACIL_COUNT=\$(
        awk '!/^>/ { sequence = \$0; count += gsub(/[Uu]/, "", sequence) }
             END { print count + 0 }' "${normalized_fasta}"
    )
    awk '/^>/ { print; next }
         { gsub(/[Uu]/, "T"); print }' \
        "${normalized_fasta}" > "\$DEEP6_INPUT"

    export CUDA_VISIBLE_DEVICES='-1'
    export OMP_NUM_THREADS="${task.cpus}"
    export TF_NUM_INTRAOP_THREADS="${task.cpus}"
    export TF_NUM_INTEROP_THREADS='1'

    if ! python "${projectDir}/bin/run_deep6.py" \
        --input "\$DEEP6_INPUT" \
        --minimum-length ${params.deep6_minlen} \
        --deep6-installation "${deep6_installation}" \
        --models "${deep6_models}" \
        --output "\$FINAL_SCORE_FILE" \
        > "\$LOG_FILE" 2>&1; then
        echo "ERROR: Deep6 prediction failed for sample '${prefix}'." >&2
        echo "Last 50 lines of \$LOG_FILE:" >&2
        tail -n 50 "\$LOG_FILE" >&2 || true
        exit 1
    fi

    if [[ ! -s "\$FINAL_SCORE_FILE" ]]; then
        echo "ERROR: Deep6 completed without the expected score file:" >&2
        echo "       \$FINAL_SCORE_FILE" >&2
        exit 1
    fi

    EXPECTED_HEADER=\$'name\tlength\tduplo\teuk\tmono\tpro\tribo\tvari'
    ACTUAL_HEADER=\$(head -n 1 "\$FINAL_SCORE_FILE" | tr -d '\r')
    if [[ "\$ACTUAL_HEADER" != "\$EXPECTED_HEADER" ]]; then
        echo "ERROR: Deep6 produced an unexpected score-table header:" >&2
        echo "       \$ACTUAL_HEADER" >&2
        exit 1
    fi

    PREDICTION_COUNT=\$(
        awk 'NR > 1 && NF { count++ } END { print count + 0 }' \
            "\$FINAL_SCORE_FILE"
    )
    DEEP6_REVISION=\$(
        awk -F '\t' 'NR == 2 { print \$9 }' "${deep6_database_metadata}"
    )
    [[ -n "\$DEEP6_REVISION" ]] || DEEP6_REVISION='unknown'

    printf 'sample_id\tinput_type\tminimum_length\tprediction_count\turacil_to_thymine_conversions\tdeep6_version\tdeep6_revision\traw_score_file\n' \
        > "\$METADATA_FILE"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "${params.deep6_minlen}" "\$PREDICTION_COUNT" \
        "\$URACIL_COUNT" '1' "\$DEEP6_REVISION" "\$FINAL_SCORE_FILE" \
        >> "\$METADATA_FILE"

    echo "Deep6 sample=${prefix} type=${type} predictions=\$PREDICTION_COUNT scores=\$FINAL_SCORE_FILE"
    """
}
