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

    RAW_DIRECTORY="${prefix}.deep6_raw"
    RAW_SCORE_FILE="\$RAW_DIRECTORY/${normalized_fasta.name}_predict_${params.deep6_minlen}bp_deep6.txt"
    FINAL_SCORE_FILE="${prefix}.deep6_scores.tsv"
    LOG_FILE="${prefix}.deep6.log"
    METADATA_FILE="${prefix}.deep6_run_metadata.tsv"

    mkdir -p "\$RAW_DIRECTORY"

    export CUDA_VISIBLE_DEVICES='-1'
    export OMP_NUM_THREADS="${task.cpus}"
    export TF_NUM_INTRAOP_THREADS="${task.cpus}"
    export TF_NUM_INTEROP_THREADS='1'

    if ! python "${deep6_installation}/Master/deep6.py" \
        -i "${normalized_fasta}" \
        -l ${params.deep6_minlen} \
        -m "${deep6_models}" \
        -o "\$RAW_DIRECTORY" \
        > "\$LOG_FILE" 2>&1; then
        echo "ERROR: Deep6 prediction failed for sample '${prefix}'." >&2
        echo "Last 50 lines of \$LOG_FILE:" >&2
        tail -n 50 "\$LOG_FILE" >&2 || true
        exit 1
    fi

    if [[ ! -s "\$RAW_SCORE_FILE" ]]; then
        echo "ERROR: Deep6 completed without the expected score file:" >&2
        echo "       \$RAW_SCORE_FILE" >&2
        exit 1
    fi

    EXPECTED_HEADER=\$'name\tlength\tduplo\teuk\tmono\tpro\tribo\tvari'
    ACTUAL_HEADER=\$(head -n 1 "\$RAW_SCORE_FILE" | tr -d '\r')
    if [[ "\$ACTUAL_HEADER" != "\$EXPECTED_HEADER" ]]; then
        echo "ERROR: Deep6 produced an unexpected score-table header:" >&2
        echo "       \$ACTUAL_HEADER" >&2
        exit 1
    fi

    cp "\$RAW_SCORE_FILE" "\$FINAL_SCORE_FILE"

    PREDICTION_COUNT=\$(
        awk 'NR > 1 && NF { count++ } END { print count + 0 }' \
            "\$FINAL_SCORE_FILE"
    )
    DEEP6_REVISION=\$(
        awk -F '\t' 'NR == 2 { print \$9 }' "${deep6_database_metadata}"
    )
    [[ -n "\$DEEP6_REVISION" ]] || DEEP6_REVISION='unknown'

    printf 'sample_id\tinput_type\tminimum_length\tprediction_count\tdeep6_version\tdeep6_revision\traw_score_file\n' \
        > "\$METADATA_FILE"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "${params.deep6_minlen}" "\$PREDICTION_COUNT" \
        '1' "\$DEEP6_REVISION" "\$FINAL_SCORE_FILE" \
        >> "\$METADATA_FILE"

    echo "Deep6 sample=${prefix} type=${type} predictions=\$PREDICTION_COUNT scores=\$FINAL_SCORE_FILE"
    """
}
