process RUN_DEEPMICROCLASS2 {

    tag "${prefix}"

    conda "${projectDir}/envs/deepmicroclass2.yml"

    cpus { Math.min(params.deepmicroclass2_cpus as int, params.max_cpus as int) }
    memory params.deepmicroclass2_memory
    time params.deepmicroclass2_time

    publishDir { "${params.outdir}/${prefix}_results/deepmicroclass2" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)
    tuple path(deepmicroclass2_installation),
          path(deepmicroclass2_installation_metadata)

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.deepmicroclass2_scores.tsv"),
          path("${prefix}.deepmicroclass2.log"),
          path("${prefix}.deepmicroclass2_run_metadata.tsv"),
          emit: results

    script:
    """
    set -euo pipefail

    if [[ "${type}" != 'dna' ]]; then
        echo "ERROR: DeepMicroClass2 received non-DNA sample '${prefix}' with type '${type}'." >&2
        exit 1
    fi

    RAW_DIRECTORY="${prefix}.deepmicroclass2_raw"
    RAW_SCORE_FILE="\$RAW_DIRECTORY/classification.tsv"
    FINAL_SCORE_FILE="${prefix}.deepmicroclass2_scores.tsv"
    LOG_FILE="${prefix}.deepmicroclass2.log"
    METADATA_FILE="${prefix}.deepmicroclass2_run_metadata.tsv"
    CLASS_THRESHOLDS='arc=0.625000,bac=0.400000,chlor=0.390625,euk=0.435547,eukvir=0.951172,mit=0.261719,pls=0.832031,prokvir=0.997925'

    case "${params.deepmicroclass2_model}" in
        8class)
            MINIMUM_LENGTH=500
            ;;
        high_precision|300bp)
            MINIMUM_LENGTH=300
            ;;
        *)
            echo "ERROR: Unsupported DeepMicroClass2 model mode '${params.deepmicroclass2_model}'." >&2
            exit 1
            ;;
    esac

    mkdir -p "\$RAW_DIRECTORY"

    # The initial integration uses CPU inference for portability. PyTorch uses
    # these limits so one task cannot silently consume every server core.
    export CUDA_VISIBLE_DEVICES=''
    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export NUMEXPR_NUM_THREADS="${task.cpus}"

    if ! python "${deepmicroclass2_installation}/predict.py" \
        --contig "${normalized_fasta}" \
        --model "${params.deepmicroclass2_model}" \
        --out_dir "\$RAW_DIRECTORY" \
        > "\$LOG_FILE" 2>&1; then
        echo "ERROR: DeepMicroClass2 prediction failed for sample '${prefix}'." >&2
        echo "Last 50 lines of \$LOG_FILE:" >&2
        tail -n 50 "\$LOG_FILE" >&2 || true
        exit 1
    fi

    if [[ ! -s "\$RAW_SCORE_FILE" ]]; then
        echo "ERROR: DeepMicroClass2 completed without classification.tsv." >&2
        exit 1
    fi

    EXPECTED_HEADER=\$'contig\tlabel\tconfidence\tarc\tbac\tchlor\teuk\teukvir\tmit\tpls\tprokvir'
    ACTUAL_HEADER=\$(head -n 1 "\$RAW_SCORE_FILE" | tr -d '\r')
    if [[ "\$ACTUAL_HEADER" != "\$EXPECTED_HEADER" ]]; then
        echo "ERROR: DeepMicroClass2 produced an unexpected score-table header:" >&2
        echo "       \$ACTUAL_HEADER" >&2
        exit 1
    fi

    cp "\$RAW_SCORE_FILE" "\$FINAL_SCORE_FILE"

    INPUT_COUNT=\$(awk 'NR > 1 && NF { count++ } END { print count + 0 }' "${header_map}")
    ELIGIBLE_COUNT=\$(awk -F '\t' -v min_length="\$MINIMUM_LENGTH" \
        'NR > 1 && \$7 >= min_length { count++ } END { print count + 0 }' \
        "${header_map}")
    PREDICTION_COUNT=\$(awk 'NR > 1 && NF { count++ } END { print count + 0 }' "\$FINAL_SCORE_FILE")

    if [[ "\$PREDICTION_COUNT" -ne "\$ELIGIBLE_COUNT" ]]; then
        echo "ERROR: DeepMicroClass2 prediction coverage mismatch for sample '${prefix}':" >&2
        echo "       eligible sequences: \$ELIGIBLE_COUNT" >&2
        echo "       prediction rows:    \$PREDICTION_COUNT" >&2
        exit 1
    fi

    INSTALLED_REVISION=\$(
        awk -F '\t' 'NR == 2 { print \$8 }' "${deepmicroclass2_installation_metadata}"
    )
    [[ -n "\$INSTALLED_REVISION" ]] || INSTALLED_REVISION='unknown'

    RUN_STATUS='completed'
    if [[ "\$PREDICTION_COUNT" -eq 0 ]]; then
        RUN_STATUS='completed_no_eligible_sequences'
    fi

    printf 'sample_id\tinput_type\tmodel_mode\tminimum_length\tinput_sequence_count\teligible_sequence_count\tprediction_count\tclass_thresholds\tdeepmicroclass2_revision\tdevice\tthreads\trun_status\tscore_file\n' \
        > "\$METADATA_FILE"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "${params.deepmicroclass2_model}" "\$MINIMUM_LENGTH" \
        "\$INPUT_COUNT" "\$ELIGIBLE_COUNT" "\$PREDICTION_COUNT" \
        "\$CLASS_THRESHOLDS" "\$INSTALLED_REVISION" 'cpu' "${task.cpus}" "\$RUN_STATUS" \
        "\$FINAL_SCORE_FILE" \
        >> "\$METADATA_FILE"

    echo "DeepMicroClass2 sample=${prefix} type=${type} predictions=\$PREDICTION_COUNT scores=\$FINAL_SCORE_FILE"
    """
}
