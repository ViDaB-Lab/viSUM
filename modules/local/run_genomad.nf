process RUN_GENOMAD {

    tag "${prefix}"

    conda "bioconda::genomad=${params.genomad_version}"

    cpus { Math.min(params.genomad_cpus as int, params.max_cpus as int) }
    memory params.genomad_memory
    time params.genomad_time

    publishDir { "${params.outdir}/${prefix}_results/genomad" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)
    path genomad_db

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.genomad"),
          path("${prefix}_virus_summary.tsv"),
          path("${prefix}_virus.fna"),
          path("${prefix}_virus_genes.tsv"),
          path("${prefix}_virus_proteins.faa"),
          path("${prefix}_plasmid_summary.tsv"),
          path("${prefix}.genomad_run_metadata.tsv"),
          emit: results

    script:
    def cleanupFlag = params.genomad_cleanup ? '--cleanup' : ''
    def calibrationEnabled = params.genomad_score_calibration instanceof Boolean
        ? params.genomad_score_calibration
        : params.genomad_score_calibration.toString().trim().equalsIgnoreCase('true')
    def calibrationFlag = calibrationEnabled ? '--enable-score-calibration' : ''
    """
    set -euo pipefail

    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export NUMEXPR_NUM_THREADS="${task.cpus}"

    ln -s "${normalized_fasta}" "${prefix}.fasta"

    genomad end-to-end \
        ${cleanupFlag} \
        ${calibrationFlag} \
        --threads ${task.cpus} \
        --splits ${params.genomad_splits} \
        "${prefix}.fasta" \
        "${prefix}.genomad" \
        "${genomad_db}"

    SUMMARY_DIR="${prefix}.genomad/${prefix}_summary"

    if [[ ! -d "\$SUMMARY_DIR" ]]; then
        echo "ERROR: geNomad completed without the expected summary directory:" >&2
        echo "       \$SUMMARY_DIR" >&2
        exit 1
    fi

    for expected_file in \
        "${prefix}_virus_summary.tsv" \
        "${prefix}_virus.fna" \
        "${prefix}_virus_genes.tsv" \
        "${prefix}_virus_proteins.faa" \
        "${prefix}_plasmid_summary.tsv"
    do
        if [[ ! -f "\$SUMMARY_DIR/\$expected_file" ]]; then
            echo "ERROR: Missing expected geNomad output: \$SUMMARY_DIR/\$expected_file" >&2
            exit 1
        fi
        cp "\$SUMMARY_DIR/\$expected_file" "\$expected_file"
    done

    VIRUS_CALL_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' \
        "${prefix}_virus_summary.tsv")
    PLASMID_CALL_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' \
        "${prefix}_plasmid_summary.tsv")
    INPUT_SEQUENCE_COUNT=\$(awk '/^>/{ count++ } END { print count + 0 }' \
        "${prefix}.fasta")
    SCORE_CALIBRATION_REQUESTED='${calibrationEnabled}'
    SCORE_CALIBRATION_APPLIED='false'

    # geNomad may decide whether calibration is usable for a particular input.
    # Record what it actually produced instead of inferring this from input size.
    TOTAL_CALL_COUNT=\$((VIRUS_CALL_COUNT + PLASMID_CALL_COUNT))
    FDR_CALL_COUNT=\$(awk -F '\t' '
        FNR == 1 {
            fdr_column = 0
            for (column = 1; column <= NF; column++) {
                if (\$column == "fdr") {
                    fdr_column = column
                    break
                }
            }
            next
        }
        fdr_column > 0 && tolower(\$fdr_column) !~ /^(|na|nan|none)$/ { count++ }
        END { print count + 0 }
    ' "${prefix}_virus_summary.tsv" "${prefix}_plasmid_summary.tsv")

    if [[ "\$FDR_CALL_COUNT" -gt 0 && \
          "\$FDR_CALL_COUNT" -ne "\$TOTAL_CALL_COUNT" ]]; then
        echo "ERROR: geNomad produced FDR values for only \$FDR_CALL_COUNT of \$TOTAL_CALL_COUNT calls." >&2
        exit 1
    fi
    if [[ "\$FDR_CALL_COUNT" -gt 0 || \
          -d "${prefix}.genomad/${prefix}_score_calibration" ]]; then
        SCORE_CALIBRATION_APPLIED='true'
    fi
    if [[ "\$VIRUS_CALL_COUNT" -eq 0 ]]; then
        RUN_STATUS='completed_no_viruses_detected'
    else
        RUN_STATUS='completed_with_virus_calls'
    fi

    printf 'sample_id\tinput_type\tinput_fasta\tinput_sequence_count\tgenomad_version\tthreads\tsplits\tcleanup\tscore_calibration_requested\tscore_calibration_applied\trun_status\tvirus_call_count\tplasmid_call_count\n' \
        > "${prefix}.genomad_run_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" \
        "${type}" \
        "${normalized_fasta.name}" \
        "\$INPUT_SEQUENCE_COUNT" \
        "\$(genomad --version 2>&1 | head -n 1)" \
        "${task.cpus}" \
        "${params.genomad_splits}" \
        "${params.genomad_cleanup}" \
        "\$SCORE_CALIBRATION_REQUESTED" \
        "\$SCORE_CALIBRATION_APPLIED" \
        "\$RUN_STATUS" \
        "\$VIRUS_CALL_COUNT" \
        "\$PLASMID_CALL_COUNT" \
        >> "${prefix}.genomad_run_metadata.tsv"
    """
}
