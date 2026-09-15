process STANDARDIZE_DEEP6 {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/deep6" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(score_table),
          path(run_metadata),
          path(header_map)

    output:
    tuple val(prefix),
          val('deep6'),
          path("${prefix}.deep6_evidence.tsv"),
          emit: evidence

    script:
    """
    python3 "${projectDir}/bin/standardize_deep6.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --header-map "${header_map}" \
        --score-table "${score_table}" \
        --run-metadata "${run_metadata}" \
        --minimum-score "${params.deep6_min_score}" \
        --median-multiplier "${params.deep6_median_multiplier}" \
        --output "${prefix}.deep6_evidence.tsv"
    """
}
