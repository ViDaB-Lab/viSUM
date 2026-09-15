process STANDARDIZE_DEEPMICROCLASS2 {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/deepmicroclass2" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(score_table),
          path(run_metadata),
          path(header_map)

    output:
    tuple val(prefix),
          val('deepmicroclass2'),
          path("${prefix}.deepmicroclass2_evidence.tsv"),
          emit: evidence

    script:
    """
    python3 "${projectDir}/bin/standardize_deepmicroclass2.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --header-map "${header_map}" \
        --score-table "${score_table}" \
        --run-metadata "${run_metadata}" \
        --output "${prefix}.deepmicroclass2_evidence.tsv"
    """
}
