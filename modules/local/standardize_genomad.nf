process STANDARDIZE_GENOMAD {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/genomad" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(virus_summary),
          path(virus_genes),
          path(plasmid_summary),
          path(run_metadata),
          path(header_map)

    output:
    tuple val(prefix),
          val('genomad'),
          path("${prefix}.genomad_evidence.tsv"),
          emit: evidence

    script:
    """
    python3 "${projectDir}/bin/standardize_genomad.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --header-map "${header_map}" \
        --virus-summary "${virus_summary}" \
        --virus-genes "${virus_genes}" \
        --plasmid-summary "${plasmid_summary}" \
        --run-metadata "${run_metadata}" \
        --output "${prefix}.genomad_evidence.tsv"
    """
}
