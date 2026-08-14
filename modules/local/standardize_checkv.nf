process STANDARDIZE_CHECKV {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/checkv" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(quality_summary),
          path(completeness),
          path(contamination),
          path(complete_genomes),
          path(run_metadata),
          path(candidate_fasta)

    output:
    tuple val(prefix),
          val('checkv'),
          path("${prefix}.checkv_evidence.tsv"),
          emit: evidence

    script:
    """
    python3 "${projectDir}/bin/standardize_checkv.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --candidate-fasta "${candidate_fasta}" \
        --quality-summary "${quality_summary}" \
        --completeness "${completeness}" \
        --contamination "${contamination}" \
        --complete-genomes "${complete_genomes}" \
        --run-metadata "${run_metadata}" \
        --output "${prefix}.checkv_evidence.tsv"
    """
}
