process STANDARDIZE_CENOTETAKER3 {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/cenotetaker3" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(virus_summary),
          path(virus_fasta),
          path(prune_summary),
          path(gene_annotations),
          path(run_metadata),
          path(header_map)

    output:
    tuple val(prefix),
          val('cenotetaker3'),
          path("${prefix}.cenotetaker3_evidence.tsv"),
          emit: evidence

    script:
    """
    python3 "${projectDir}/bin/standardize_cenotetaker3.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --header-map "${header_map}" \
        --virus-summary "${virus_summary}" \
        --virus-fasta "${virus_fasta}" \
        --prune-summary "${prune_summary}" \
        --gene-annotations "${gene_annotations}" \
        --run-metadata "${run_metadata}" \
        --output "${prefix}.cenotetaker3_evidence.tsv"
    """
}
