process STANDARDIZE_VIRSORTER2 {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/virsorter2" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(score_table),
          path(boundary_table),
          path(viral_fasta),
          path(run_metadata),
          path(normalized_fasta),
          path(header_map)

    output:
    tuple val(prefix),
          val('virsorter2'),
          path("${prefix}.virsorter2_evidence.tsv"),
          emit: evidence
    tuple val(prefix),
          path("${prefix}.virsorter2_boundary_audit.tsv"),
          emit: boundary_audit

    script:
    """
    # Adapter policy v1.2: preserve unscored lt2gene calls as review-only evidence.
    python3 "${projectDir}/bin/standardize_virsorter2.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --header-map "${header_map}" \
        --score-table "${score_table}" \
        --boundary-table "${boundary_table}" \
        --viral-fasta "${viral_fasta}" \
        --run-metadata "${run_metadata}" \
        --normalized-fasta "${normalized_fasta}" \
        --boundary-audit "${prefix}.virsorter2_boundary_audit.tsv" \
        --output "${prefix}.virsorter2_evidence.tsv"
    """
}
