process STANDARDIZE_TESORTER {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/tesorter" },
        mode: 'copy',
        pattern: '*.tesorter_*.tsv'

    input:
    tuple val(prefix),
          val(type),
          path(region_map),
          path(classifications),
          path(domains),
          path(domain_gff),
          path(log),
          path(run_metadata)

    output:
    tuple val(prefix),
          val('tesorter'),
          path("${prefix}.tesorter_evidence.tsv"),
          emit: evidence
    tuple val(prefix),
          path("${prefix}.tesorter_standardization_audit.tsv"),
          emit: audit

    script:
    """
    python3 "${projectDir}/bin/standardize_tesorter.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --region-map "${region_map}" \
        --classifications "${classifications}" \
        --domains "${domains}" \
        --run-metadata "${run_metadata}" \
        --output-evidence "${prefix}.tesorter_evidence.tsv" \
        --output-audit "${prefix}.tesorter_standardization_audit.tsv"
    """
}
