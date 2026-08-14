process REFINE_PROVIRAL_REGIONS {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/refinement" }, mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(candidate_fasta),
          path(discovery_audit),
          path(discovery_summary),
          val(evidence_file_count),
          path(evidence_files),
          val(allow_ct3_only_refinement)

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.refined_candidates.fasta"),
          path("${prefix}.provirus_region_map.tsv"),
          path("${prefix}.provirus_boundary_audit.tsv"),
          path("${prefix}.provirus_refinement_summary.tsv"),
          emit: refined

    script:
    def evidenceList = evidence_files instanceof List ? evidence_files : [evidence_files]
    def evidenceArguments = evidence_file_count > 0
        ? evidenceList.collect { evidence -> "'${evidence}'" }.join(' ')
        : ''
    def allowCt3OnlyArgument = allow_ct3_only_refinement
        ? '--allow-ct3-only-refinement'
        : ''
    """
    python3 "${projectDir}/bin/refine_proviral_regions.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --candidate-fasta "${candidate_fasta}" \
        --evidence ${evidenceArguments} \
        ${allowCt3OnlyArgument} \
        --output-fasta "${prefix}.refined_candidates.fasta" \
        --output-map "${prefix}.provirus_region_map.tsv" \
        --output-audit "${prefix}.provirus_boundary_audit.tsv" \
        --output-summary "${prefix}.provirus_refinement_summary.tsv"
    """
}
