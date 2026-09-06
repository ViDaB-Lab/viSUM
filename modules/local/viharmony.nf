process VIHARMONY {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'
    cpus 1
    memory params.harmonizer_memory
    time params.harmonizer_time

    publishDir { "${params.outdir}/${prefix}_results/viharmony" }, mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(normalized_fasta),
          path(header_map),
          path(discovery_gate),
          path(refined_fasta),
          path(region_map),
          val(evidence_file_count),
          path(evidence_files),
          val(group_file_count),
          path(vcontact3_groups),
          path(ictv_msl),
          val(audit_mode),
          val(vcontact3_min_taxonomy_length)

    output:
    tuple val(prefix),
          path("${prefix}.final.normalized.fasta"),
          path("${prefix}.final.original_ids.fasta"),
          path("${prefix}.final_metadata.tsv"),
          path("${prefix}.review_queue.tsv"),
          path("${prefix}.sequence_disposition.tsv"),
          path("${prefix}.sequence_map.tsv"),
          path("${prefix}.harmonizer_manifest.json"),
          path("${prefix}.database_candidates.fasta"),
          path("${prefix}.database_candidates.tsv"),
          path("${prefix}.all_candidates.fasta"),
          path("${prefix}.review_candidates.fasta"),
          path("${prefix}.provisional.normalized.fasta"),
          path("${prefix}.provisional.original_ids.fasta"),
          path("${prefix}.provisional_metadata.tsv"),
          emit: results
    tuple val(prefix),
          path("${prefix}.*_audit.tsv.gz"),
          optional: true,
          emit: audits

    script:
    def evidenceList = evidence_files instanceof List ? evidence_files : [evidence_files]
    def evidenceArguments = evidence_file_count > 0
        ? evidenceList.collect { evidence -> "'${evidence}'" }.join(' ')
        : ''
    def groupList = vcontact3_groups instanceof List ? vcontact3_groups : [vcontact3_groups]
    def groupArguments = group_file_count > 0
        ? groupList.collect { groups -> "'${groups}'" }.join(' ')
        : ''
    """
    # viHARMONY decision policy v0.6: refined regions remain subject to
    # viral-origin and mobile-element conflict adjudication.
    python3 "${projectDir}/bin/run_viharmony.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --normalized-fasta "${normalized_fasta}" \
        --header-map "${header_map}" \
        --discovery-gate "${discovery_gate}" \
        --refined-fasta "${refined_fasta}" \
        --region-map "${region_map}" \
        --ictv-msl "${ictv_msl}" \
        --evidence ${evidenceArguments} \
        --vcontact3-groups ${groupArguments} \
        --vcontact3-min-taxonomy-length "${vcontact3_min_taxonomy_length}" \
        --audit-mode "${audit_mode}" \
        --output-prefix "${prefix}"
    """
}
