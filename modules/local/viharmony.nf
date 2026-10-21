process VIHARMONY {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'
    cpus { Math.min(params.harmonizer_cpus as int, params.max_cpus as int) }
    memory params.harmonizer_memory
    time params.harmonizer_time

    // Referenced task.ext values are native Nextflow cache dependencies.
    // A script-local computed value alone is NOT a dependency fingerprint.
    ext harmonizer_source_digest: {
        def sources = ['run_viharmony.py', 'evidence_schema.py', 'ictv_reconciliation.py', 'te_routing.py', 'reporting_bundle.py']
        def sourceBytes = sources.collect { name -> file("${projectDir}/bin/${name}").text }.join('\n').getBytes('UTF-8')
        java.security.MessageDigest.getInstance('SHA-256').digest(sourceBytes).encodeHex().toString()
    }

    // A resumed reporting change must refresh the sample's published files,
    // including when restoring a previously cached rank cap.
    publishDir { "${params.outdir}/${prefix}_results/viharmony" }, mode: 'copy', overwrite: true, failOnError: true

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
          path(ictv_reconciliation),
          path(ictv_reconciliation_manifest),
          val(audit_mode),
          val(vcontact3_min_taxonomy_length)

    output:
    tuple val(prefix),
          path("audit/legacy/${prefix}.final.normalized.fasta"),
          path("audit/legacy/${prefix}.final.original_ids.fasta"),
          path("audit/${prefix}.final_metadata.tsv"),
          path("audit/${prefix}.review_queue.tsv"),
          path("audit/${prefix}.sequence_disposition.tsv"),
          path("audit/${prefix}.sequence_map.tsv"),
          path("audit/${prefix}.harmonizer_manifest.json"),
          path("audit/legacy/${prefix}.database_candidates.fasta"),
          path("audit/legacy/${prefix}.database_candidates.tsv"),
          path("audit/legacy/${prefix}.all_candidates.fasta"),
          path("audit/legacy/${prefix}.provisional.normalized.fasta"),
          path("audit/legacy/${prefix}.provisional.original_ids.fasta"),
          path("audit/legacy/${prefix}.provisional_metadata.tsv"),
          emit: results
    tuple val(prefix),
          path("audit/${prefix}.*_audit.tsv.gz"),
          optional: true,
          emit: audits

    tuple val(prefix),
          path("TE_candidates/${prefix}.TE_candidates.fasta"),
          path("TE_candidates/${prefix}.TE_candidates.summary.tsv"),
          path("${prefix}.viral_candidates.fasta"),
          path("${prefix}.viral_candidates.summary.tsv"),
          emit: candidate_reports

    tuple val(prefix),
          path("${prefix}.taxonomy.tsv"),
          path("quality/${prefix}.checkv.tsv"),
          path("groups/${prefix}.vcontact3.tsv"),
          path("review/${prefix}.review.tsv"),
          emit: reporting_tables
    tuple val(prefix),
          path("review/by_reason/${prefix}.review.*.tsv"),
          optional: true,
          emit: review_views

    script:
    def evidenceList = evidence_files instanceof List ? evidence_files : [evidence_files]
    def maxRank = params.harmonizer_max_rank.toString().trim().toLowerCase()
    def pairFlag = type == 'rna' && params.rna_pair_homology_floor.toString().equalsIgnoreCase('true')
        ? '--rna-pair-homology-floor' : ''
    def evidenceArguments = evidence_file_count > 0
        ? evidenceList.collect { evidence -> "'${evidence}'" }.join(' ')
        : ''
    def groupList = vcontact3_groups instanceof List ? vcontact3_groups : [vcontact3_groups]
    def groupArguments = group_file_count > 0
        ? groupList.collect { groups -> "'${groups}'" }.join(' ')
        : ''
    """
    # viHARMONY source fingerprint: ${task.ext.harmonizer_source_digest}
    # TE-routing-v1 with clarified evidence scope; saved-input parity verified separately.
    # viHARMONY decision policy v1.4: audited RNA CheckV-only conflict exception;
    # VS2 hallmark-only calls are review-only;
    # do not bypass the RNA pair floor or project unlocalized calls onto children.
    # Global MSL taxonomy resolution and
    # conservative strict species reporting; disputed CT3/viCAT boundaries require
    # regional support; parent-only cellular/plasmid evidence remains context,
    # and strongly supported virus/plasmid hybrids retain both interpretations.
    python3 "${projectDir}/bin/run_viharmony.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --normalized-fasta "${normalized_fasta}" \
        --header-map "${header_map}" \
        --discovery-gate "${discovery_gate}" \
        --refined-fasta "${refined_fasta}" \
        --region-map "${region_map}" \
        --ictv-msl "${ictv_msl}" \
        --ictv-reconciliation "${ictv_reconciliation}" \
        --evidence ${evidenceArguments} \
        --vcontact3-groups ${groupArguments} \
        --vcontact3-min-taxonomy-length "${vcontact3_min_taxonomy_length}" \
        --audit-mode "${audit_mode}" \
        --max-rank "${maxRank}" \
        --workers "${task.cpus}" \
        ${pairFlag} \
        --output-prefix "${prefix}"
    """
}
