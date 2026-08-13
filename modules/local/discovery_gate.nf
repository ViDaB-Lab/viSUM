process DISCOVERY_GATE {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/discovery_gate" }, mode: 'copy'

    input:
    tuple val(prefix), val(type), path(fasta), path(header_map), path(evidence_files)

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.discovery_candidates.fasta"),
          path("${prefix}.discovery_gate.tsv"),
          path("${prefix}.discovery_gate_summary.tsv"),
          emit: candidates
    tuple val(prefix),
          val(type),
          path("${prefix}.discovery_noncandidates.fasta"),
          emit: noncandidates

    script:
    def evidenceList = evidence_files instanceof List ? evidence_files : [evidence_files]
    def evidenceArguments = evidenceList
        .collect { evidence -> "'${evidence}'" }
        .join(' ')
    """
    python3 "${projectDir}/bin/discovery_gate.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --fasta "${fasta}" \
        --header-map "${header_map}" \
        --evidence ${evidenceArguments} \
        --output-candidates "${prefix}.discovery_candidates.fasta" \
        --output-noncandidates "${prefix}.discovery_noncandidates.fasta" \
        --output-audit "${prefix}.discovery_gate.tsv" \
        --output-summary "${prefix}.discovery_gate_summary.tsv"
    """
}
