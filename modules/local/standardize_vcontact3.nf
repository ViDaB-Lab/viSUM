process STANDARDIZE_VCONTACT3 {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'
    cpus 1
    memory params.vcontact3_standardizer_memory
    time params.vcontact3_standardizer_time

    publishDir { "${params.outdir}/${prefix}_results/vcontact3" },
        mode: 'copy',
        pattern: '*.vcontact3_*.tsv'

    input:
    tuple val(prefix),
          val(type),
          path(provirus_region_map),
          path(assignments),
          path(performance_metrics),
          path(vcontact3_logs),
          path(run_metadata)

    output:
    tuple val(prefix),
          val('vcontact3'),
          path("${prefix}.vcontact3_evidence.tsv"),
          emit: evidence
    tuple val(prefix),
          path("${prefix}.vcontact3_group_membership.tsv"),
          emit: groups

    script:
    def assignmentFiles = assignments instanceof List ? assignments : [assignments]
    def assignmentArgs = assignmentFiles
        .collect { assignment -> "--assignments \"${assignment}\"" }
        .join(' ')
    """
    set -euo pipefail

    python3 "${projectDir}/bin/standardize_vcontact3.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --region-map "${provirus_region_map}" \
        ${assignmentArgs} \
        --run-metadata "${run_metadata}" \
        --output-evidence "${prefix}.vcontact3_evidence.tsv" \
        --output-groups "${prefix}.vcontact3_group_membership.tsv"
    """
}
