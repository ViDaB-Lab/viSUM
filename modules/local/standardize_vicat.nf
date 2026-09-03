process STANDARDIZE_VICAT {

    tag "${prefix}"
    conda "${projectDir}/envs/vicat.yml"
    cpus { Math.min(params.vicat_standardizer_cpus as int, params.max_cpus as int) }
    memory params.vicat_standardizer_memory
    time params.vicat_standardizer_time

    publishDir { "${params.outdir}/${prefix}_results/vicat" },
        mode: 'copy',
        pattern: '*.vicat_*'

    input:
    tuple val(prefix), val(type), path(orf_map), path(header_map),
          path(prepared_hits), path(prepared_metadata)

    output:
    tuple val(prefix), val('vicat'), path("${prefix}.vicat_evidence.tsv"), emit: evidence
    tuple val(prefix), path("${prefix}.vicat_orf_evidence.tsv"), emit: loci
    tuple val(prefix), path("${prefix}.vicat_cluster_evidence.tsv"), emit: clusters
    tuple val(prefix), val('vicat_context'), path("${prefix}.vicat_context.tsv"), emit: context
    tuple val(prefix), path("${prefix}.vicat_reference_audit.tsv"), emit: audit

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/standardize_vicat.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --orf-map "${orf_map}" \
        --prepared-hits "${prepared_hits}" \
        --prepared-metadata "${prepared_metadata}" \
        --header-map "${header_map}" \
        --threads "${task.cpus}" \
        --orf-taxonomy-support "${params.vicat_orf_taxonomy_support}" \
        --contig-taxonomy-support "${params.vicat_contig_taxonomy_support}" \
        --locus-overlap "${params.vicat_locus_overlap}" \
        --competitive-min-margin "${params.vicat_competitive_min_margin}" \
        --cluster-min-viral-loci "${params.vicat_cluster_min_viral_loci}" \
        --cluster-max-neutral-gap "${params.vicat_cluster_max_neutral_gap}" \
        --output-loci "${prefix}.vicat_orf_evidence.tsv" \
        --output-clusters "${prefix}.vicat_cluster_evidence.tsv" \
        --output-context "${prefix}.vicat_context.tsv" \
        --output-audit "${prefix}.vicat_reference_audit.tsv" \
        --output-evidence "${prefix}.vicat_evidence.tsv"
    """
}
