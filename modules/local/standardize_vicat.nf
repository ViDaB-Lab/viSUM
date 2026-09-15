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
    path standardizer_script

    output:
    tuple val(prefix), val('vicat'), path("${prefix}.vicat_evidence.tsv"), emit: evidence
    tuple val(prefix), path("${prefix}.vicat_orf_evidence.tsv"), emit: loci
    tuple val(prefix), path("${prefix}.vicat_cluster_evidence.tsv"), emit: clusters
    tuple val(prefix), val('vicat_context'), path("${prefix}.vicat_context.tsv"), emit: context
    tuple val(prefix), val('vicat'), path("${prefix}.vicat_provirus_evidence.tsv"), emit: provirus
    tuple val(prefix), path("${prefix}.vicat_reference_audit.tsv"), emit: audit

    script:
    def configuredMinViralLoci = params.vicat_cluster_min_viral_loci.toString().trim()
    def effectiveMinViralLoci = configuredMinViralLoci.equalsIgnoreCase('auto') ?
        (type == 'rna' ? 1 : 2) : configuredMinViralLoci.toInteger()
    """
    # viCAT evidence policy v3: type-aware clean-locus threshold, guarded DNA
    # single-locus rescue, and a two-locus minimum in cellular context.
    set -euo pipefail

    python "${standardizer_script}" \
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
        --cluster-min-viral-loci "${effectiveMinViralLoci}" \
        --cluster-max-neutral-gap "${params.vicat_cluster_max_neutral_gap}" \
        --dna-single-locus-rescue "${params.vicat_dna_single_locus_rescue}" \
        --output-loci "${prefix}.vicat_orf_evidence.tsv" \
        --output-clusters "${prefix}.vicat_cluster_evidence.tsv" \
        --output-context "${prefix}.vicat_context.tsv" \
        --output-provirus-evidence "${prefix}.vicat_provirus_evidence.tsv" \
        --output-audit "${prefix}.vicat_reference_audit.tsv" \
        --output-evidence "${prefix}.vicat_evidence.tsv"
    """
}
