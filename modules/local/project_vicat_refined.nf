process PROJECT_VICAT_REFINED {

    tag "${prefix}"
    conda "${projectDir}/envs/vicat.yml"
    cpus 1
    memory params.vicat_standardizer_memory
    time params.vicat_standardizer_time

    publishDir { "${params.outdir}/${prefix}_results/vicat" },
        mode: 'copy',
        pattern: '*.vicat_refined_*'

    input:
    tuple val(prefix), val(type), path(region_map), path(vicat_loci)

    output:
    tuple val(prefix), val('vicat'), path("${prefix}.vicat_refined_evidence.tsv"), emit: evidence
    tuple val(prefix), path("${prefix}.vicat_refined_projection.tsv"), emit: projection

    script:
    """
    python "${projectDir}/bin/project_vicat_refined.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --loci "${vicat_loci}" \
        --region-map "${region_map}" \
        --contig-taxonomy-support "${params.vicat_contig_taxonomy_support}" \
        --cluster-min-viral-loci "${params.vicat_cluster_min_viral_loci}" \
        --cluster-max-neutral-gap "${params.vicat_cluster_max_neutral_gap}" \
        --output-evidence "${prefix}.vicat_refined_evidence.tsv" \
        --output-projection "${prefix}.vicat_refined_projection.tsv"
    """
}
