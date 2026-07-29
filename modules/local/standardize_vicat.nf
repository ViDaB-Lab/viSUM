process STANDARDIZE_VICAT {

    tag "${prefix}"
    conda "${projectDir}/envs/vicat.yml"
    cpus 2
    memory params.vicat_standardizer_memory
    time params.vicat_standardizer_time

    publishDir { "${params.outdir}/${prefix}_results/vicat" },
        mode: 'copy',
        pattern: '*.vicat_*'

    input:
    tuple val(prefix), val(type), path(orf_map), path(header_map), path(diamond)
    tuple path(vicat_database), path(vicat_database_metadata)

    output:
    tuple val(prefix), val('vicat'), path("${prefix}.vicat_evidence.tsv"), emit: evidence
    tuple val(prefix), path("${prefix}.vicat_orf_evidence.tsv"), emit: loci

    script:
    """
    set -euo pipefail

    if [[ -s "${vicat_database}/IMGVR5_UViG.vicat_taxonomy_lookup.parquet" ]]; then
        LOOKUP="${vicat_database}/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
    elif [[ -s "${vicat_database}/taxonomy/IMGVR5_UViG.vicat_taxonomy_lookup.parquet" ]]; then
        LOOKUP="${vicat_database}/taxonomy/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
    else
        echo "ERROR: viCAT taxonomy lookup was not found in ${vicat_database}." >&2
        exit 1
    fi

    python "${projectDir}/bin/standardize_vicat.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --orf-map "${orf_map}" \
        --diamond "${diamond}" \
        --taxonomy-lookup "\$LOOKUP" \
        --header-map "${header_map}" \
        --taxonomy-support "${params.vicat_taxonomy_support}" \
        --locus-overlap "${params.vicat_locus_overlap}" \
        --output-loci "${prefix}.vicat_orf_evidence.tsv" \
        --output-evidence "${prefix}.vicat_evidence.tsv"
    """
}
