process PREPARE_VICAT_HITS {

    tag "${prefix}"
    conda "${projectDir}/envs/vicat.yml"
    cpus { Math.min(params.vicat_standardizer_cpus as int, params.max_cpus as int) }
    memory params.vicat_standardizer_memory
    time params.vicat_standardizer_time

    publishDir { "${params.outdir}/${prefix}_results/vicat" },
        mode: 'copy',
        pattern: '*.vicat_hit_preparation.tsv'

    input:
    tuple val(prefix), val(type), path(orf_map), path(header_map), path(diamond)
    tuple path(vicat_database), path(vicat_database_metadata)

    output:
    tuple val(prefix), val(type), path(orf_map), path(header_map),
          path("${prefix}.vicat_prepared_hits.parquet"),
          path("${prefix}.vicat_hit_preparation.tsv"),
          emit: results

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

    MANIFEST_ARGS=()
    if [[ -s "${vicat_database}/vicat_competitive_reference_manifest.parquet" ]]; then
        MANIFEST_ARGS=(--reference-manifest "${vicat_database}/vicat_competitive_reference_manifest.parquet")
    fi

    python "${projectDir}/bin/prepare_vicat_hits.py" \
        --diamond "${diamond}" \
        --taxonomy-lookup "\$LOOKUP" \
        "\${MANIFEST_ARGS[@]}" \
        --threads "${task.cpus}" \
        --memory-limit "${task.memory}" \
        --temp-directory "\$PWD/duckdb_tmp" \
        --output "${prefix}.vicat_prepared_hits.parquet" \
        --output-metadata "${prefix}.vicat_hit_preparation.tsv"
    """
}
