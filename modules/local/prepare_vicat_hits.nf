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
    tuple val(prefix), val(type), path(orf_map), path(header_map),
          path(viral_diamond), path(nonviral_diamond)
    tuple path(reference_subset), path(reference_subset_metadata)
    tuple path(vicat_nonviral_database), path(vicat_nonviral_database_metadata)

    output:
    tuple val(prefix), val(type), path(orf_map), path(header_map),
          path("${prefix}.vicat_prepared_hits.parquet"),
          path("${prefix}.vicat_hit_preparation.tsv"),
          emit: results

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/prepare_vicat_hits.py" \
        --viral-diamond "${viral_diamond}" \
        --viral-reference-subset "${reference_subset}" \
        --viral-reference-subset-metadata "${reference_subset_metadata}" \
        --nonviral-diamond "${nonviral_diamond}" \
        --nonviral-metadata "${vicat_nonviral_database}/vicat_nonviral_representative_metadata.parquet" \
        --threads "${task.cpus}" \
        --memory-limit "${task.memory}" \
        --temp-directory "\$PWD/duckdb_tmp" \
        --output "${prefix}.vicat_prepared_hits.parquet" \
        --output-metadata "${prefix}.vicat_hit_preparation.tsv"
    """
}
