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
    tuple path(reference_subset), path(reference_subset_metadata)

    output:
    tuple val(prefix), val(type), path(orf_map), path(header_map),
          path("${prefix}.vicat_prepared_hits.parquet"),
          path("${prefix}.vicat_hit_preparation.tsv"),
          emit: results

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/prepare_vicat_hits.py" \
        --diamond "${diamond}" \
        --reference-subset "${reference_subset}" \
        --reference-subset-metadata "${reference_subset_metadata}" \
        --threads "${task.cpus}" \
        --memory-limit "${task.memory}" \
        --temp-directory "\$PWD/duckdb_tmp" \
        --output "${prefix}.vicat_prepared_hits.parquet" \
        --output-metadata "${prefix}.vicat_hit_preparation.tsv"
    """
}
