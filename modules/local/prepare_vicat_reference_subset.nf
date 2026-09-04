process PREPARE_VICAT_REFERENCE_SUBSET {

    tag "run-wide reference subset"
    conda "${projectDir}/envs/vicat.yml"
    cpus { Math.min(params.vicat_standardizer_cpus as int, params.max_cpus as int) }
    memory params.vicat_standardizer_memory
    time params.vicat_standardizer_time

    publishDir { "${params.outdir}/database_setup" },
        mode: 'copy',
        pattern: 'vicat_reference_subset_metadata.tsv'

    input:
    path diamonds
    tuple path(vicat_database), path(vicat_database_metadata)

    output:
    tuple path('vicat_reference_subset.parquet'),
          path('vicat_reference_subset_metadata.tsv'),
          emit: subset

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

    python "${projectDir}/bin/prepare_vicat_reference_subset.py" \
        --diamond ${diamonds} \
        --taxonomy-lookup "\$LOOKUP" \
        --threads "${task.cpus}" \
        --memory-limit "${task.memory}" \
        --temp-directory "\$PWD/duckdb_tmp" \
        --output vicat_reference_subset.parquet \
        --output-metadata vicat_reference_subset_metadata.tsv
    """
}
