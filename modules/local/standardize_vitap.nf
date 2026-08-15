process STANDARDIZE_VITAP {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'
    cpus 1
    memory params.vitap_standardizer_memory
    time params.vitap_standardizer_time

    publishDir { "${params.outdir}/${prefix}_results/vitap" },
        mode: 'copy',
        pattern: '*.vitap_*.*'

    input:
    tuple val(prefix),
          val(type),
          path(provirus_region_map),
          path(best_lineages),
          path(all_lineages),
          path(uniref_fallback),
          path(vitap_log),
          path(run_metadata)
    tuple path(vitap_database), path(vitap_database_metadata)

    output:
    tuple val(prefix),
          val('vitap'),
          path("${prefix}.vitap_evidence.tsv"),
          emit: evidence
    tuple val(prefix),
          path("${prefix}.vitap_standardization_audit.tsv"),
          emit: audit

    script:
    """
    set -euo pipefail

    VMR_CSV=\$(find -L "${vitap_database}" -maxdepth 1 -type f \
        -name '*.csv' -size +0c | sort | head -n 1)
    if [[ -z "\$VMR_CSV" || ! -s "\$VMR_CSV" ]]; then
        echo "ERROR: No VMR taxonomy CSV was found in ${vitap_database}." >&2
        exit 1
    fi

    python3 "${projectDir}/bin/standardize_vitap.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --region-map "${provirus_region_map}" \
        --best-lineages "${best_lineages}" \
        --all-lineages "${all_lineages}" \
        --uniref-fallback "${uniref_fallback}" \
        --run-metadata "${run_metadata}" \
        --vmr "\$VMR_CSV" \
        --output-evidence "${prefix}.vitap_evidence.tsv" \
        --output-audit "${prefix}.vitap_standardization_audit.tsv"
    """
}
