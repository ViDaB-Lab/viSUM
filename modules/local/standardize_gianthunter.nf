process STANDARDIZE_GIANTHUNTER {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/gianthunter" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(prediction_table),
          path(gene_annotations),
          path(run_metadata),
          path(header_map)
    path ictv_csv

    output:
    tuple val(prefix),
          val('gianthunter'),
          path("${prefix}.gianthunter_evidence.tsv"),
          emit: evidence

    script:
    """
    python3 "${projectDir}/bin/standardize_gianthunter.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --header-map "${header_map}" \
        --prediction-table "${prediction_table}" \
        --gene-annotations "${gene_annotations}" \
        --run-metadata "${run_metadata}" \
        --ictv-csv "${ictv_csv}" \
        --output "${prefix}.gianthunter_evidence.tsv"
    """
}
