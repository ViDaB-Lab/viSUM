process STANDARDIZE_VIRBOT {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}_results/virbot" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(score_table),
          path(virus_fasta),
          path(run_metadata),
          path(header_map)
    path ictv_csv

    output:
    tuple val(prefix),
          val('virbot'),
          path("${prefix}.virbot_evidence.tsv"),
          emit: evidence

    script:
    """
    python3 "${projectDir}/bin/standardize_virbot.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --header-map "${header_map}" \
        --score-table "${score_table}" \
        --virus-fasta "${virus_fasta}" \
        --run-metadata "${run_metadata}" \
        --ictv-csv "${ictv_csv}" \
        --output "${prefix}.virbot_evidence.tsv"
    """
}
