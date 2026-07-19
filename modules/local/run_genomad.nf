process RUN_GENOMAD {

    tag "${prefix}"

    conda "bioconda::genomad=${params.genomad_version}"

    cpus params.threads
    memory params.memory
    time params.time

    publishDir { "${params.outdir}/${prefix}_results/genomad" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)
    path genomad_db

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.genomad"),
          path("${prefix}_virus_summary.tsv"),
          path("${prefix}_virus.fna"),
          path("${prefix}_virus_genes.tsv"),
          path("${prefix}_virus_proteins.faa"),
          path("${prefix}_plasmid_summary.tsv"),
          path("${prefix}.genomad_run_metadata.tsv"),
          emit: results

    script:
    def cleanupFlag = params.genomad_cleanup ? '--cleanup' : ''
    """
    set -euo pipefail

    ln -s "${normalized_fasta}" "${prefix}.fasta"

    genomad end-to-end \
        ${cleanupFlag} \
        --threads ${task.cpus} \
        --splits ${params.genomad_splits} \
        "${prefix}.fasta" \
        "${prefix}.genomad" \
        "${genomad_db}"

    SUMMARY_DIR="${prefix}.genomad/${prefix}_summary"

    if [[ ! -d "\$SUMMARY_DIR" ]]; then
        echo "ERROR: geNomad completed without the expected summary directory:" >&2
        echo "       \$SUMMARY_DIR" >&2
        exit 1
    fi

    for expected_file in \
        "${prefix}_virus_summary.tsv" \
        "${prefix}_virus.fna" \
        "${prefix}_virus_genes.tsv" \
        "${prefix}_virus_proteins.faa" \
        "${prefix}_plasmid_summary.tsv"
    do
        if [[ ! -f "\$SUMMARY_DIR/\$expected_file" ]]; then
            echo "ERROR: Missing expected geNomad output: \$SUMMARY_DIR/\$expected_file" >&2
            exit 1
        fi
        cp "\$SUMMARY_DIR/\$expected_file" "\$expected_file"
    done

    printf 'sample_id\tinput_type\tinput_fasta\tgenomad_version\tthreads\tsplits\tcleanup\n' \
        > "${prefix}.genomad_run_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" \
        "${type}" \
        "${normalized_fasta.name}" \
        "\$(genomad --version 2>&1 | head -n 1)" \
        "${task.cpus}" \
        "${params.genomad_splits}" \
        "${params.genomad_cleanup}" \
        >> "${prefix}.genomad_run_metadata.tsv"
    """
}
