process RUN_VIRSORTER2 {

    tag "${prefix}"

    conda "bioconda::virsorter=${params.virsorter2_version}"

    cpus { Math.min(params.virsorter2_cpus as int, params.max_cpus as int) }
    memory params.virsorter2_memory
    time params.virsorter2_time

    publishDir { "${params.outdir}/${prefix}_results/virsorter2" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)
    path virsorter2_db

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.virsorter2"),
          path("${prefix}.virsorter2_final-viral-score.tsv"),
          path("${prefix}.virsorter2_final-viral-boundary.tsv"),
          path("${prefix}.virsorter2_final-viral-combined.fasta"),
          path("${prefix}.virsorter2_run_metadata.tsv"),
          emit: results

    script:
    def groups = type == 'rna' ? params.vs2_groups_rna : params.vs2_groups_dna

    """
    set -euo pipefail

    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export NUMEXPR_NUM_THREADS="${task.cpus}"

    ln -s "${normalized_fasta}" "${prefix}.fasta"

    virsorter run \
        -w "${prefix}.virsorter2" \
        -i "${prefix}.fasta" \
        -d "${virsorter2_db}" \
        --include-groups "${groups}" \
        --min-length ${params.vs2_min_length} \
        --min-score ${params.vs2_min_score} \
        --keep-original-seq \
        -j ${task.cpus} \
        all

    SCORE_FILE="${prefix}.virsorter2/final-viral-score.tsv"
    BOUNDARY_FILE="${prefix}.virsorter2/final-viral-boundary.tsv"
    FASTA_FILE="${prefix}.virsorter2/final-viral-combined.fa"

    for expected_file in "\$SCORE_FILE" "\$BOUNDARY_FILE" "\$FASTA_FILE"; do
        if [[ ! -f "\$expected_file" ]]; then
            echo "ERROR: VirSorter2 completed without expected output: \$expected_file" >&2
            exit 1
        fi
    done

    cp "\$SCORE_FILE" "${prefix}.virsorter2_final-viral-score.tsv"
    cp "\$BOUNDARY_FILE" "${prefix}.virsorter2_final-viral-boundary.tsv"
    cp "\$FASTA_FILE" "${prefix}.virsorter2_final-viral-combined.fasta"

    VIRUS_CALL_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' \
        "${prefix}.virsorter2_final-viral-score.tsv")
    if [[ "\$VIRUS_CALL_COUNT" -eq 0 ]]; then
        RUN_STATUS='completed_no_viruses_detected'
    else
        RUN_STATUS='completed_with_virus_calls'
    fi

    VIRSORTER2_VERSION=\$(
        python3 -c 'import virsorter; print(virsorter.__version__)' 2>/dev/null \
            || true
    )
    if [[ -z "\$VIRSORTER2_VERSION" ]]; then
        VIRSORTER2_VERSION='unknown'
    fi

    printf 'sample_id\tinput_type\tinput_fasta\tvirsorter2_version\tdatabase_path\tclassifier_groups\tmin_length\tmin_score\tkeep_original_sequence\tthreads\trun_status\tvirus_call_count\n' \
        > "${prefix}.virsorter2_run_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" \
        "${type}" \
        "${normalized_fasta.name}" \
        "\$VIRSORTER2_VERSION" \
        "${virsorter2_db}" \
        "${groups}" \
        "${params.vs2_min_length}" \
        "${params.vs2_min_score}" \
        'true' \
        "${task.cpus}" \
        "\$RUN_STATUS" \
        "\$VIRUS_CALL_COUNT" \
        >> "${prefix}.virsorter2_run_metadata.tsv"
    """
}
