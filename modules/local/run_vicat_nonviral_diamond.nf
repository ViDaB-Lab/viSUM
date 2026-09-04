process RUN_VICAT_NONVIRAL_DIAMOND {

    tag "${prefix}"
    conda "${projectDir}/envs/vicat.yml"
    cpus { Math.min(params.vicat_cpus as int, params.max_cpus as int) }
    memory params.vicat_memory
    time params.vicat_time

    publishDir { "${params.outdir}/${prefix}_results/vicat" },
        mode: 'copy',
        pattern: '*.vicat_nonviral_*',
        saveAs: { filename -> filename.contains('.raw.') ? null : filename }

    input:
    tuple val(prefix), val(type), path(orfs_faa), path(orf_map), path(header_map)
    tuple path(vicat_nonviral_database), path(vicat_nonviral_database_metadata)

    output:
    tuple val(prefix), val(type), path(orf_map), path(header_map),
          path("${prefix}.vicat_nonviral_diamond.tsv"),
          emit: results
    path("${prefix}.vicat_nonviral_run_metadata.tsv"), emit: metadata
    path("${prefix}.vicat_nonviral_diamond.log"), emit: log

    script:
    """
    set -euo pipefail

    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export NUMEXPR_NUM_THREADS="${task.cpus}"

    DB="${vicat_nonviral_database}/vicat_nonviral.dmnd"
    [[ -s "\$DB" ]] || {
        echo "ERROR: viCAT nonviral DIAMOND database was not found: \$DB" >&2
        exit 1
    }

    RAW="${prefix}.vicat_nonviral_diamond.raw.tsv"
    LOG="${prefix}.vicat_nonviral_diamond.log"
    ORF_COUNT=\$(grep -c '^>' "${orfs_faa}" || true)

    if [[ "\$ORF_COUNT" -eq 0 ]]; then
        : > "\$RAW"
        printf 'viCAT nonviral DIAMOND skipped: no predicted ORFs.\n' > "\$LOG"
    else
        diamond blastp --query "${orfs_faa}" --db "\$DB" --out "\$RAW" \
            --outfmt 6 qseqid sseqid pident length qlen slen qstart qend \
            sstart send evalue bitscore qcovhsp scovhsp \
            --${params.vicat_diamond_sensitivity} \
            --min-score "${params.vicat_min_bitscore}" \
            --query-cover "${params.vicat_min_query_cover}" \
            --top "${params.vicat_top_percent}" \
            --block-size "${params.vicat_block_size}" \
            --index-chunks "${params.vicat_index_chunks}" \
            --threads "${task.cpus}" > "\$LOG" 2>&1
    fi

    printf 'qseqid\tsseqid\tpident\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovhsp\tscovhsp\n' \
        > "${prefix}.vicat_nonviral_diamond.tsv"
    cat "\$RAW" >> "${prefix}.vicat_nonviral_diamond.tsv"
    HIT_COUNT=\$(wc -l < "\$RAW")
    HIT_ORFS=\$(awk -F '\t' 'NF {seen[\$1]=1} END {print length(seen)+0}' "\$RAW")

    printf 'sample_id\tinput_type\tdatabase_role\torf_count\thit_orf_count\talignment_count\tdatabase_path\tsensitivity\tminimum_bitscore\tminimum_query_cover\ttop_percent\tblock_size\tindex_chunks\tthreads\n' \
        > "${prefix}.vicat_nonviral_run_metadata.tsv"
    printf '%s\t%s\tnonviral\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$ORF_COUNT" "\$HIT_ORFS" "\$HIT_COUNT" \
        "\$DB" "${params.vicat_diamond_sensitivity}" \
        "${params.vicat_min_bitscore}" "${params.vicat_min_query_cover}" \
        "${params.vicat_top_percent}" "${params.vicat_block_size}" \
        "${params.vicat_index_chunks}" "${task.cpus}" \
        >> "${prefix}.vicat_nonviral_run_metadata.tsv"

    echo "VICAT_NONVIRAL_DIAMOND sample=${prefix} orfs=\$ORF_COUNT hit_orfs=\$HIT_ORFS alignments=\$HIT_COUNT"
    """
}
