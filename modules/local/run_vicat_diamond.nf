process RUN_VICAT_DIAMOND {

    tag "${prefix}"
    conda "${projectDir}/envs/vicat.yml"
    cpus params.vicat_threads
    memory params.vicat_memory
    time params.vicat_time

    publishDir { "${params.outdir}/${prefix}_results/vicat" },
        mode: 'copy',
        pattern: "${prefix}.vicat_{diamond.tsv,diamond.log,run_metadata.tsv}"

    input:
    tuple val(prefix), val(type), path(orfs_faa), path(orf_map), path(header_map), path(orf_metadata), path(orf_log)
    tuple path(vicat_database), path(vicat_database_metadata)

    output:
    tuple val(prefix), val(type), path(orf_map), path(header_map),
          path("${prefix}.vicat_diamond.tsv"),
          path("${prefix}.vicat_run_metadata.tsv"),
          path("${prefix}.vicat_diamond.log"),
          emit: results

    script:
    """
    set -euo pipefail

    if [[ -s "${vicat_database}/IMGVR5_UViG_representatives.dmnd" ]]; then
        DB="${vicat_database}/IMGVR5_UViG_representatives.dmnd"
    elif [[ -s "${vicat_database}/diamond/IMGVR5_UViG_representatives.dmnd" ]]; then
        DB="${vicat_database}/diamond/IMGVR5_UViG_representatives.dmnd"
    else
        echo "ERROR: viCAT DIAMOND database was not found in ${vicat_database}." >&2
        exit 1
    fi

    FINAL="${prefix}.vicat_diamond.tsv"
    RAW="${prefix}.vicat_diamond.raw.tsv"
    LOG="${prefix}.vicat_diamond.log"
    HEADER=$'qseqid\tsseqid\tpident\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovhsp\tscovhsp'
    ORF_COUNT=\$(grep -c '^>' "${orfs_faa}" || true)

    if [[ "\$ORF_COUNT" -eq 0 ]]; then
        : > "\$RAW"
        printf 'viCAT DIAMOND skipped: no predicted ORFs.\n' > "\$LOG"
    else
        diamond blastp \
            --query "${orfs_faa}" \
            --db "\$DB" \
            --out "\$RAW" \
            --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
            --${params.vicat_diamond_sensitivity} \
            --evalue "${params.vicat_evalue}" \
            --min-score "${params.vicat_min_bitscore}" \
            --query-cover "${params.vicat_min_query_cover}" \
            --top "${params.vicat_top_percent}" \
            --threads "${task.cpus}" \
            > "\$LOG" 2>&1
    fi

    printf '%s\n' "\$HEADER" > "\$FINAL"
    cat "\$RAW" >> "\$FINAL"
    HIT_COUNT=\$(wc -l < "\$RAW")
    HIT_ORFS=\$(awk -F '\t' 'NF {seen[\$1]=1} END {print length(seen)+0}' "\$RAW")

    printf 'sample_id\tinput_type\torf_count\thit_orf_count\talignment_count\tdatabase_path\tsensitivity\tevalue\tminimum_bitscore\tminimum_query_cover\ttop_percent\n' \
        > "${prefix}.vicat_run_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$ORF_COUNT" "\$HIT_ORFS" "\$HIT_COUNT" \
        "\$DB" "${params.vicat_diamond_sensitivity}" "${params.vicat_evalue}" \
        "${params.vicat_min_bitscore}" "${params.vicat_min_query_cover}" \
        "${params.vicat_top_percent}" >> "${prefix}.vicat_run_metadata.tsv"

    echo "VICAT_DIAMOND sample=${prefix} orfs=\$ORF_COUNT hit_orfs=\$HIT_ORFS alignments=\$HIT_COUNT"
    """
}
