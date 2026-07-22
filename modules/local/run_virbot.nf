process RUN_VIRBOT {

    tag "${prefix}"

    conda "${projectDir}/envs/virbot.yml"

    cpus params.threads
    memory params.virbot_memory
    time params.virbot_time

    publishDir { "${params.outdir}/${prefix}_results/virbot" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)
    tuple path(virbot_installation),
          path(virbot_database),
          path(virbot_installation_metadata),
          val(sensitive_mode),
          val(taxa_mode)

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.virbot_scores.csv"),
          path("${prefix}.virbot_virus_sequences.fasta"),
          path("${prefix}.virbot.log"),
          path("${prefix}.virbot_run_metadata.tsv"),
          emit: results

    script:
    """
    set -euo pipefail

    if [[ "${type}" != 'rna' ]]; then
        echo "ERROR: VirBot received non-RNA sample '${prefix}' with type '${type}'." >&2
        exit 1
    fi

    if [[ ! -s "${virbot_database}/VirBot.hmm" ]]; then
        echo "ERROR: VirBot's staged HMM database is missing or empty." >&2
        exit 1
    fi

    RAW_DIRECTORY="${prefix}.virbot_raw"
    RAW_SCORE_FILE="\$RAW_DIRECTORY/pos_contig_score.csv"
    RAW_FASTA_FILE="\$RAW_DIRECTORY/output.vb.fasta"
    FINAL_SCORE_FILE="${prefix}.virbot_scores.csv"
    FINAL_FASTA_FILE="${prefix}.virbot_virus_sequences.fasta"
    LOG_FILE="${prefix}.virbot.log"
    METADATA_FILE="${prefix}.virbot_run_metadata.tsv"

    SENSITIVE_FLAG=''
    if [[ "${sensitive_mode}" == 'true' ]]; then
        SENSITIVE_FLAG='--sen'
    fi

    if ! python "${virbot_installation}/virbot/VirBot.py" \
        --input "${normalized_fasta}" \
        --output "\$RAW_DIRECTORY" \
        --taxa "${taxa_mode}" \
        --threads "${task.cpus}" \
        \$SENSITIVE_FLAG \
        > "\$LOG_FILE" 2>&1; then
        echo "ERROR: VirBot prediction failed for sample '${prefix}'." >&2
        echo "Last 50 lines of \$LOG_FILE:" >&2
        tail -n 50 "\$LOG_FILE" >&2 || true
        exit 1
    fi

    if [[ ! -s "\$RAW_SCORE_FILE" ]]; then
        echo "ERROR: VirBot completed without pos_contig_score.csv." >&2
        exit 1
    fi
    if [[ ! -f "\$RAW_FASTA_FILE" ]]; then
        echo "ERROR: VirBot completed without output.vb.fasta." >&2
        exit 1
    fi

    EXPECTED_HEADER='Contig_acc,RNA-viral_gene_content,Encoded_proteins_num,Likely_taxa'
    ACTUAL_HEADER=\$(head -n 1 "\$RAW_SCORE_FILE" | tr -d '\r')
    if [[ "\$ACTUAL_HEADER" != "\$EXPECTED_HEADER" ]]; then
        echo "ERROR: VirBot produced an unexpected score-table header:" >&2
        echo "       \$ACTUAL_HEADER" >&2
        exit 1
    fi

    cp "\$RAW_SCORE_FILE" "\$FINAL_SCORE_FILE"
    cp "\$RAW_FASTA_FILE" "\$FINAL_FASTA_FILE"

    INPUT_COUNT=\$(awk 'NR > 1 && NF { count++ } END { print count + 0 }' "${header_map}")
    POSITIVE_COUNT=\$(awk 'NR > 1 && NF { count++ } END { print count + 0 }' "\$FINAL_SCORE_FILE")
    FASTA_POSITIVE_COUNT=\$(awk '/^>/ { count++ } END { print count + 0 }' "\$FINAL_FASTA_FILE")
    if [[ "\$POSITIVE_COUNT" -ne "\$FASTA_POSITIVE_COUNT" ]]; then
        echo "ERROR: VirBot score/FASTA call counts disagree for sample '${prefix}':" >&2
        echo "       score rows:     \$POSITIVE_COUNT" >&2
        echo "       FASTA records:  \$FASTA_POSITIVE_COUNT" >&2
        exit 1
    fi
    INSTALLED_REVISION=\$(
        awk -F '\t' 'NR == 2 { print \$8 }' "${virbot_installation_metadata}"
    )
    VIRBOT_VERSION=\$(
        awk -F '\t' 'NR == 2 { print \$9 }' "${virbot_installation_metadata}"
    )
    [[ -n "\$INSTALLED_REVISION" ]] || INSTALLED_REVISION='unknown'
    [[ -n "\$VIRBOT_VERSION" ]] || VIRBOT_VERSION='unknown'

    RUN_STATUS='completed'
    if [[ "\$POSITIVE_COUNT" -eq 0 ]]; then
        RUN_STATUS='completed_no_virus_calls'
    fi

    printf 'sample_id\tinput_type\tinput_sequence_count\tpositive_sequence_count\tsensitive_mode\ttaxa_mode\tvirbot_version\tvirbot_revision\trun_status\tscore_file\tvirus_fasta\n' \
        > "\$METADATA_FILE"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$INPUT_COUNT" "\$POSITIVE_COUNT" \
        "${sensitive_mode}" "${taxa_mode}" "\$VIRBOT_VERSION" \
        "\$INSTALLED_REVISION" "\$RUN_STATUS" "\$FINAL_SCORE_FILE" \
        "\$FINAL_FASTA_FILE" \
        >> "\$METADATA_FILE"

    echo "VirBot sample=${prefix} type=${type} positives=\$POSITIVE_COUNT scores=\$FINAL_SCORE_FILE"
    """
}
