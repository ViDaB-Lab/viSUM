process RUN_CHECKV {

    tag "${prefix}"

    conda "${projectDir}/envs/checkv.yml"

    cpus { Math.min(params.checkv_cpus as int, params.max_cpus as int) }
    memory params.checkv_memory
    time params.checkv_time

    publishDir { "${params.outdir}/${prefix}_results/checkv" },
        mode: 'copy'

    input:
    tuple val(prefix),
          val(type),
          path(candidate_fasta),
          path(discovery_audit),
          path(discovery_summary)
    tuple path(checkv_database), path(checkv_database_metadata)

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.checkv_quality_summary.tsv"),
          path("${prefix}.checkv_completeness.tsv"),
          path("${prefix}.checkv_contamination.tsv"),
          path("${prefix}.checkv_complete_genomes.tsv"),
          path("${prefix}.checkv_proviruses.fna"),
          path("${prefix}.checkv.log"),
          path("${prefix}.checkv_run_metadata.tsv"),
          emit: results
    path("${prefix}.checkv_raw"), emit: raw_output

    script:
    """
    set -euo pipefail

    RAW_DIRECTORY="${prefix}.checkv_raw"
    QUALITY_SUMMARY="${prefix}.checkv_quality_summary.tsv"
    COMPLETENESS="${prefix}.checkv_completeness.tsv"
    CONTAMINATION="${prefix}.checkv_contamination.tsv"
    COMPLETE_GENOMES="${prefix}.checkv_complete_genomes.tsv"
    PROVIRUSES="${prefix}.checkv_proviruses.fna"
    LOG_FILE="${prefix}.checkv.log"
    METADATA_FILE="${prefix}.checkv_run_metadata.tsv"

    mkdir -p "\$RAW_DIRECTORY"

    INPUT_COUNT=\$(awk '/^>/ { count++ } END { print count + 0 }' "${candidate_fasta}")

    if [[ "\$INPUT_COUNT" -eq 0 ]]; then
        printf '%s\n' \
            'contig_id	contig_length	provirus	proviral_length	gene_count	viral_genes	host_genes	checkv_quality	miuvig_quality	completeness	completeness_method	complete_genome_type	contamination	kmer_freq	warnings' \
            > "\$RAW_DIRECTORY/quality_summary.tsv"
        printf '%s\n' \
            'contig_id	contig_length	proviral_length	aai_expected_length	aai_completeness	aai_confidence	aai_error	aai_num_hits	aai_top_hit	aai_id	aai_af	hmm_completeness_lower	hmm_completeness_upper	hmm_hits' \
            > "\$RAW_DIRECTORY/completeness.tsv"
        printf '%s\n' \
            'contig_id	contig_length	total_genes	viral_genes	host_genes	provirus	proviral_length	host_length	region_types	region_lengths	region_coords_bp	region_coords_genes	region_viral_genes	region_host_genes' \
            > "\$RAW_DIRECTORY/contamination.tsv"
        printf '%s\n' \
            'contig_id	contig_length	prediction_type	confidence_level	confidence_reason	repeat_length	repeat_count' \
            > "\$RAW_DIRECTORY/complete_genomes.tsv"
        : > "\$RAW_DIRECTORY/proviruses.fna"
        printf 'CheckV skipped: the discovery gate produced no candidate sequences.\n' \
            > "\$LOG_FILE"
        RUN_STATUS='skipped_no_discovery_candidates'
    else
        if ! checkv end_to_end \
            "${candidate_fasta}" \
            "\$RAW_DIRECTORY" \
            -d "${checkv_database}" \
            -t "${task.cpus}" \
            > "\$LOG_FILE" 2>&1; then
            echo "ERROR: CheckV analysis failed for sample '${prefix}'." >&2
            echo "Last 50 lines of \$LOG_FILE:" >&2
            tail -n 50 "\$LOG_FILE" >&2 || true
            exit 1
        fi
        RUN_STATUS='completed'
    fi

    for expected_file in \
        quality_summary.tsv \
        completeness.tsv \
        contamination.tsv \
        complete_genomes.tsv; do
        if [[ ! -s "\$RAW_DIRECTORY/\$expected_file" ]]; then
            echo "ERROR: CheckV completed without required output: \$expected_file" >&2
            exit 1
        fi
    done

    # Some successful runs contain no predicted proviruses. Preserve a valid
    # empty file instead of treating that biological outcome as an error.
    if [[ ! -f "\$RAW_DIRECTORY/proviruses.fna" ]]; then
        : > "\$RAW_DIRECTORY/proviruses.fna"
    fi

    QUALITY_ROW_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' \
        "\$RAW_DIRECTORY/quality_summary.tsv")
    if [[ "\$QUALITY_ROW_COUNT" -ne "\$INPUT_COUNT" ]]; then
        echo "ERROR: CheckV quality-summary coverage mismatch for sample '${prefix}':" >&2
        echo "       candidate sequences: \$INPUT_COUNT" >&2
        echo "       quality rows:        \$QUALITY_ROW_COUNT" >&2
        exit 1
    fi

    cp "\$RAW_DIRECTORY/quality_summary.tsv" "\$QUALITY_SUMMARY"
    cp "\$RAW_DIRECTORY/completeness.tsv" "\$COMPLETENESS"
    cp "\$RAW_DIRECTORY/contamination.tsv" "\$CONTAMINATION"
    cp "\$RAW_DIRECTORY/complete_genomes.tsv" "\$COMPLETE_GENOMES"
    cp "\$RAW_DIRECTORY/proviruses.fna" "\$PROVIRUSES"

    # Keep the primary CheckV reports but remove regenerable intermediates.
    rm -rf "\$RAW_DIRECTORY/tmp"

    PROVIRUS_COUNT=\$(awk -F '\t' \
        'NR > 1 && \$3 == "Yes" { count++ } END { print count + 0 }' \
        "\$QUALITY_SUMMARY")
    DETERMINED_QUALITY_COUNT=\$(awk -F '\t' \
        'NR > 1 && \$8 != "Not-determined" && \$8 != "" { count++ } END { print count + 0 }' \
        "\$QUALITY_SUMMARY")
    CHECKV_VERSION=\$(checkv --version 2>&1 | head -n 1)
    DATABASE_RELEASE=\$(awk -F '\t' 'NR == 2 { print \$5 }' \
        "${checkv_database_metadata}")

    printf 'sample_id\tinput_type\tcandidate_sequence_count\tquality_summary_row_count\tdetermined_quality_count\tprovirus_count\tthreads\tcheckv_version\tdatabase_path\tdatabase_release\trun_status\tinput_fasta\n' \
        > "\$METADATA_FILE"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$INPUT_COUNT" "\$QUALITY_ROW_COUNT" \
        "\$DETERMINED_QUALITY_COUNT" "\$PROVIRUS_COUNT" "${task.cpus}" \
        "\$CHECKV_VERSION" "${checkv_database}" "\$DATABASE_RELEASE" \
        "\$RUN_STATUS" "${candidate_fasta.name}" \
        >> "\$METADATA_FILE"

    echo "CHECKV sample=${prefix} type=${type} candidates=\$INPUT_COUNT quality_rows=\$QUALITY_ROW_COUNT proviruses=\$PROVIRUS_COUNT"
    """
}
