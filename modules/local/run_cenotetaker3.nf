process RUN_CENOTETAKER3 {

    tag "${prefix}"

    // Match the database module exactly so Nextflow reuses one CT3 environment.
    conda "bioconda::cenote-taker3=${params.ct3_version} conda-forge::wget"

    cpus { Math.min(params.ct3_cpus as int, params.max_cpus as int) }
    memory params.ct3_memory
    time params.ct3_time

    publishDir { "${params.outdir}/${prefix}_results/cenotetaker3" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)
    path cenotetaker3_db

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.cenotetaker3_virus_summary.tsv"),
          path("${prefix}.cenotetaker3_virus_sequences.fna"),
          path("${prefix}.cenotetaker3_virus_AA.faa"),
          path("${prefix}.cenotetaker3_prune_summary.tsv"),
          path("${prefix}.cenotetaker3_gene_annotations.tsv"),
          path("${prefix}.cenotetaker3_run_arguments.txt"),
          path("${prefix}.cenotetaker3.log"),
          path("${prefix}.cenotetaker3_run_metadata.tsv"),
          emit: results

    script:
    def runTitle = 'ct3_run'
    def moleculeType = type == 'rna' ? 'RNA' : 'DNA'
    def pruneProphage = type == 'rna' ? params.ct3_prune_prophage_rna : params.ct3_prune_prophage_dna
    def minimumInputLength = Math.min(
        params.ct3_minlen_circ as Integer,
        params.ct3_minlen_linear as Integer
    )

    """
    set -euo pipefail

    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export NUMEXPR_NUM_THREADS="${task.cpus}"

    cenotetaker3 \
        -c "${normalized_fasta}" \
        -r "${runTitle}" \
        -p "${pruneProphage}" \
        -t ${task.cpus} \
        -wd . \
        --cenote-dbs "${cenotetaker3_db}" \
        --hmmscan_dbs "${params.ct3_hmm_db_version}" \
        --molecule_type "${moleculeType}" \
        --caller "${params.ct3_caller}" \
        -hh "${params.ct3_hh}" \
        --taxdb "${params.ct3_taxdb}" \
        -db ${params.ct3_domaindb} \
        --minimum_length_circular ${params.ct3_minlen_circ} \
        --circ_minimum_hallmark_genes ${params.ct3_circ_minhall} \
        --minimum_length_linear ${params.ct3_minlen_linear} \
        --lin_minimum_hallmark_genes ${params.ct3_linear_minhall} \
        --wrap "${params.ct3_wrap}" \
        --genbank "${params.ct3_genbank}"

    CT3_DIR="${runTitle}"
    SUMMARY_FILE="\$CT3_DIR/${runTitle}_virus_summary.tsv"
    VIRUS_FASTA="\$CT3_DIR/${runTitle}_virus_sequences.fna"
    VIRUS_PROTEINS="\$CT3_DIR/${runTitle}_virus_AA.faa"
    PRUNE_SUMMARY="\$CT3_DIR/${runTitle}_prune_summary.tsv"
    GENE_ANNOTATIONS="\$CT3_DIR/final_genes_to_contigs_annotation_summary.tsv"
    RUN_ARGUMENTS="\$CT3_DIR/run_arguments.txt"
    CT3_LOG="\$CT3_DIR/${runTitle}_cenotetaker.log"
    FILTERED_FASTA="\$CT3_DIR/${runTitle}.contigs_over_${minimumInputLength}nt.fasta"
    CONTIGS_TO_KEEP="\$CT3_DIR/ct_processing/contigs_to_keep.txt"
    HALLMARK_COUNTS="\$CT3_DIR/ct_processing/hallmarks_per_orig_contigs.tsv"
    CONTIGS_OVER_THRESHOLD="\$CT3_DIR/ct_processing/contigs_over_threshold.txt"
    TERMINAL_REPEAT_SUMMARY="\$CT3_DIR/ct_processing/threshold_contigs_terminal_repeat_summary.tsv"

    for audit_file in "\$RUN_ARGUMENTS" "\$CT3_LOG"; do
        if [[ ! -f "\$audit_file" ]]; then
            echo "ERROR: Cenote-Taker 3 did not create required audit file: \$audit_file" >&2
            exit 1
        fi
    done

    if [[ -f "\$SUMMARY_FILE" ]]; then
        for call_file in "\$VIRUS_FASTA" "\$VIRUS_PROTEINS" "\$GENE_ANNOTATIONS"; do
            if [[ ! -f "\$call_file" ]]; then
                echo "ERROR: Cenote-Taker 3 created a summary but omitted: \$call_file" >&2
                exit 1
            fi
        done

        cp "\$SUMMARY_FILE" "${prefix}.cenotetaker3_virus_summary.tsv"
        cp "\$VIRUS_FASTA" "${prefix}.cenotetaker3_virus_sequences.fna"
        cp "\$VIRUS_PROTEINS" "${prefix}.cenotetaker3_virus_AA.faa"
        cp "\$GENE_ANNOTATIONS" "${prefix}.cenotetaker3_gene_annotations.tsv"

        if [[ -f "\$PRUNE_SUMMARY" ]]; then
            cp "\$PRUNE_SUMMARY" "${prefix}.cenotetaker3_prune_summary.tsv"
        else
            printf 'contig\tcontig_length\tchunk_length\tchunk_name\tchunk_start\tchunk_stop\n' \
                > "${prefix}.cenotetaker3_prune_summary.tsv"
        fi
    else
        NO_VIRUS_REASON=''

        if [[ -f "\$FILTERED_FASTA" && ! -s "\$FILTERED_FASTA" ]]; then
            NO_VIRUS_REASON='no_sequences_met_the_minimum_length'
        elif [[ -f "\$CONTIGS_TO_KEEP" && ! -s "\$CONTIGS_TO_KEEP" && -s "\$HALLMARK_COUNTS" ]]; then
            NO_VIRUS_REASON='no_sequences_met_the_hallmark_requirement'
        elif [[ -f "\$CONTIGS_OVER_THRESHOLD" && ! -s "\$CONTIGS_OVER_THRESHOLD" && -s "\$TERMINAL_REPEAT_SUMMARY" ]]; then
            NO_VIRUS_REASON='no_sequences_met_the_final_thresholds'
        fi

        if [[ -z "\$NO_VIRUS_REASON" ]]; then
            echo "ERROR: Cenote-Taker 3 did not create its final summary and no valid zero-virus checkpoint was found." >&2
            echo "Review \$CT3_LOG and the task work directory; this is treated as a tool failure, not a biological zero." >&2
            exit 1
        fi

        printf 'contig\tinput_name\torganism\tvirus_seq_length\tend_feature\tgene_count\tvirion_hallmark_count\trep_hallmark_count\tRDRP_hallmark_count\tvirion_hallmark_genes\trep_hallmark_genes\tRDRP_hallmark_genes\ttaxonomy_hierarchy\tORF_caller\tgcode\tavg_read_depth\n' \
            > "${prefix}.cenotetaker3_virus_summary.tsv"
        : > "${prefix}.cenotetaker3_virus_sequences.fna"
        : > "${prefix}.cenotetaker3_virus_AA.faa"
        printf 'contig\tcontig_length\tchunk_length\tchunk_name\tchunk_start\tchunk_stop\n' \
            > "${prefix}.cenotetaker3_prune_summary.tsv"
        printf 'contig\tgene_start\tgene_stop\tgene_name\tgene_orient\tcontig_length\tdtr_seq\tevidence_acession\tevidence_description\tEvidence_source\tvscore_category\tchunk_name\tchunk_length\tchunk_start\tchunk_stop\n' \
            > "${prefix}.cenotetaker3_gene_annotations.tsv"

        echo "Cenote-Taker 3 completed with zero virus calls: \$NO_VIRUS_REASON"
    fi

    cp "\$RUN_ARGUMENTS" "${prefix}.cenotetaker3_run_arguments.txt"
    cp "\$CT3_LOG" "${prefix}.cenotetaker3.log"

    VIRUS_CALL_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' \
        "${prefix}.cenotetaker3_virus_summary.tsv")
    if [[ "\$VIRUS_CALL_COUNT" -eq 0 ]]; then
        RUN_STATUS='completed_no_viruses_detected'
    else
        RUN_STATUS='completed_with_virus_calls'
    fi

    CT3_VERSION=\$(
        python3 -c "from importlib.metadata import version; print(version('cenotetaker3'))" \
            2>/dev/null || true
    )
    if [[ -z "\$CT3_VERSION" ]]; then
        CT3_VERSION='unknown'
    fi

    printf 'sample_id\tinput_type\tinput_fasta\tcenotetaker3_version\tdatabase_path\thmm_database_version\tmolecule_type\tprune_prophage\thallmark_databases\tcaller\ttaxonomy_database\thhsuite_tool\tminimum_length_circular\tminimum_hallmarks_circular\tminimum_length_linear\tminimum_hallmarks_linear\twrap\tgenbank\tthreads\trun_status\tvirus_call_count\n' \
        > "${prefix}.cenotetaker3_run_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" \
        "${type}" \
        "${normalized_fasta.name}" \
        "\$CT3_VERSION" \
        "${cenotetaker3_db}" \
        "${params.ct3_hmm_db_version}" \
        "${moleculeType}" \
        "${pruneProphage}" \
        "${params.ct3_domaindb}" \
        "${params.ct3_caller}" \
        "${params.ct3_taxdb}" \
        "${params.ct3_hh}" \
        "${params.ct3_minlen_circ}" \
        "${params.ct3_circ_minhall}" \
        "${params.ct3_minlen_linear}" \
        "${params.ct3_linear_minhall}" \
        "${params.ct3_wrap}" \
        "${params.ct3_genbank}" \
        "${task.cpus}" \
        "\$RUN_STATUS" \
        "\$VIRUS_CALL_COUNT" \
        >> "${prefix}.cenotetaker3_run_metadata.tsv"
    """
}
