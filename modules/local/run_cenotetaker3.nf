process RUN_CENOTETAKER3 {

    tag "${prefix}"

    // Match the database module exactly so Nextflow reuses one CT3 environment.
    conda "bioconda::cenote-taker3=${params.ct3_version} conda-forge::wget"

    cpus params.threads
    memory params.memory
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

    """
    set -euo pipefail

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

    for expected_file in \
        "\$SUMMARY_FILE" \
        "\$VIRUS_FASTA" \
        "\$VIRUS_PROTEINS" \
        "\$GENE_ANNOTATIONS" \
        "\$RUN_ARGUMENTS" \
        "\$CT3_LOG"
    do
        if [[ ! -f "\$expected_file" ]]; then
            echo "ERROR: Cenote-Taker 3 completed without expected output: \$expected_file" >&2
            echo "Review \$CT3_LOG and the task work directory for the underlying CT3 message." >&2
            exit 1
        fi
    done

    cp "\$SUMMARY_FILE" "${prefix}.cenotetaker3_virus_summary.tsv"
    cp "\$VIRUS_FASTA" "${prefix}.cenotetaker3_virus_sequences.fna"
    cp "\$VIRUS_PROTEINS" "${prefix}.cenotetaker3_virus_AA.faa"
    cp "\$GENE_ANNOTATIONS" "${prefix}.cenotetaker3_gene_annotations.tsv"
    cp "\$RUN_ARGUMENTS" "${prefix}.cenotetaker3_run_arguments.txt"
    cp "\$CT3_LOG" "${prefix}.cenotetaker3.log"

    if [[ -f "\$PRUNE_SUMMARY" ]]; then
        cp "\$PRUNE_SUMMARY" "${prefix}.cenotetaker3_prune_summary.tsv"
    else
        printf 'contig\tcontig_length\tchunk_length\tchunk_name\tchunk_start\tchunk_stop\n' \
            > "${prefix}.cenotetaker3_prune_summary.tsv"
    fi

    CT3_VERSION=\$(
        python3 -c "from importlib.metadata import version; print(version('cenote-taker3'))" \
            2>/dev/null || true
    )
    if [[ -z "\$CT3_VERSION" ]]; then
        CT3_VERSION='unknown'
    fi

    printf 'sample_id\tinput_type\tinput_fasta\tcenotetaker3_version\tdatabase_path\thmm_database_version\tmolecule_type\tprune_prophage\thallmark_databases\tcaller\ttaxonomy_database\thhsuite_tool\tminimum_length_circular\tminimum_hallmarks_circular\tminimum_length_linear\tminimum_hallmarks_linear\twrap\tgenbank\tthreads\n' \
        > "${prefix}.cenotetaker3_run_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
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
        >> "${prefix}.cenotetaker3_run_metadata.tsv"
    """
}
