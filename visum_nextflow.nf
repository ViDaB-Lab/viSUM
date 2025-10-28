#!/usr/bin/env nextflow
/*
 * ==========================================================
 *  viSUM (DSL2)
 *  Author: <Your Name>
 *  Project: <Project Name>
 *  Date: <Date>
 * ==========================================================
 */

nextflow.enable.dsl=2

/*
 * ----------------------------------------------------------
 *                  INPUT VALIDATION
 * ----------------------------------------------------------
 * Verify input paths and create output directory if needed.
 */
workflow.onComplete {
    log.info "Pipeline complete! Results saved to: ${params.output_dir}"
}

if (!file(params.input_dir).exists()) {
    error "❌ Input directory does not exist: ${params.input_dir}"
}

/*
 * ----------------------------------------------------------
 *                  PROCESS DEFINITIONS
 * ----------------------------------------------------------
 */

process kraken_cell_filtering {
    tag "$sample_id"

    publishDir "${params.output_dir}/kraken", mode: 'copy'

    input:
    path sample_file
    val sample_id

    output:
    path "${sample_id}.out"

    script:
    """
    kraken2 --db /work/databases/kraken2/k2pluspftpall/ --use-names --threads 10 --report ${base}_kraken_report.txt --output ${base}_kraken_classifications.txt  $x
    """
}

process kraken_processing {
    tag "$sample_id"

    publishDir "${params.output_dir}/kraken", mode: 'copy'

    input:
    path sample_file
    val sample_id

    output:
    path "${sample_id}.out"

    script:
    """
    for x in *_kraken_classifications.txt; do
        base=${x%_kraken_classifications.txt}
        fasta="${base}_all_seqs_10_25_nr_.fasta"   # adjust to match your actual FASTA naming
        awk '$1=="U"{print $2}' $x > ${base}_unclassified_ids.txt
        grep "Virus" $x > ${base}_viral_ids.txt
        grep "virus" $x >> ${base}_viral_ids.txt
        cat ${base}_unclassified_ids.txt ${base}_viral_ids.txt | sort -u > ${base}_keep_ids.txt
        seqtk subseq ${fasta} ${base}_keep_ids.txt > ${base}_viral_unclassified.fasta
    done
    """
}

/*
 * ----------------------------------------------------------
 *                  WORKFLOW DEFINITION
 * ----------------------------------------------------------
 */

workflow MAIN {
    /*
     * Get list of files to process
     */
    Channel
        .fromPath("${params.input_dir}/*.fasta")
        .map { file -> tuple(file.simpleName, file) }
        |>
        EXAMPLE_PROCESS
}

/*
 * ----------------------------------------------------------
 *                  END OF FILE
 * ----------------------------------------------------------
 */


