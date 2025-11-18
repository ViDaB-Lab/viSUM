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

process kraken_cell_filtering {

    publishDir "${params.output_dir}/kraken", mode: 'copy'

    input:
        tuple val(sampleID), file(fasta) 

    output:
        path "${sampleID}.out"

    script:
    """
    echo "${fasta}   ${sampleID}" > "${sampleID}.out"
    """
}

workflow {
    /*
     * Get list of files to process
     */
    file_ch = channel.fromPath("${params.input_dir}/*.fasta")
        .map { file -> tuple(file.simpleName, file) }
        .view()
        

    kraken_cell_filtering(file_ch)
}


