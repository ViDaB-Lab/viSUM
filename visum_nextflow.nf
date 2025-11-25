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

process ORF_prediction {

    conda 'bioconda::pyrodigal-gv'

    publishDir "${params.output_dir}/ORF_prediction", mode: 'copy'

    input:
        tuple val(sampleID), file(fasta) 

    output:
        tuple val(sampleID), file("*.faa"), file("*.ffn"), file("*.gff")

    script:
    """
    pyrodigal-gv -i $fasta -a ${params.sampleID}_proteins.faa -d ${params.sampleID}_genes.ffn -f gff -o ${params.sampleID}_orfs.gff -j ${params.threads}
    """
}

process genomad {

    conda 'bioconda::genomad'

    publishDir "${params.output_dir}/genomad", mode: 'copy'

    input:
        tuple val(sampleID), file(fasta), file(faa), file(ffn), file(gff)

    output:
        tuple val(sampleID),

    script:
    """
    
    """
}

workflow {
    /*
     * Get list of files to process
     */
    file_ch = channel.fromPath("${params.input_dir}/*.fasta")
        .map { file -> tuple(file.sampleName, file) }
        .view()
        
    // RUN ORF prediction
    orfs_ch = ORF_prediction(file_ch)
    // JOIN file_ch with ORF prediction outputs
    pred_input = file_ch
        .join(orfs_ch)
        .map { sampleID, fasta, faa, ffn, gff ->
            tuple(sampleID, fasta, faa, ffn, gff)
        }
    // RUN genomad
    genomad(pred_input)
    // RUN genomad


}


