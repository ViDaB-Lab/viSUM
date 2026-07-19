process NORMALIZE_FASTA {

    tag "${prefix}"

    conda 'conda-forge::python=3.11'

    publishDir { "${params.outdir}/${prefix}/prep" }, mode: 'copy'

    input:
    tuple val(prefix), val(type), path(fasta)

    output:
    tuple val(prefix),
          val(type),
          path("${prefix}.normalized.fasta"),
          path("${prefix}.header_map.tsv"),
          emit: normalized_records

    script:
    """
    python3 "${projectDir}/bin/normalize_fasta.py" \
        --input "${fasta}" \
        --sample-id "${prefix}" \
        --output-fasta "${prefix}.normalized.fasta" \
        --output-map "${prefix}.header_map.tsv"
    """
}
