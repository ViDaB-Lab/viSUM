process VALIDATE_PREPARED_ASSETS {
    tag 'read-only preflight'
    conda 'conda-forge::python=3.11'
    cache false
    cpus 1
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: 'asset_validation.json'

    input:
    tuple val(samples), val(encoded_spec)

    output:
    val samples, emit: samples
    path 'asset_validation.json', emit: report

    script:
    """
    python3 -B "${projectDir}/bin/validate_prepared_assets.py" \
        --spec-base64 "${encoded_spec}" --report asset_validation.json
    """
}
