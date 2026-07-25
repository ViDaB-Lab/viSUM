process PREDICT_VICAT_ORFS {

    tag "${prefix}"
    conda "${projectDir}/envs/vicat.yml"
    cpus params.vicat_orf_threads
    memory params.vicat_orf_memory
    time params.vicat_orf_time

    publishDir { "${params.outdir}/${prefix}_results/vicat" },
        mode: 'copy',
        pattern: '*.vicat_*'

    input:
    tuple val(prefix), val(type), path(normalized_fasta), path(header_map)

    output:
    tuple val(prefix), val(type),
          path("${prefix}.vicat_orfs.faa"),
          path("${prefix}.vicat_orf_map.tsv"),
          path(header_map),
          path("${prefix}.vicat_orf_run_metadata.tsv"),
          path("${prefix}.vicat_orf_prediction.log"),
          emit: orfs

    script:
    """
    set -euo pipefail

    LOG="${prefix}.vicat_orf_prediction.log"
    : > "\$LOG"

    pyrodigal-gv -p meta -i "${normalized_fasta}" \
        -a gv.faa -f gff -o gv.gff -q -j "${task.cpus}" \
        >> "\$LOG" 2>&1

    rv_args=()
    if [[ "${type}" == 'rna' ]]; then
        pyrodigal-rv -p meta -i "${normalized_fasta}" \
            -a rv.faa -f gff -o rv.gff -q -j "${task.cpus}" \
            >> "\$LOG" 2>&1
        rv_args=(--rv-proteins rv.faa --rv-gff rv.gff)
    fi

    python "${projectDir}/bin/combine_vicat_orfs.py" \
        --sample-id "${prefix}" \
        --input-type "${type}" \
        --gv-proteins gv.faa \
        --gv-gff gv.gff \
        "\${rv_args[@]}" \
        --output-proteins "${prefix}.vicat_orfs.faa" \
        --output-map "${prefix}.vicat_orf_map.tsv" \
        >> "\$LOG" 2>&1

    GV_ORFS=\$(grep -c '^>' gv.faa || true)
    RV_ORFS=0
    [[ "${type}" == 'rna' ]] && RV_ORFS=\$(grep -c '^>' rv.faa || true)
    COMBINED_ORFS=\$(grep -c '^>' "${prefix}.vicat_orfs.faa" || true)

    printf 'sample_id\tinput_type\tpyrodigal_gv_orfs\tpyrodigal_rv_orfs\tcombined_orfs\tcallers\n' \
        > "${prefix}.vicat_orf_run_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$GV_ORFS" "\$RV_ORFS" "\$COMBINED_ORFS" \
        "\$([[ "${type}" == 'rna' ]] && echo 'pyrodigal-gv,pyrodigal-rv' || echo 'pyrodigal-gv')" \
        >> "${prefix}.vicat_orf_run_metadata.tsv"

    echo "VICAT_ORFS sample=${prefix} type=${type} gv=\$GV_ORFS rv=\$RV_ORFS combined=\$COMBINED_ORFS"
    """
}
