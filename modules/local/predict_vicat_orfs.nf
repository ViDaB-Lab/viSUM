process PREDICT_VICAT_ORFS {

    tag "${prefix}"
    conda "${projectDir}/envs/vicat.yml"
    cpus { Math.min(params.vicat_orf_cpus as int, params.max_cpus as int) }
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
          emit: orfs
    path("${prefix}.vicat_orf_run_metadata.tsv"), emit: metadata
    path("${prefix}.vicat_orf_prediction.log"), emit: log

    script:
    """
    set -euo pipefail

    export OMP_NUM_THREADS="${task.cpus}"
    export MKL_NUM_THREADS="${task.cpus}"
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export NUMEXPR_NUM_THREADS="${task.cpus}"

    LOG="${prefix}.vicat_orf_prediction.log"
    : > "\$LOG"

    python "${projectDir}/bin/predict_vicat_orfs_batched.py" --caller gv --input "${normalized_fasta}" \
        --proteins gv.faa --gff gv.gff --threads "${task.cpus}" \
        >> "\$LOG" 2>&1

    rv_args=()
    if [[ "${type}" == 'rna' ]]; then
        python "${projectDir}/bin/predict_vicat_orfs_batched.py" --caller rv --input "${normalized_fasta}" \
            --proteins rv.faa --gff rv.gff --threads "${task.cpus}" \
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

    printf 'sample_id\tinput_type\tpyrodigal_gv_orfs\tpyrodigal_rv_orfs\tcombined_orfs\tcallers\tthreads\n' \
        > "${prefix}.vicat_orf_run_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$GV_ORFS" "\$RV_ORFS" "\$COMBINED_ORFS" \
        "\$([[ "${type}" == 'rna' ]] && echo 'pyrodigal-gv,pyrodigal-rv' || echo 'pyrodigal-gv')" \
        "${task.cpus}" \
        >> "${prefix}.vicat_orf_run_metadata.tsv"

    echo "VICAT_ORFS sample=${prefix} type=${type} gv=\$GV_ORFS rv=\$RV_ORFS combined=\$COMBINED_ORFS"
    """
}
