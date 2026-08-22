process RUN_TESORTER {

    tag "${prefix}"

    conda "${projectDir}/envs/tesorter.yml"

    cpus { Math.min(params.tesorter_cpus as int, params.max_cpus as int) }
    memory params.tesorter_memory
    time params.tesorter_time

    publishDir { "${params.outdir}/${prefix}_results/tesorter" },
        mode: 'copy'

    input:
    tuple val(prefix), val(type), path(refined_fasta), path(region_map),
          path(boundary_audit), path(refinement_summary)

    output:
    tuple val(prefix),
          val(type),
          path(region_map),
          path("${prefix}.tesorter.rexdb.cls.tsv"),
          path("${prefix}.tesorter.rexdb.dom.tsv"),
          path("${prefix}.tesorter.rexdb.dom.gff3"),
          path("${prefix}.tesorter.log"),
          path("${prefix}.tesorter_run_metadata.tsv"),
          emit: results

    script:
    """
    set -euo pipefail

    OUTPUT_PREFIX="${prefix}.tesorter"
    CLASSIFICATIONS="${prefix}.tesorter.rexdb.cls.tsv"
    DOMAINS="${prefix}.tesorter.rexdb.dom.tsv"
    DOMAIN_GFF="${prefix}.tesorter.rexdb.dom.gff3"
    LOG_FILE="${prefix}.tesorter.log"
    METADATA_FILE="${prefix}.tesorter_run_metadata.tsv"

    INPUT_COUNT=\$(grep -c '^>' "${refined_fasta}" || true)
    TESORTER_VERSION=\$(TEsorter --version 2>&1 | tail -n 1 | tr -d '\r')

    if [[ "\$INPUT_COUNT" -eq 0 ]]; then
        : > "\$CLASSIFICATIONS"
        : > "\$DOMAINS"
        : > "\$DOMAIN_GFF"
        printf 'TEsorter skipped: refined candidate FASTA was empty.\n' > "\$LOG_FILE"
        RUN_STATUS='skipped_empty_input'
    else
        if ! TEsorter "${refined_fasta}" \
            -db rexdb \
            -pre "\$OUTPUT_PREFIX" \
            -p "${task.cpus}" \
            > "\$LOG_FILE" 2>&1; then
            echo "ERROR: TEsorter failed for sample '${prefix}'." >&2
            tail -n 50 "\$LOG_FILE" >&2 || true
            exit 1
        fi

        # Successful no-hit runs may leave one or more report files empty or
        # absent. Preserve that valid outcome as explicit empty evidence.
        [[ -f "\$CLASSIFICATIONS" ]] || : > "\$CLASSIFICATIONS"
        [[ -f "\$DOMAINS" ]] || : > "\$DOMAINS"
        [[ -f "\$DOMAIN_GFF" ]] || : > "\$DOMAIN_GFF"

        CLASSIFICATION_COUNT=\$(awk 'NF { count++ } END { print count + 0 }' "\$CLASSIFICATIONS")
        if [[ "\$CLASSIFICATION_COUNT" -eq 0 ]]; then
            RUN_STATUS='completed_no_te_classifications'
        else
            RUN_STATUS='completed_with_te_classifications'
        fi
    fi

    CLASSIFICATION_COUNT=\$(awk 'NF { count++ } END { print count + 0 }' "\$CLASSIFICATIONS")
    DOMAIN_COUNT=\$(awk 'NF && \$1 !~ /^#/ { count++ } END { print count + 0 }' "\$DOMAINS")

    printf 'sample_id\tinput_type\tinput_sequence_count\tte_classification_count\tdomain_row_count\ttesorter_version\tdatabase\tthreads\trun_status\tinput_fasta\n' > "\$METADATA_FILE"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${prefix}" "${type}" "\$INPUT_COUNT" "\$CLASSIFICATION_COUNT" \
        "\$DOMAIN_COUNT" "\$TESORTER_VERSION" 'rexdb' "${task.cpus}" \
        "\$RUN_STATUS" "${refined_fasta}" >> "\$METADATA_FILE"

    echo "TESORTER sample=${prefix} type=${type} classifications=\$CLASSIFICATION_COUNT domains=\$DOMAIN_COUNT"
    """
}

