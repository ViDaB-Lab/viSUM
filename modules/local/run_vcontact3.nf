process RUN_VCONTACT3 {

    tag "${prefix}"

    conda "${projectDir}/envs/vcontact3.yml"

    cpus { Math.min(params.vcontact3_cpus as int, params.max_cpus as int) }
    memory { params.vcontact3_memory }
    time { params.vcontact3_time }

    publishDir { "${params.outdir}/${prefix}_results/vcontact3" },
        mode: 'copy',
        pattern: '*.vcontact3_*.csv'
    publishDir { "${params.outdir}/${prefix}_results/vcontact3" },
        mode: 'copy',
        pattern: '*.vcontact3_*.tsv'
    publishDir { "${params.outdir}/${prefix}_results/vcontact3" },
        mode: 'copy',
        pattern: '*.vcontact3_*.log'

    input:
    tuple val(prefix),
          val(type),
          path(refined_fasta),
          path(provirus_region_map),
          path(provirus_boundary_audit),
          path(provirus_refinement_summary)
    tuple path(vcontact3_database), path(vcontact3_database_metadata)
    val database_domain_mode

    output:
    tuple val(prefix),
          val(type),
          path(provirus_region_map),
          path("${prefix}.vcontact3_*_final_assignments.csv"),
          path("${prefix}.vcontact3_*_performance_metrics.csv"),
          path("${prefix}.vcontact3_*.log"),
          path("${prefix}.vcontact3_run_metadata.tsv"),
          emit: results
    path("${prefix}.vcontact3_raw"), emit: raw_output

    script:
    """
    set -euo pipefail

    RAW_DIRECTORY="${prefix}.vcontact3_raw"
    METADATA_FILE="${prefix}.vcontact3_run_metadata.tsv"
    mkdir -p "\$RAW_DIRECTORY"

    case "${database_domain_mode}" in
        both)
            DOMAINS=(prokaryotes eukaryotes)
            ;;
        prokaryotes|eukaryotes)
            DOMAINS=("${database_domain_mode}")
            ;;
        *)
            echo "ERROR: Unsupported vConTACT3 database domain mode: ${database_domain_mode}" >&2
            exit 1
            ;;
    esac

    INPUT_COUNT=\$(awk '/^>/ { count++ } END { print count + 0 }' "${refined_fasta}")
    DATABASE_VERSION=\$(awk -F '\t' 'NR == 2 { print \$6 }' "${vcontact3_database_metadata}")
    RESOLVED_DATABASE_PATH=\$(readlink -f "${vcontact3_database}")
    VCONTACT3_VERSION=\$(python -c "from importlib.metadata import version; print(version('vcontact3'))")

    if [[ ! "\$DATABASE_VERSION" =~ ^[0-9]{3}\$ ]]; then
        echo "ERROR: Could not read a three-digit installed vConTACT3 database version from ${vcontact3_database_metadata}." >&2
        exit 1
    fi
    if [[ -z "\$RESOLVED_DATABASE_PATH" || ! -d "\$RESOLVED_DATABASE_PATH" ]]; then
        echo "ERROR: Could not resolve the vConTACT3 database path from ${vcontact3_database}." >&2
        exit 1
    fi

    printf 'sample_id\tinput_type\tdatabase_domain\trefined_sequence_count\tassignment_row_count\tuser_assignment_row_count\tthreads\tvcontact3_version\tdatabase_path\tdatabase_version\trun_status\tinput_fasta\n' \
        > "\$METADATA_FILE"

    for DOMAIN in "\${DOMAINS[@]}"; do
        DOMAIN_OUTPUT="\$RAW_DIRECTORY/\$DOMAIN"
        ASSIGNMENTS="${prefix}.vcontact3_\${DOMAIN}_final_assignments.csv"
        METRICS="${prefix}.vcontact3_\${DOMAIN}_performance_metrics.csv"
        LOG_FILE="${prefix}.vcontact3_\${DOMAIN}.log"
        mkdir -p "\$DOMAIN_OUTPUT"

        if [[ "\$INPUT_COUNT" -eq 0 ]]; then
            mkdir -p "\$DOMAIN_OUTPUT/exports"
            printf 'Genome,GenomeName,Reference,Size_Kb\n' \
                > "\$DOMAIN_OUTPUT/exports/final_assignments.csv"
            printf 'index,realm,rank\n' \
                > "\$DOMAIN_OUTPUT/exports/performance_metrics.csv"
            printf 'vConTACT3 skipped: provirus refinement produced no candidate sequences.\n' \
                > "\$LOG_FILE"
            RUN_STATUS='skipped_no_refined_candidates'
        else
            if ! vcontact3 run \
                --nucleotide "${refined_fasta}" \
                --output "\$DOMAIN_OUTPUT" \
                --db-path "\$RESOLVED_DATABASE_PATH" \
                --db-version "\$DATABASE_VERSION" \
                --db-domain "\$DOMAIN" \
                --pyrodigal-gv \
                --threads "${task.cpus}" \
                --no-progress \
                > "\$LOG_FILE" 2>&1; then
                echo "ERROR: vConTACT3 \$DOMAIN analysis failed for sample '${prefix}'." >&2
                echo "Last 50 lines of \$LOG_FILE:" >&2
                tail -n 50 "\$LOG_FILE" >&2 || true
                exit 1
            fi
            RUN_STATUS='completed'
        fi

        for REQUIRED_FILE in final_assignments.csv performance_metrics.csv; do
            if [[ ! -s "\$DOMAIN_OUTPUT/exports/\$REQUIRED_FILE" ]]; then
                echo "ERROR: vConTACT3 \$DOMAIN run did not produce \$REQUIRED_FILE." >&2
                exit 1
            fi
        done

        cp "\$DOMAIN_OUTPUT/exports/final_assignments.csv" "\$ASSIGNMENTS"
        cp "\$DOMAIN_OUTPUT/exports/performance_metrics.csv" "\$METRICS"

        ASSIGNMENT_COUNT=\$(awk 'NR > 1 { count++ } END { print count + 0 }' "\$ASSIGNMENTS")
        USER_ASSIGNMENT_COUNT=\$(python -c "import csv,sys; rows=csv.DictReader(open(sys.argv[1], newline='', encoding='utf-8-sig')); print(sum(str(row.get('Reference', '')).strip().lower() in {'false','0'} for row in rows))" "\$ASSIGNMENTS")

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${prefix}" "${type}" "\$DOMAIN" "\$INPUT_COUNT" \
            "\$ASSIGNMENT_COUNT" "\$USER_ASSIGNMENT_COUNT" "${task.cpus}" \
            "\$VCONTACT3_VERSION" "\$RESOLVED_DATABASE_PATH" \
            "\$DATABASE_VERSION" "\$RUN_STATUS" "${refined_fasta.name}" \
            >> "\$METADATA_FILE"
    done

    echo "VCONTACT3 sample=${prefix} type=${type} domains=${database_domain_mode} refined=\$INPUT_COUNT"
    """
}
