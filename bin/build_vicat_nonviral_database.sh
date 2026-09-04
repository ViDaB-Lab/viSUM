#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: build_vicat_nonviral_database.sh \
  --classified-dir NONVIRAL_CLASSIFIED_DIR \
  --destination NONVIRAL_DATABASE_DIR --work-dir BUILD_DIR \
  --threads N [--min-seq-id 0.90 --coverage 0.80]

Clusters each viCAT nonviral reference class independently, preserves class
labels and cluster membership, and builds a standalone nonviral DIAMOND
database. An existing work directory is resumed only when its input hashes,
parameters, and builder scripts match the current invocation.
EOF
    exit 2
}

CLASSIFIED_DIR=''
DESTINATION=''
WORK_DIR=''
THREADS=''
MIN_SEQ_ID='0.90'
COVERAGE='0.80'

while [[ $# -gt 0 ]]; do
    case "$1" in
        --classified-dir) CLASSIFIED_DIR="$2"; shift 2 ;;
        --destination) DESTINATION="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --min-seq-id) MIN_SEQ_ID="$2"; shift 2 ;;
        --coverage) COVERAGE="$2"; shift 2 ;;
        *) usage ;;
    esac
done

for variable in CLASSIFIED_DIR DESTINATION WORK_DIR THREADS; do
    [[ -n "${!variable}" ]] || usage
done
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: --threads must be a positive integer" >&2
    exit 2
}
[[ ! -e "$DESTINATION" ]] || {
    echo "ERROR: Refusing to replace existing destination: $DESTINATION" >&2
    exit 1
}

for command_name in python mmseqs diamond pigz sha256sum gzip awk grep cmp; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: Required viCAT build command is unavailable: $command_name" >&2
        exit 1
    }
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METADATA="$CLASSIFIED_DIR/vicat_nonviral_reference_metadata.tsv.gz"
CLASSIFICATION_SUMMARY="$CLASSIFIED_DIR/vicat_nonviral_classification_summary.tsv"
CLASSES=(
    cellular_chromosome
    cellular_unplaced
    plasmid
    plastid
    mitochondrial
    shared_nonviral
)

[[ -s "$METADATA" ]] || {
    echo "ERROR: Classified metadata not found: $METADATA" >&2
    exit 1
}
[[ -s "$CLASSIFICATION_SUMMARY" ]] || {
    echo "ERROR: Classification summary not found: $CLASSIFICATION_SUMMARY" >&2
    exit 1
}
gzip -t "$METADATA"
for reference_class in "${CLASSES[@]}"; do
    source_fasta="$CLASSIFIED_DIR/${reference_class}.faa.gz"
    [[ -s "$source_fasta" ]] || {
        echo "ERROR: Classified FASTA not found: $source_fasta" >&2
        exit 1
    }
    gzip -t "$source_fasta"
done

mkdir -p \
    "$WORK_DIR/checkpoints" \
    "$WORK_DIR/config" \
    "$WORK_DIR/input" \
    "$WORK_DIR/mmseqs" \
    "$WORK_DIR/membership" \
    "$WORK_DIR/representatives" \
    "$WORK_DIR/package"

CURRENT_CONFIG="$WORK_DIR/config/build_config.current.tsv"
SAVED_CONFIG="$WORK_DIR/config/build_config.tsv"
rm -f "$CURRENT_CONFIG"
{
    printf 'parameter\tvalue\n'
    printf 'min_seq_id\t%s\n' "$MIN_SEQ_ID"
    printf 'coverage\t%s\n' "$COVERAGE"
    printf 'class_count\t%s\n' "${#CLASSES[@]}"
    printf 'builder_sha256\t%s\n' "$(sha256sum "$0" | awk '{print $1}')"
    printf 'metadata_helper_sha256\t%s\n' \
        "$(sha256sum "$SCRIPT_DIR/prepare_vicat_nonviral_cluster_metadata.py" | awk '{print $1}')"
    printf 'metadata_sha256\t%s\n' "$(sha256sum "$METADATA" | awk '{print $1}')"
    for reference_class in "${CLASSES[@]}"; do
        source_fasta="$CLASSIFIED_DIR/${reference_class}.faa.gz"
        printf '%s_sha256\t%s\n' "$reference_class" \
            "$(sha256sum "$source_fasta" | awk '{print $1}')"
    done
} > "$CURRENT_CONFIG"

if [[ -s "$SAVED_CONFIG" ]]; then
    cmp -s "$CURRENT_CONFIG" "$SAVED_CONFIG" || {
        echo "ERROR: Existing work directory belongs to different inputs, parameters, or builder code." >&2
        echo "Use a new --work-dir or remove the old work directory after verifying its path." >&2
        exit 1
    }
    rm -f "$CURRENT_CONFIG"
else
    mv "$CURRENT_CONFIG" "$SAVED_CONFIG"
fi

checkpoint_done() {
    [[ -s "$WORK_DIR/checkpoints/$1.complete" ]]
}

mark_complete() {
    date -Iseconds > "$WORK_DIR/checkpoints/$1.complete"
}

for reference_class in "${CLASSES[@]}"; do
    class_label="${reference_class^^}"
    source_fasta="$CLASSIFIED_DIR/${reference_class}.faa.gz"
    input_fasta="$WORK_DIR/input/${reference_class}.faa"
    sequence_db="$WORK_DIR/mmseqs/${reference_class}"
    cluster_db="$WORK_DIR/mmseqs/${reference_class}_clusters"
    representative_db="$WORK_DIR/mmseqs/${reference_class}_representatives"
    class_tmp="$WORK_DIR/mmseqs/${reference_class}_tmp"
    membership="$WORK_DIR/membership/${reference_class}.tsv"
    representative_fasta="$WORK_DIR/representatives/${reference_class}.faa.gz"

    if ! checkpoint_done "${reference_class}_unpacked"; then
        echo "UNPACK class=$class_label"
        rm -f "$input_fasta" "$input_fasta.tmp"
        gzip -cd "$source_fasta" > "$input_fasta.tmp"
        [[ "$(grep -c '^>' "$input_fasta.tmp")" -gt 0 ]] || {
            echo "ERROR: No sequences found for $class_label" >&2
            exit 1
        }
        mv "$input_fasta.tmp" "$input_fasta"
        mark_complete "${reference_class}_unpacked"
    fi

    if ! checkpoint_done "${reference_class}_createdb"; then
        echo "CREATE_DB class=$class_label"
        rm -f "${sequence_db}"*
        mmseqs createdb "$input_fasta" "$sequence_db" \
            --dbtype 1 --shuffle 0 --createdb-mode 0
        mmseqs dbtype "$sequence_db" | grep -q 'Aminoacid'
        mark_complete "${reference_class}_createdb"
    fi

    if ! checkpoint_done "${reference_class}_clustered"; then
        echo "CLUSTER class=$class_label min_seq_id=$MIN_SEQ_ID coverage=$COVERAGE"
        rm -f "${cluster_db}"* "${representative_db}"* "$membership" \
            "${representative_fasta%.gz}" "$representative_fasta"
        rm -rf "$class_tmp"
        mmseqs linclust "$sequence_db" "$cluster_db" "$class_tmp" \
            --min-seq-id "$MIN_SEQ_ID" --cov-mode 0 -c "$COVERAGE" \
            --threads "$THREADS"
        mmseqs createtsv "$sequence_db" "$sequence_db" "$cluster_db" "$membership"
        mmseqs createsubdb "$cluster_db" "$sequence_db" "$representative_db" \
            --subdb-mode 1
        mmseqs convert2fasta "$representative_db" "${representative_fasta%.gz}"
        pigz -f -p "$THREADS" "${representative_fasta%.gz}"
        gzip -t "$representative_fasta"

        input_count="$(grep -c '^>' "$input_fasta")"
        member_count="$(awk 'END {print NR}' "$membership")"
        representative_count="$(gzip -cd "$representative_fasta" | grep -c '^>')"
        [[ "$member_count" == "$input_count" ]] || {
            echo "ERROR: $class_label membership count $member_count != input count $input_count" >&2
            exit 1
        }
        [[ "$representative_count" -gt 0 && "$representative_count" -le "$input_count" ]] || {
            echo "ERROR: Invalid representative count for $class_label: $representative_count" >&2
            exit 1
        }
        echo "CLUSTERED class=$class_label input=$input_count representatives=$representative_count"
        mark_complete "${reference_class}_clustered"
    else
        echo "RESUME class=$class_label status=clustered"
    fi
done

COMBINED_MEMBERSHIP="$WORK_DIR/package/vicat_nonviral_cluster_membership.tsv.gz"
COMBINED_FASTA="$WORK_DIR/package/vicat_nonviral_representatives.faa.gz"
if ! checkpoint_done combined_references; then
    membership_tmp="$WORK_DIR/package/vicat_nonviral_cluster_membership.tsv.tmp"
    fasta_tmp="$WORK_DIR/package/vicat_nonviral_representatives.faa.tmp"
    rm -f "$membership_tmp" "$COMBINED_MEMBERSHIP" "$fasta_tmp" "$COMBINED_FASTA"
    printf 'reference_class\trepresentative_id\tmember_id\n' > "$membership_tmp"
    : > "$fasta_tmp"
    for reference_class in "${CLASSES[@]}"; do
        class_label="${reference_class^^}"
        awk -v class_label="$class_label" 'BEGIN {OFS="\t"} {print class_label, $1, $2}' \
            "$WORK_DIR/membership/${reference_class}.tsv" >> "$membership_tmp"
        gzip -cd "$WORK_DIR/representatives/${reference_class}.faa.gz" >> "$fasta_tmp"
    done
    pigz -p "$THREADS" -c "$membership_tmp" > "$COMBINED_MEMBERSHIP"
    pigz -p "$THREADS" -c "$fasta_tmp" > "$COMBINED_FASTA"
    rm -f "$membership_tmp" "$fasta_tmp"
    gzip -t "$COMBINED_MEMBERSHIP" "$COMBINED_FASTA"
    mark_complete combined_references
fi

REPRESENTATIVE_METADATA_TSV="$WORK_DIR/package/vicat_nonviral_representative_metadata.tsv.gz"
REPRESENTATIVE_METADATA_PARQUET="$WORK_DIR/package/vicat_nonviral_representative_metadata.parquet"
CLUSTER_SUMMARY="$WORK_DIR/package/vicat_nonviral_cluster_summary.tsv"
if ! checkpoint_done representative_metadata; then
    rm -f "$REPRESENTATIVE_METADATA_TSV" "$REPRESENTATIVE_METADATA_PARQUET" \
        "$CLUSTER_SUMMARY" "$WORK_DIR/metadata.duckdb"
    python "$SCRIPT_DIR/prepare_vicat_nonviral_cluster_metadata.py" \
        --metadata "$METADATA" \
        --membership "$COMBINED_MEMBERSHIP" \
        --output-tsv "$REPRESENTATIVE_METADATA_TSV" \
        --output-parquet "$REPRESENTATIVE_METADATA_PARQUET" \
        --output-summary "$CLUSTER_SUMMARY" \
        --work-database "$WORK_DIR/metadata.duckdb" \
        --threads "$THREADS"
    gzip -t "$REPRESENTATIVE_METADATA_TSV"
    rm -f "$WORK_DIR/metadata.duckdb"
    mark_complete representative_metadata
fi

NONVIRAL_DIAMOND="$WORK_DIR/package/vicat_nonviral.dmnd"
if ! checkpoint_done diamond_database; then
    rm -f "$NONVIRAL_DIAMOND"
    diamond makedb --in "$COMBINED_FASTA" \
        --db "${NONVIRAL_DIAMOND%.dmnd}" --threads "$THREADS"
    diamond dbinfo --db "$NONVIRAL_DIAMOND" >/dev/null
    mark_complete diamond_database
fi

representative_count="$(gzip -cd "$COMBINED_FASTA" | grep -c '^>')"
diamond_count="$(diamond dbinfo --db "$NONVIRAL_DIAMOND" | awk '$1 == "Sequences" {print $2}')"
metadata_count="$(python - "$REPRESENTATIVE_METADATA_PARQUET" <<'PY'
import duckdb
import sys
connection = duckdb.connect()
print(connection.execute("SELECT count(*) FROM read_parquet(?)", [sys.argv[1]]).fetchone()[0])
PY
)"
[[ "$diamond_count" == "$representative_count" ]] || {
    echo "ERROR: DIAMOND count $diamond_count != representative FASTA count $representative_count" >&2
    exit 1
}
[[ "$metadata_count" == "$representative_count" ]] || {
    echo "ERROR: Metadata count $metadata_count != representative FASTA count $representative_count" >&2
    exit 1
}

mkdir -p "$WORK_DIR/package/classes"
for reference_class in "${CLASSES[@]}"; do
    cp "$WORK_DIR/representatives/${reference_class}.faa.gz" \
        "$WORK_DIR/package/classes/${reference_class}.representatives.faa.gz"
done
cp "$METADATA" "$WORK_DIR/package/vicat_nonviral_source_metadata.tsv.gz"
cp "$CLASSIFICATION_SUMMARY" \
    "$WORK_DIR/package/vicat_nonviral_classification_summary.tsv"

BUILD_METADATA="$WORK_DIR/package/vicat_nonviral_database_metadata.tsv"
source_protein_count="$(awk -F '\t' 'NR > 1 {total += $2} END {print total + 0}' "$CLUSTER_SUMMARY")"
{
    printf 'field\tvalue\n'
    printf 'database_name\tviCAT class-aware nonviral references\n'
    printf 'source_proteins\t%s\n' "$source_protein_count"
    printf 'representatives\t%s\n' "$representative_count"
    printf 'clustering_min_seq_id\t%s\n' "$MIN_SEQ_ID"
    printf 'clustering_bidirectional_coverage\t%s\n' "$COVERAGE"
    printf 'minimum_protein_length\tnone\n'
    printf 'threads_used\t%s\n' "$THREADS"
    printf 'reference_classes\t%s\n' "$(IFS=,; echo "${CLASSES[*]^^}")"
    printf 'build_completed\t%s\n' "$(date -Iseconds)"
} > "$BUILD_METADATA"

cp "$SAVED_CONFIG" "$WORK_DIR/package/vicat_nonviral_build_config.tsv"
(
    cd "$WORK_DIR/package"
    sha256sum \
        vicat_nonviral.dmnd \
        vicat_nonviral_representatives.faa.gz \
        vicat_nonviral_cluster_membership.tsv.gz \
        vicat_nonviral_representative_metadata.tsv.gz \
        vicat_nonviral_representative_metadata.parquet \
        vicat_nonviral_source_metadata.tsv.gz \
        vicat_nonviral_classification_summary.tsv \
        vicat_nonviral_cluster_summary.tsv \
        vicat_nonviral_database_metadata.tsv \
        vicat_nonviral_build_config.tsv \
        classes/*.representatives.faa.gz > SHA256SUMS
    sha256sum --check SHA256SUMS
    date -Iseconds > .visum_db_complete
)

mkdir -p "$(dirname "$DESTINATION")"
mv "$WORK_DIR/package" "$DESTINATION"
echo "viCAT nonviral database build completed: $DESTINATION"
