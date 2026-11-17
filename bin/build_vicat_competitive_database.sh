#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: build_vicat_competitive_database.sh \
  --base-database VIRAL_DATABASE_DIR \
  --viral-representatives IMGVR5_UViG_representatives.faa.gz \
  --cellular-proteins cellular_proteins.faa[.gz] \
  --cellular-metadata cellular_proteins.tsv \
  --destination DATABASE_DIR --work-dir BUILD_DIR --threads N \
  [--cellular-min-seq-id 0.90 --cellular-coverage 0.80 \
   --cellular-min-protein-length 50]

The base database must contain the established viCAT viral DIAMOND database
and taxonomy lookup. Cellular metadata must be tab-delimited and contain
protein_id, cellular_group, and source_accession columns.
EOF
    exit 2
}

BASE_DATABASE=''
VIRAL_REPRESENTATIVES=''
CELLULAR_PROTEINS=''
CELLULAR_METADATA=''
DESTINATION=''
WORK_DIR=''
THREADS=''
CELLULAR_MIN_SEQ_ID='0.90'
CELLULAR_COVERAGE='0.80'
CELLULAR_MIN_PROTEIN_LENGTH='50'

while [[ $# -gt 0 ]]; do
    case "$1" in
        --base-database) BASE_DATABASE="$2"; shift 2 ;;
        --viral-representatives) VIRAL_REPRESENTATIVES="$2"; shift 2 ;;
        --cellular-proteins) CELLULAR_PROTEINS="$2"; shift 2 ;;
        --cellular-metadata) CELLULAR_METADATA="$2"; shift 2 ;;
        --destination) DESTINATION="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --cellular-min-seq-id) CELLULAR_MIN_SEQ_ID="$2"; shift 2 ;;
        --cellular-coverage) CELLULAR_COVERAGE="$2"; shift 2 ;;
        --cellular-min-protein-length) CELLULAR_MIN_PROTEIN_LENGTH="$2"; shift 2 ;;
        *) usage ;;
    esac
done

for variable in BASE_DATABASE VIRAL_REPRESENTATIVES CELLULAR_PROTEINS \
    CELLULAR_METADATA DESTINATION WORK_DIR THREADS; do
    [[ -n "${!variable}" ]] || usage
done

for command_name in python mmseqs diamond pigz sha256sum gzip seqkit; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: Required viCAT build command is unavailable: $command_name" >&2
        exit 1
    }
done

locate_viral_runtime() {
    local candidate="$1"
    if [[ -s "$candidate/IMGVR5_UViG_representatives.dmnd" && \
          -s "$candidate/IMGVR5_UViG.vicat_taxonomy_lookup.parquet" ]]; then
        VIRAL_DIAMOND="$candidate/IMGVR5_UViG_representatives.dmnd"
        VIRAL_TAXONOMY="$candidate/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
        VIRAL_SUMMARY="$candidate/vicat_taxonomy_build_summary.tsv"
        return 0
    fi
    if [[ -s "$candidate/diamond/IMGVR5_UViG_representatives.dmnd" && \
          -s "$candidate/taxonomy/IMGVR5_UViG.vicat_taxonomy_lookup.parquet" ]]; then
        VIRAL_DIAMOND="$candidate/diamond/IMGVR5_UViG_representatives.dmnd"
        VIRAL_TAXONOMY="$candidate/taxonomy/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
        VIRAL_SUMMARY="$candidate/taxonomy/IMGVR5_UViG.vicat_taxonomy_build_summary.tsv"
        return 0
    fi
    return 1
}

locate_viral_runtime "$BASE_DATABASE" || {
    echo "ERROR: Base viCAT viral database failed validation: $BASE_DATABASE" >&2
    exit 1
}
diamond dbinfo --db "$VIRAL_DIAMOND" >/dev/null
[[ -f "$VIRAL_REPRESENTATIVES" ]] || {
    echo "ERROR: Viral representative FASTA not found: $VIRAL_REPRESENTATIVES" >&2
    exit 1
}
[[ -f "$CELLULAR_PROTEINS" ]] || {
    echo "ERROR: Cellular protein FASTA not found: $CELLULAR_PROTEINS" >&2
    exit 1
}
[[ -f "$CELLULAR_METADATA" ]] || {
    echo "ERROR: Cellular metadata not found: $CELLULAR_METADATA" >&2
    exit 1
}
[[ ! -e "$DESTINATION" ]] || {
    echo "ERROR: Refusing to replace existing destination: $DESTINATION" >&2
    exit 1
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR/mmseqs" "$WORK_DIR/runtime" "$WORK_DIR/manifest_work"

FILTERED="$WORK_DIR/cellular.filtered.faa"
seqkit seq --remove-gaps --min-len "$CELLULAR_MIN_PROTEIN_LENGTH" \
    "$CELLULAR_PROTEINS" > "$FILTERED"
[[ "$(grep -c '^>' "$FILTERED")" -gt 0 ]] || {
    echo "ERROR: No cellular proteins remain after length filtering." >&2
    exit 1
}

CELL_DB="$WORK_DIR/mmseqs/cellular"
CLUSTER_DB="$WORK_DIR/mmseqs/cellular_clusters"
REP_DB="$WORK_DIR/mmseqs/cellular_representatives"
REP_FASTA="$WORK_DIR/cellular_representatives.faa"
mmseqs createdb "$FILTERED" "$CELL_DB" --dbtype 1 --shuffle 0 --createdb-mode 0
mmseqs linclust "$CELL_DB" "$CLUSTER_DB" "$WORK_DIR/mmseqs/tmp" \
    --min-seq-id "$CELLULAR_MIN_SEQ_ID" --cov-mode 0 \
    -c "$CELLULAR_COVERAGE" --threads "$THREADS"
mmseqs createsubdb "$CLUSTER_DB" "$CELL_DB" "$REP_DB" --subdb-mode 1
mmseqs convert2fasta "$REP_DB" "$REP_FASTA"
pigz -p "$THREADS" "$REP_FASTA"
REP_FASTA="${REP_FASTA}.gz"

COMBINED_FASTA="$WORK_DIR/vicat_viral_cellular.faa"
MANIFEST="$WORK_DIR/runtime/vicat_competitive_reference_manifest.parquet"
SUMMARY="$WORK_DIR/runtime/vicat_competitive_database_summary.tsv"
python "$SCRIPT_DIR/prepare_vicat_competitive_references.py" \
    --viral-representatives "$VIRAL_REPRESENTATIVES" \
    --viral-taxonomy-lookup "$VIRAL_TAXONOMY" \
    --cellular-representatives "$REP_FASTA" \
    --cellular-metadata "$CELLULAR_METADATA" \
    --output-fasta "$COMBINED_FASTA" \
    --output-manifest "$MANIFEST" \
    --output-summary "$SUMMARY" \
    --work-dir "$WORK_DIR/manifest_work"

diamond makedb --in "$COMBINED_FASTA" \
    --db "$WORK_DIR/runtime/vicat_viral_cellular" --threads "$THREADS"
diamond dbinfo --db "$WORK_DIR/runtime/vicat_viral_cellular.dmnd" >/dev/null

# Retain the existing viral-only runtime so current viCAT analysis remains
# backward compatible until competitive ORF scoring is enabled separately.
cp --reflink=auto "$VIRAL_DIAMOND" \
    "$WORK_DIR/runtime/IMGVR5_UViG_representatives.dmnd"
cp --reflink=auto "$VIRAL_TAXONOMY" \
    "$WORK_DIR/runtime/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
if [[ -s "$VIRAL_SUMMARY" ]]; then
    cp "$VIRAL_SUMMARY" "$WORK_DIR/runtime/vicat_taxonomy_build_summary.tsv"
fi

viral_count="$(diamond dbinfo --db "$VIRAL_DIAMOND" | awk '$1 == "Sequences" {print $2}')"
cellular_count="$(awk -F '\t' '$1 == "cellular_representatives" {print $2}' "$SUMMARY")"
combined_count="$(diamond dbinfo --db "$WORK_DIR/runtime/vicat_viral_cellular.dmnd" \
    | awk '$1 == "Sequences" {print $2}')"
expected_combined=$((viral_count + cellular_count))
[[ "$combined_count" == "$expected_combined" ]] || {
    echo "ERROR: Competitive database contains $combined_count sequences; expected $expected_combined" >&2
    exit 1
}

{
    printf 'field\tvalue\n'
    printf 'database_name\tviCAT labeled viral-plus-cellular references\n'
    printf 'base_viral_database\t%s\n' "$BASE_DATABASE"
    printf 'viral_representatives\t%s\n' "$viral_count"
    printf 'cellular_representatives\t%s\n' "$cellular_count"
    printf 'combined_references\t%s\n' "$combined_count"
    printf 'cellular_clustering_min_seq_id\t%s\n' "$CELLULAR_MIN_SEQ_ID"
    printf 'cellular_clustering_bidirectional_coverage\t%s\n' "$CELLULAR_COVERAGE"
    printf 'cellular_minimum_protein_length\t%s\n' "$CELLULAR_MIN_PROTEIN_LENGTH"
    printf 'reference_labels\tVIRAL|,CELLULAR|\n'
    printf 'build_completed\t%s\n' "$(date -Iseconds)"
} > "$WORK_DIR/runtime/vicat_database_metadata.tsv"

(
    cd "$WORK_DIR/runtime"
    checksum_files=(
        IMGVR5_UViG_representatives.dmnd
        IMGVR5_UViG.vicat_taxonomy_lookup.parquet
        vicat_viral_cellular.dmnd
        vicat_competitive_reference_manifest.parquet
        vicat_competitive_database_summary.tsv
        vicat_database_metadata.tsv
    )
    [[ -s vicat_taxonomy_build_summary.tsv ]] \
        && checksum_files+=(vicat_taxonomy_build_summary.tsv)
    sha256sum "${checksum_files[@]}" > SHA256SUMS
    sha256sum --check SHA256SUMS
    date -Iseconds > .visum_db_complete
)

mkdir -p "$(dirname "$DESTINATION")"
mv "$WORK_DIR/runtime" "$DESTINATION"
rm -rf "$WORK_DIR"
echo "viCAT competitive database build completed: $DESTINATION"
