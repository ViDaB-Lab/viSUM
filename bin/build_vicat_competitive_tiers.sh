#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: build_vicat_competitive_tiers.sh \
  --base-database VIRAL_DATABASE_DIR \
  --viral-representatives IMGVR5_UViG_representatives.faa.gz \
  --cellular-proteins cellular_candidates.faa[.gz] \
  --cellular-metadata cellular_candidates.tsv[.gz] \
  --destination-root DATABASE_ROOT --work-root BUILD_ROOT --threads N \
  [--tiers small:250000,medium:750000,large:1500000 --seed 41 \
   --cellular-min-seq-id 0.90 --cellular-coverage 0.80 \
   --cellular-min-protein-length 50]

This benchmark-only wrapper creates nested cellular tiers once, then builds
one labeled viral-plus-cellular viCAT DIAMOND database for each tier.
EOF
    exit 2
}

BASE_DATABASE='' VIRAL_REPRESENTATIVES='' CELLULAR_PROTEINS=''
CELLULAR_METADATA='' DESTINATION_ROOT='' WORK_ROOT='' THREADS=''
TIERS='small:250000,medium:750000,large:1500000' SEED='41'
CELLULAR_MIN_SEQ_ID='0.90' CELLULAR_COVERAGE='0.80'
CELLULAR_MIN_PROTEIN_LENGTH='50'
while [[ $# -gt 0 ]]; do
    case "$1" in
        --base-database) BASE_DATABASE="$2"; shift 2 ;;
        --viral-representatives) VIRAL_REPRESENTATIVES="$2"; shift 2 ;;
        --cellular-proteins) CELLULAR_PROTEINS="$2"; shift 2 ;;
        --cellular-metadata) CELLULAR_METADATA="$2"; shift 2 ;;
        --destination-root) DESTINATION_ROOT="$2"; shift 2 ;;
        --work-root) WORK_ROOT="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --tiers) TIERS="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --cellular-min-seq-id) CELLULAR_MIN_SEQ_ID="$2"; shift 2 ;;
        --cellular-coverage) CELLULAR_COVERAGE="$2"; shift 2 ;;
        --cellular-min-protein-length) CELLULAR_MIN_PROTEIN_LENGTH="$2"; shift 2 ;;
        *) usage ;;
    esac
done
for variable in BASE_DATABASE VIRAL_REPRESENTATIVES CELLULAR_PROTEINS CELLULAR_METADATA DESTINATION_ROOT WORK_ROOT THREADS; do
    [[ -n "${!variable}" ]] || usage
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_INPUTS="$WORK_ROOT/tier_inputs"
mkdir -p "$TIER_INPUTS" "$DESTINATION_ROOT"
python "$SCRIPT_DIR/prepare_vicat_cellular_tiers.py" \
    --cellular-proteins "$CELLULAR_PROTEINS" \
    --cellular-metadata "$CELLULAR_METADATA" \
    --output-dir "$TIER_INPUTS" --tiers "$TIERS" --seed "$SEED"

IFS=',' read -ra TIER_ITEMS <<< "$TIERS"
for item in "${TIER_ITEMS[@]}"; do
    name="${item%%:*}"
    bash "$SCRIPT_DIR/build_vicat_competitive_database.sh" \
        --base-database "$BASE_DATABASE" \
        --viral-representatives "$VIRAL_REPRESENTATIVES" \
        --cellular-proteins "$TIER_INPUTS/cellular_${name}.faa.gz" \
        --cellular-metadata "$TIER_INPUTS/cellular_${name}.metadata.tsv.gz" \
        --destination "$DESTINATION_ROOT/$name" \
        --work-dir "$WORK_ROOT/build_${name}" --threads "$THREADS" \
        --cellular-min-seq-id "$CELLULAR_MIN_SEQ_ID" \
        --cellular-coverage "$CELLULAR_COVERAGE" \
        --cellular-min-protein-length "$CELLULAR_MIN_PROTEIN_LENGTH"
done

cp "$TIER_INPUTS/cellular_tier_summary.tsv" "$DESTINATION_ROOT/"
echo "Built viCAT competitive tiers under: $DESTINATION_ROOT"
