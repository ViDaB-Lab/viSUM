#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: build_vicat_database.sh \
  --proteins IMGVR5_UViG.faa.gz --metadata IMGVR5_UViG.tsv.gz \
  --destination DATABASE_DIR --work-dir BUILD_DIR \
  --threads N --memory-limit 300GB \
  --protein-sha256 HASH --metadata-sha256 HASH \
  --expected-uvigs N --expected-proteins N --expected-representatives N
EOF
    exit 2
}

PROTEINS=''
METADATA=''
DESTINATION=''
WORK_DIR=''
THREADS=''
MEMORY_LIMIT=''
PROTEIN_SHA256=''
METADATA_SHA256=''
EXPECTED_UVIGS=''
EXPECTED_PROTEINS=''
EXPECTED_REPRESENTATIVES=''

while [[ $# -gt 0 ]]; do
    case "$1" in
        --proteins) PROTEINS="$2"; shift 2 ;;
        --metadata) METADATA="$2"; shift 2 ;;
        --destination) DESTINATION="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --memory-limit) MEMORY_LIMIT="$2"; shift 2 ;;
        --protein-sha256) PROTEIN_SHA256="$2"; shift 2 ;;
        --metadata-sha256) METADATA_SHA256="$2"; shift 2 ;;
        --expected-uvigs) EXPECTED_UVIGS="$2"; shift 2 ;;
        --expected-proteins) EXPECTED_PROTEINS="$2"; shift 2 ;;
        --expected-representatives) EXPECTED_REPRESENTATIVES="$2"; shift 2 ;;
        *) usage ;;
    esac
done

for value in PROTEINS METADATA DESTINATION WORK_DIR THREADS MEMORY_LIMIT \
    PROTEIN_SHA256 METADATA_SHA256 EXPECTED_UVIGS EXPECTED_PROTEINS \
    EXPECTED_REPRESENTATIVES; do
    [[ -n "${!value}" ]] || usage
done

for command_name in python mmseqs diamond pigz sha256sum gzip; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: Required viCAT build command is unavailable: $command_name" >&2
        exit 1
    }
done

[[ -f "$PROTEINS" ]] || { echo "ERROR: Protein source not found: $PROTEINS" >&2; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: Metadata source not found: $METADATA" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$WORK_DIR" "$WORK_DIR/checkpoints" "$WORK_DIR/mmseqs" \
    "$WORK_DIR/dedup" "$WORK_DIR/taxonomy" "$WORK_DIR/taxonomy_work" \
    "$WORK_DIR/runtime"

checkpoint_done() { [[ -f "$WORK_DIR/checkpoints/$1.complete" ]]; }
mark_complete() { date -Iseconds > "$WORK_DIR/checkpoints/$1.complete"; }

if ! checkpoint_done source_validation; then
    gzip -t "$PROTEINS"
    gzip -t "$METADATA"
    printf '%s  %s\n' "$PROTEIN_SHA256" "$PROTEINS" | sha256sum --check -
    printf '%s  %s\n' "$METADATA_SHA256" "$METADATA" | sha256sum --check -
    mark_complete source_validation
fi

COMPACT_METADATA="$WORK_DIR/dedup/IMGVR5_UViG.vicat_metadata.tsv.gz"
if ! checkpoint_done compact_metadata; then
    rm -f "$COMPACT_METADATA"
    python "$SCRIPT_DIR/prepare_vicat_metadata.py" \
        --input "$METADATA" \
        --output "$COMPACT_METADATA" \
        --expected-rows "$EXPECTED_UVIGS"
    mark_complete compact_metadata
fi

SEQ_DB="$WORK_DIR/mmseqs/IMGVR5_UViG"
if ! checkpoint_done mmseqs_createdb; then
    rm -f "${SEQ_DB}"*
    mmseqs createdb "$PROTEINS" "$SEQ_DB" \
        --dbtype 1 --shuffle 0 --createdb-mode 0 --write-lookup 1
    mmseqs dbtype "$SEQ_DB" | grep -q 'Aminoacid'
    mark_complete mmseqs_createdb
fi

HASH_DB="$WORK_DIR/mmseqs/IMGVR5_UViG_clust100"
if ! checkpoint_done mmseqs_clusthash; then
    rm -f "${HASH_DB}"*
    mmseqs clusthash "$SEQ_DB" "$HASH_DB" --min-seq-id 1.0 --threads "$THREADS"
    mark_complete mmseqs_clusthash
fi

CLUSTER_DB="$WORK_DIR/mmseqs/IMGVR5_UViG_clusters100"
if ! checkpoint_done mmseqs_clusters; then
    rm -f "${CLUSTER_DB}"*
    mmseqs clust "$SEQ_DB" "$HASH_DB" "$CLUSTER_DB" \
        --cluster-mode 1 --threads "$THREADS"
    mmseqs dbtype "$CLUSTER_DB" | grep -q 'Clustering'
    mark_complete mmseqs_clusters
fi

MEMBERS="$WORK_DIR/dedup/IMGVR5_UViG_clusters100.members.tsv.gz"
if ! checkpoint_done cluster_members; then
    rm -f "$MEMBERS"
    mmseqs createtsv "$SEQ_DB" "$SEQ_DB" "$CLUSTER_DB" /dev/stdout \
        | pigz -p "$THREADS" > "$MEMBERS"
    gzip -t "$MEMBERS"
    mark_complete cluster_members
fi

REP_DB="$WORK_DIR/mmseqs/IMGVR5_UViG_representatives"
REP_FASTA="$WORK_DIR/dedup/IMGVR5_UViG_representatives.faa.gz"
if ! checkpoint_done representative_fasta; then
    rm -f "${REP_DB}"* "$REP_FASTA" "${REP_FASTA%.gz}"
    mmseqs createsubdb "$CLUSTER_DB" "$SEQ_DB" "$REP_DB" --subdb-mode 1
    mmseqs convert2fasta "$REP_DB" "${REP_FASTA%.gz}"
    pigz -p "$THREADS" "${REP_FASTA%.gz}"
    gzip -t "$REP_FASTA"
    mark_complete representative_fasta
fi

TAXONOMY_LOOKUP="$WORK_DIR/taxonomy/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
TAXONOMY_SUMMARY="$WORK_DIR/taxonomy/IMGVR5_UViG.vicat_taxonomy_build_summary.tsv"
if ! checkpoint_done taxonomy_lookup; then
    rm -f "$WORK_DIR/taxonomy/IMGVR5_UViG.vicat_taxonomy_"*.parquet \
        "$TAXONOMY_SUMMARY"
    taxonomy_resume=()
    [[ -f "$WORK_DIR/taxonomy_work/IMGVR5_UViG.vicat_taxonomy_build.duckdb" ]] \
        && taxonomy_resume=(--resume)
    python "$SCRIPT_DIR/build_vicat_taxonomy_lookup.py" \
        --cluster-members "$MEMBERS" \
        --metadata "$COMPACT_METADATA" \
        --output-dir "$WORK_DIR/taxonomy" \
        --work-dir "$WORK_DIR/taxonomy_work" \
        --threads "$THREADS" \
        --memory-limit "$MEMORY_LIMIT" \
        --prefix IMGVR5_UViG \
        --expected-member-count "$EXPECTED_PROTEINS" \
        --expected-representative-count "$EXPECTED_REPRESENTATIVES" \
        "${taxonomy_resume[@]}"
    mark_complete taxonomy_lookup
fi

DIAMOND_DB="$WORK_DIR/runtime/IMGVR5_UViG_representatives.dmnd"
if ! checkpoint_done diamond_database; then
    rm -f "$DIAMOND_DB"
    diamond makedb --in "$REP_FASTA" \
        --db "${DIAMOND_DB%.dmnd}" --threads "$THREADS"
    diamond dbinfo --db "$DIAMOND_DB" >/dev/null
    mark_complete diamond_database
fi

cp "$TAXONOMY_LOOKUP" "$WORK_DIR/runtime/IMGVR5_UViG.vicat_taxonomy_lookup.parquet"
cp "$TAXONOMY_SUMMARY" "$WORK_DIR/runtime/vicat_taxonomy_build_summary.tsv"

diamond_count="$(diamond dbinfo --db "$DIAMOND_DB" | awk '$1 == "Sequences" {print $2}')"
lookup_count="$(python - "$WORK_DIR/runtime/IMGVR5_UViG.vicat_taxonomy_lookup.parquet" <<'PY'
import duckdb
import sys
connection = duckdb.connect()
print(connection.execute(
    "SELECT count(*) FROM read_parquet(?)", [sys.argv[1]]
).fetchone()[0])
PY
)"
[[ "$diamond_count" == "$EXPECTED_REPRESENTATIVES" ]] || {
    echo "ERROR: DIAMOND sequence count is $diamond_count; expected $EXPECTED_REPRESENTATIVES" >&2
    exit 1
}
[[ "$lookup_count" == "$EXPECTED_REPRESENTATIVES" ]] || {
    echo "ERROR: Taxonomy lookup count is $lookup_count; expected $EXPECTED_REPRESENTATIVES" >&2
    exit 1
}

{
    printf 'field\tvalue\n'
    printf 'database_name\tviCAT MetaVR5 exact-protein representatives\n'
    printf 'metavr_release\tIMGVR5\n'
    printf 'source_protein_sha256\t%s\n' "$PROTEIN_SHA256"
    printf 'source_metadata_sha256\t%s\n' "$METADATA_SHA256"
    printf 'source_uvig_count\t%s\n' "$EXPECTED_UVIGS"
    printf 'source_protein_count\t%s\n' "$EXPECTED_PROTEINS"
    printf 'representative_protein_count\t%s\n' "$EXPECTED_REPRESENTATIVES"
    printf 'deduplication\tMMseqs2 exact identity; connected components\n'
    printf 'taxonomy_aggregation\tvOTU-balanced conservative LCA\n'
    printf 'mmseqs_version\t%s\n' "$(mmseqs version)"
    printf 'diamond_version\t%s\n' "$(diamond version | awk '{print $3}')"
    printf 'build_completed\t%s\n' "$(date -Iseconds)"
} > "$WORK_DIR/runtime/vicat_database_metadata.tsv"

(
    cd "$WORK_DIR/runtime"
    sha256sum IMGVR5_UViG_representatives.dmnd \
        IMGVR5_UViG.vicat_taxonomy_lookup.parquet \
        vicat_taxonomy_build_summary.tsv \
        vicat_database_metadata.tsv > SHA256SUMS
    sha256sum --check SHA256SUMS
    date -Iseconds > .visum_db_complete
)

if [[ -e "$DESTINATION" ]]; then
    echo "ERROR: Refusing to replace existing destination: $DESTINATION" >&2
    exit 1
fi
mkdir -p "$(dirname "$DESTINATION")"
mv "$WORK_DIR/runtime" "$DESTINATION"

# Only validated runtime files survive a successful build.
rm -rf "$WORK_DIR"
echo "viCAT database build completed: $DESTINATION"
