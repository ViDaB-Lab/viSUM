#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: build_vicat_database.sh \
  --proteins IMGVR5_UViG.faa.gz --metadata IMGVR5_UViG.tsv.gz \
  --destination DATABASE_DIR --work-dir BUILD_DIR \
  --threads N --memory-limit 300GB \
  --protein-sha256 HASH --metadata-sha256 HASH \
  --expected-uvigs N --expected-proteins N --expected-representatives N \
  [--cellular-proteins cellular_proteins.faa.gz \
   --cellular-metadata cellular_proteins.tsv \
   --cellular-min-seq-id 0.90 --cellular-coverage 0.80 \
   --cellular-min-protein-length 50]

Cellular metadata must be tab-delimited and contain protein_id,
cellular_group, and source_accession columns. When cellular inputs are
provided, the runtime contains both the established viral-only database and
a separately labeled viral-plus-cellular database for competitive scoring.
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
CELLULAR_PROTEINS=''
CELLULAR_METADATA=''
CELLULAR_MIN_SEQ_ID='0.90'
CELLULAR_COVERAGE='0.80'
CELLULAR_MIN_PROTEIN_LENGTH='50'

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
        --cellular-proteins) CELLULAR_PROTEINS="$2"; shift 2 ;;
        --cellular-metadata) CELLULAR_METADATA="$2"; shift 2 ;;
        --cellular-min-seq-id) CELLULAR_MIN_SEQ_ID="$2"; shift 2 ;;
        --cellular-coverage) CELLULAR_COVERAGE="$2"; shift 2 ;;
        --cellular-min-protein-length) CELLULAR_MIN_PROTEIN_LENGTH="$2"; shift 2 ;;
        *) usage ;;
    esac
done

if [[ -n "$CELLULAR_PROTEINS" || -n "$CELLULAR_METADATA" ]]; then
    [[ -n "$CELLULAR_PROTEINS" && -n "$CELLULAR_METADATA" ]] || {
        echo "ERROR: --cellular-proteins and --cellular-metadata must be supplied together." >&2
        exit 2
    }
fi

for value in PROTEINS METADATA DESTINATION WORK_DIR THREADS MEMORY_LIMIT \
    PROTEIN_SHA256 METADATA_SHA256 EXPECTED_UVIGS EXPECTED_PROTEINS \
    EXPECTED_REPRESENTATIVES; do
    [[ -n "${!value}" ]] || usage
done

for command_name in python mmseqs diamond pigz sha256sum gzip seqkit; do
    command -v "$command_name" >/dev/null 2>&1 || {
        echo "ERROR: Required viCAT build command is unavailable: $command_name" >&2
        exit 1
    }
done

[[ -f "$PROTEINS" ]] || { echo "ERROR: Protein source not found: $PROTEINS" >&2; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: Metadata source not found: $METADATA" >&2; exit 1; }
if [[ -n "$CELLULAR_PROTEINS" ]]; then
    [[ -f "$CELLULAR_PROTEINS" ]] || {
        echo "ERROR: Cellular protein source not found: $CELLULAR_PROTEINS" >&2
        exit 1
    }
    [[ -f "$CELLULAR_METADATA" ]] || {
        echo "ERROR: Cellular metadata source not found: $CELLULAR_METADATA" >&2
        exit 1
    }
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Every invocation verifies the pinned sources, including checkpoint resumes.
printf '%s  %s\n' "$PROTEIN_SHA256" "$PROTEINS" | sha256sum --check -
printf '%s  %s\n' "$METADATA_SHA256" "$METADATA" | sha256sum --check -
if [[ -d "$WORK_DIR/checkpoints" && ! -s "$WORK_DIR/build_config.tsv" ]]; then
    echo 'ERROR: Legacy viCAT build work lacks a provenance stamp; use a new work directory.' >&2
    exit 1
fi
mkdir -p "$WORK_DIR" "$WORK_DIR/checkpoints" "$WORK_DIR/mmseqs" \
    "$WORK_DIR/dedup" "$WORK_DIR/taxonomy" "$WORK_DIR/taxonomy_work" \
    "$WORK_DIR/runtime"

{
    printf 'parameter\tvalue\n'
    printf 'protein_sha256\t%s\nmetadata_sha256\t%s\n' "$PROTEIN_SHA256" "$METADATA_SHA256"
    printf 'uvigs\t%s\nproteins\t%s\nrepresentatives\t%s\n' "$EXPECTED_UVIGS" "$EXPECTED_PROTEINS" "$EXPECTED_REPRESENTATIVES"
    printf 'cellular_identity\t%s\ncellular_coverage\t%s\ncellular_min_length\t%s\n' "$CELLULAR_MIN_SEQ_ID" "$CELLULAR_COVERAGE" "$CELLULAR_MIN_PROTEIN_LENGTH"
    for helper in build_vicat_database.sh prepare_vicat_metadata.py build_vicat_taxonomy_lookup.py prepare_vicat_competitive_references.py; do
        printf '%s\t%s\n' "$helper" "$(sha256sum "$SCRIPT_DIR/$helper" | awk '{print $1}')"
    done
    if [[ -n "$CELLULAR_PROTEINS" ]]; then
        printf 'cellular_proteins\t%s\ncellular_metadata\t%s\n' \
            "$(sha256sum "$CELLULAR_PROTEINS" | awk '{print $1}')" \
            "$(sha256sum "$CELLULAR_METADATA" | awk '{print $1}')"
    fi
} > "$WORK_DIR/build_config.current.tsv"
if [[ -s "$WORK_DIR/build_config.tsv" ]]; then
    cmp -s "$WORK_DIR/build_config.current.tsv" "$WORK_DIR/build_config.tsv" || {
        echo 'ERROR: viCAT build inputs, parameters or code changed; use a new work directory.' >&2
        exit 1
    }
else
    mv "$WORK_DIR/build_config.current.tsv" "$WORK_DIR/build_config.tsv"
fi

checkpoint_done() { [[ -f "$WORK_DIR/checkpoints/$1.complete" ]]; }
mark_complete() { date -Iseconds > "$WORK_DIR/checkpoints/$1.complete"; }

if ! checkpoint_done source_validation; then
    gzip -t "$PROTEINS"
    gzip -t "$METADATA"
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

COMPETITIVE_DB=''
COMPETITIVE_MANIFEST=''
COMPETITIVE_SUMMARY=''
CELLULAR_REPRESENTATIVE_COUNT='0'
if [[ -n "$CELLULAR_PROTEINS" ]]; then
    CELLULAR_FILTERED="$WORK_DIR/dedup/cellular.filtered.faa"
    CELLULAR_DB="$WORK_DIR/mmseqs/cellular_filtered"
    CELLULAR_CLUSTER_DB="$WORK_DIR/mmseqs/cellular_clusters"
    CELLULAR_REP_DB="$WORK_DIR/mmseqs/cellular_representatives"
    CELLULAR_REP_FASTA="$WORK_DIR/dedup/cellular_representatives.faa.gz"

    if ! checkpoint_done cellular_filter; then
        rm -f "$CELLULAR_FILTERED"
        seqkit seq --remove-gaps --min-len "$CELLULAR_MIN_PROTEIN_LENGTH" \
            "$CELLULAR_PROTEINS" > "$CELLULAR_FILTERED"
        [[ "$(grep -c '^>' "$CELLULAR_FILTERED")" -gt 0 ]] || {
            echo "ERROR: No cellular proteins remain after length filtering." >&2
            exit 1
        }
        mark_complete cellular_filter
    fi

    if ! checkpoint_done cellular_cluster; then
        rm -f "${CELLULAR_DB}"* "${CELLULAR_CLUSTER_DB}"* "${CELLULAR_REP_DB}"*
        mmseqs createdb "$CELLULAR_FILTERED" "$CELLULAR_DB" \
            --dbtype 1 --shuffle 0 --createdb-mode 0
        mmseqs linclust "$CELLULAR_DB" "$CELLULAR_CLUSTER_DB" \
            "$WORK_DIR/mmseqs/cellular_tmp" \
            --min-seq-id "$CELLULAR_MIN_SEQ_ID" \
            --cov-mode 0 -c "$CELLULAR_COVERAGE" --threads "$THREADS"
        mmseqs createsubdb "$CELLULAR_CLUSTER_DB" "$CELLULAR_DB" \
            "$CELLULAR_REP_DB" --subdb-mode 1
        mmseqs convert2fasta "$CELLULAR_REP_DB" "${CELLULAR_REP_FASTA%.gz}"
        pigz -p "$THREADS" "${CELLULAR_REP_FASTA%.gz}"
        gzip -t "$CELLULAR_REP_FASTA"
        mark_complete cellular_cluster
    fi

    COMBINED_FASTA="$WORK_DIR/dedup/vicat_viral_cellular.faa"
    COMPETITIVE_MANIFEST="$WORK_DIR/runtime/vicat_competitive_reference_manifest.parquet"
    COMPETITIVE_SUMMARY="$WORK_DIR/runtime/vicat_competitive_database_summary.tsv"
    if ! checkpoint_done competitive_references; then
        rm -f "$COMBINED_FASTA" "$COMPETITIVE_MANIFEST" "$COMPETITIVE_SUMMARY"
        python "$SCRIPT_DIR/prepare_vicat_competitive_references.py" \
            --viral-representatives "$REP_FASTA" \
            --viral-taxonomy-lookup "$TAXONOMY_LOOKUP" \
            --cellular-representatives "$CELLULAR_REP_FASTA" \
            --cellular-metadata "$CELLULAR_METADATA" \
            --output-fasta "$COMBINED_FASTA" \
            --output-manifest "$COMPETITIVE_MANIFEST" \
            --output-summary "$COMPETITIVE_SUMMARY" \
            --work-dir "$WORK_DIR/competitive_manifest_work"
        mark_complete competitive_references
    fi

    COMPETITIVE_DB="$WORK_DIR/runtime/vicat_viral_cellular.dmnd"
    if ! checkpoint_done competitive_diamond_database; then
        rm -f "$COMPETITIVE_DB"
        diamond makedb --in "$COMBINED_FASTA" \
            --db "${COMPETITIVE_DB%.dmnd}" --threads "$THREADS"
        diamond dbinfo --db "$COMPETITIVE_DB" >/dev/null
        mark_complete competitive_diamond_database
    fi

    CELLULAR_REPRESENTATIVE_COUNT="$(awk -F '\t' \
        '$1 == "cellular_representatives" {print $2}' "$COMPETITIVE_SUMMARY")"
    combined_expected=$((EXPECTED_REPRESENTATIVES + CELLULAR_REPRESENTATIVE_COUNT))
    combined_diamond_count="$(diamond dbinfo --db "$COMPETITIVE_DB" \
        | awk '$1 == "Sequences" {print $2}')"
    manifest_count="$(python - "$COMPETITIVE_MANIFEST" <<'PY'
import duckdb
import sys
print(duckdb.connect().execute(
    "SELECT count(*) FROM read_parquet(?)", [sys.argv[1]]
).fetchone()[0])
PY
)"
    [[ "$combined_diamond_count" == "$combined_expected" ]] || {
        echo "ERROR: Competitive DIAMOND count is $combined_diamond_count; expected $combined_expected" >&2
        exit 1
    }
    [[ "$manifest_count" == "$combined_expected" ]] || {
        echo "ERROR: Competitive manifest count is $manifest_count; expected $combined_expected" >&2
        exit 1
    }
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
    printf 'competitive_database_built\t%s\n' "$([[ -n "$COMPETITIVE_DB" ]] && echo true || echo false)"
    printf 'cellular_representative_count\t%s\n' "$CELLULAR_REPRESENTATIVE_COUNT"
    printf 'cellular_clustering_min_seq_id\t%s\n' "$CELLULAR_MIN_SEQ_ID"
    printf 'cellular_clustering_bidirectional_coverage\t%s\n' "$CELLULAR_COVERAGE"
    printf 'cellular_minimum_protein_length\t%s\n' "$CELLULAR_MIN_PROTEIN_LENGTH"
    printf 'mmseqs_version\t%s\n' "$(mmseqs version)"
    printf 'diamond_version\t%s\n' "$(diamond version | awk '{print $3}')"
    printf 'build_completed\t%s\n' "$(date -Iseconds)"
} > "$WORK_DIR/runtime/vicat_database_metadata.tsv"

(
    cd "$WORK_DIR/runtime"
    cp "$WORK_DIR/build_config.tsv" vicat_build_config.tsv
    checksum_files=(vicat_build_config.tsv IMGVR5_UViG_representatives.dmnd \
        IMGVR5_UViG.vicat_taxonomy_lookup.parquet \
        vicat_taxonomy_build_summary.tsv \
        vicat_database_metadata.tsv)
    if [[ -s vicat_viral_cellular.dmnd ]]; then
        checksum_files+=(vicat_viral_cellular.dmnd \
            vicat_competitive_reference_manifest.parquet \
            vicat_competitive_database_summary.tsv)
    fi
    sha256sum "${checksum_files[@]}" > SHA256SUMS
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
