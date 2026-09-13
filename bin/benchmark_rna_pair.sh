#!/usr/bin/env bash
# Reuse the established benchmark work cache, publishing into a NEW directory.
# Run as a child script: bash bin/benchmark_rna_pair.sh MANIFEST INDIR OUTDIR
set -euo pipefail
if [[ $# -ne 3 ]]; then
    printf 'Usage: bash %s MANIFEST INDIR NEW_OUTDIR\n' "$0" >&2
    exit 2
fi
PROJECT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
MANIFEST=$(realpath -- "$1")
INDIR=$(realpath -- "$2")
OUTDIR=$(realpath -m -- "$3")
if [[ -e "$OUTDIR" ]]; then
    printf 'Refusing to overwrite existing output: %s\n' "$OUTDIR" >&2
    exit 2
fi
mkdir -p -- "$OUTDIR"
cd -- "$PROJECT"
bash "$PROJECT/visum" -ansi-log true \
    -c "$PROJECT/visum.config" -profile local_full \
    -work-dir "$PROJECT/benchmark/part1/work/vicat_nonviral_release_candidate_v1" -resume \
    --prefix_many "$MANIFEST" --indir "$INDIR" --outdir "$OUTDIR" \
    --run_genomad true --run_virsorter2 true --run_cenotetaker3 true \
    --run_deep6 true --run_deepmicroclass2 true --run_virbot true \
    --run_gianthunter true --run_vicat true --run_checkv true \
    --run_tesorter true --run_vitap true --run_vcontact3 true \
    --vicat_db "$PROJECT/databases/vicat_build" \
    --vicat_nonviral_db "$PROJECT/databases/vicat/nonviral_classaware_v1" \
    --vitap_db /work/databases/vitapdb/vitapupdate/DB_VMR-MSL41 \
    --vcontact3_db "$PROJECT/databases/vcontact3/releases/232" \
    --allow_ct3_only_refinement false --vicat_provirus_min_overlap 0.5 \
    --harmonizer_audit full --rna_pair_homology_floor true \
    --max_cpus 52 --max_memory '400 GB' \
    --genomad_cpus 4 --virsorter2_cpus 4 --ct3_cpus 4 --deep6_cpus 1 \
    --deepmicroclass2_cpus 4 --virbot_cpus 4 --gianthunter_cpus 4 \
    --vicat_orf_cpus 4 --vicat_cpus 16 --vicat_standardizer_cpus 16 \
    --checkv_cpus 4 --tesorter_cpus 4 --vitap_cpus 8 --vcontact3_cpus 2 \
    --vicat_memory '48 GB' --vicat_standardizer_memory '64 GB' \
    --vicat_time 24h --vicat_standardizer_time 24h
