#!/usr/bin/env bash
# Run on Linux. Downloads real databases after the download-free routing test.
set -euo pipefail
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${NEXTFLOW_BIN:-}" ]]; then
    nextflow_bin="$NEXTFLOW_BIN"
elif command -v nextflow >/dev/null 2>&1; then
    nextflow_bin="$(command -v nextflow)"
else
    nextflow_bin="$PROJECT/../nextflow"
fi
[[ -x "$nextflow_bin" ]] || { echo "ERROR: Nextflow executable missing: $nextflow_bin" >&2; exit 1; }
command -v conda >/dev/null
export NEXTFLOW_BIN="$nextflow_bin"
export PATH="$(dirname "$nextflow_bin"):$PATH"
export NXF_OPTS="${NXF_OPTS:-} -XX:ActiveProcessorCount=${SETUP_CPUS:-16}"
cd "$PROJECT"
python -m unittest discover -s tests -p 'test_setup_nextflow_integration.py' -v

# A new directory is mandatory; never reuse the production installation.
TEST_ROOT="$(mktemp -d "${SETUP_TEST_PARENT:-$PROJECT/..}/visum-clean-setup.XXXXXXXX")"
export VISUM_SETUP_TEST_ROOT="$TEST_ROOT"
printf 'Clean setup test directory: %s\n' "$TEST_ROOT"
printf '%s\n' 'viCAT excluded: a separate source-build test requires its input files.'
git rev-parse HEAD > "$TEST_ROOT/commit.txt"
cat > "$TEST_ROOT/isolated.config" <<'CONFIG'
conda.cacheDir = "${System.getenv('VISUM_SETUP_TEST_ROOT')}/conda"
workDir = "${System.getenv('VISUM_SETUP_TEST_ROOT')}/work"
CONFIG
cd "$TEST_ROOT"
for pass in install reuse; do
    bash "$PROJECT/visum" -c "$PROJECT/visum.config" -c "$TEST_ROOT/isolated.config" \
        --setup --dbdir "$TEST_ROOT/databases" --tooldir "$TEST_ROOT/tools" \
        --outdir "$TEST_ROOT/$pass" \
        --run_vicat false --run_vitap true --run_vcontact3 true --run_tesorter false \
        --max_cpus "${SETUP_CPUS:-16}" --max_memory "${SETUP_MEMORY:-128 GB}" \
        -ansi-log false -with-trace "$TEST_ROOT/$pass.trace.tsv" \
        2>&1 | tee "$TEST_ROOT/$pass.log"
    cp .nextflow.log "$TEST_ROOT/$pass.nextflow.log"
done
python - "$TEST_ROOT" <<'PY'
import csv
from pathlib import Path
import sys
root = Path(sys.argv[1])
expected = {
    'PREPARE_GENOMAD_DATABASE', 'PREPARE_VIRSORTER2_DATABASE',
    'PREPARE_CENOTETAKER3_DATABASE', 'PREPARE_DEEP6_DATABASE',
    'PREPARE_DEEPMICROCLASS2', 'PREPARE_VIRBOT_DATABASE',
    'PREPARE_GIANTHUNTER_DATABASE', 'PREPARE_CHECKV_DATABASE',
    'PREPARE_VITAP_DATABASE', 'PREPARE_VCONTACT3_DATABASE',
}
for run in ('install', 'reuse'):
    with (root / f'{run}.trace.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    names = {r['name'].split(' (')[0].split(':')[-1] for r in rows}
    if names != expected or len(rows) != len(expected):
        raise SystemExit(f'FAIL {run}: unexpected process set: {names}')
    if any(r['status'] != 'COMPLETED' for r in rows):
        raise SystemExit(f'FAIL {run}: preparation did not execute successfully')
    print(f'PASS {run}: all ten selected preparations completed; no analyses')
print('PASS setup smoke test. viCAT and analysis-runtime smoke tests remain separate.')
print(f'Logs and provenance: {root}')
PY
