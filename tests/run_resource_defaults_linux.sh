#!/usr/bin/env bash
# Frozen-snapshot probes only. No biological workflow, setup, resume or downloads.
set -euo pipefail
[[ "$(uname -s)" == Linux ]] || { printf 'ERROR: Run on Linux.\n' >&2; exit 1; }
SNAPSHOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
[[ -f "$SNAPSHOT_DIR/RESOURCE_TEST_SNAPSHOT.json" ]] || {
    printf 'ERROR: Use the new resource-test snapshot, not installed viSUM.\n' >&2; exit 1;
}

# These existing executables are read/copy sources only; never install into them.
TEST_PYTHON=/home/alex/programs/visumproject/viSUM/.conda/env-6b66a2dc35344c2815b881411fa02025/bin/python
NEXTFLOW_SOURCE=/home/alex/programs/visumproject/nextflow
FRAMEWORK_SOURCE=/home/alex/.nextflow/framework/26.04.6
[[ -x "$TEST_PYTHON" ]] || { printf 'ERROR: Existing Python missing: %s\n' "$TEST_PYTHON" >&2; exit 1; }
[[ -f "$NEXTFLOW_SOURCE" ]] || { printf 'ERROR: Nextflow launcher missing: %s\n' "$NEXTFLOW_SOURCE" >&2; exit 1; }
[[ -s "$FRAMEWORK_SOURCE/nextflow-26.04.6-one.jar" ]] || {
    printf 'ERROR: Nextflow 26.04.6 runtime missing: %s\n' "$FRAMEWORK_SOURCE" >&2; exit 1;
}
command -v java >/dev/null || { printf 'ERROR: Java is not on PATH.\n' >&2; exit 1; }
unset PYTHONHOME PYTHONPATH PYTHONOPTIMIZE NXF_WORK NXF_ASSETS NXF_CACHE_DIR NXF_DISABLE_CACHE \
      NXF_PARAMS_FILE NXF_DEFAULT_DSL NXF_EXECUTOR NXF_CONFIG_FILE NXF_JVM_ARGS \
      JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS
export PYTHONDONTWRITEBYTECODE=1
"$TEST_PYTHON" -B "$SNAPSHOT_DIR/tests/verify_resource_defaults_linux.py" --verify-only

RUN_DIR="$(mktemp -d "$SNAPSHOT_DIR/../resource-defaults-run.XXXXXXXX")"
mkdir -p "$RUN_DIR/bin" "$RUN_DIR/nextflow-home/framework" "$RUN_DIR/tmp" \
         "$RUN_DIR/logs" "$RUN_DIR/artifacts" "$RUN_DIR/xdg-cache"
cp "$NEXTFLOW_SOURCE" "$RUN_DIR/bin/nextflow"
chmod u+x "$RUN_DIR/bin/nextflow"
cp -a "$FRAMEWORK_SOURCE" "$RUN_DIR/nextflow-home/framework/"
cp "$SNAPSHOT_DIR/RESOURCE_TEST_SNAPSHOT.json" "$RUN_DIR/"
export PATH="$RUN_DIR/bin:$PATH"
export NXF_HOME="$RUN_DIR/nextflow-home"
export NXF_PLUGINS_DIR="$RUN_DIR/nextflow-home/plugins"
export NXF_VER=26.04.6 NXF_OFFLINE=true NXF_DISABLE_CHECK_LATEST=true
export NXF_SYNTAX_PARSER=v2
export NXF_OPTS='-Xms128m -Xmx1g -XX:ActiveProcessorCount=2'
export TMPDIR="$RUN_DIR/tmp" XDG_CACHE_HOME="$RUN_DIR/xdg-cache"
cd "$RUN_DIR"
printf '\nIsolated resource-test directory: %s\n' "$RUN_DIR"
printf 'Only config resolution and local printf probes; no analyses or database access.\n'
printf 'Unfinished checks ONLY: Docker, Singularity; defaults, parameters and site probes.\n'

TEST_STATUS=0
if ! nice -n 10 nextflow -version 2>&1 | tee "$RUN_DIR/logs/runtime.log"; then
    TEST_STATUS=1
fi
if (( TEST_STATUS == 0 )); then
    if ! nice -n 10 "$TEST_PYTHON" -B "$SNAPSHOT_DIR/tests/verify_resource_defaults_linux.py" \
        --run-dir "$RUN_DIR" 2>&1 | tee "$RUN_DIR/logs/probes.log"; then
        TEST_STATUS=1
    fi
fi
printf '%s\n' "$TEST_STATUS" > "$RUN_DIR/logs/exit-status.txt"
files=(logs artifacts RESOURCE_TEST_SNAPSHOT.json)
[[ ! -f "$RUN_DIR/report.json" ]] || files+=(report.json)
tar -czf "$RUN_DIR/resource_defaults_review_20261015_v2.tar.gz" -C "$RUN_DIR" "${files[@]}"
printf '\nTest exit status: %s\nSend this review archive: %s\n' \
    "$TEST_STATUS" "$RUN_DIR/resource_defaults_review_20261015_v2.tar.gz"
exit "$TEST_STATUS"
