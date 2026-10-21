#!/usr/bin/env python3
"""Isolated main-entry-point check. No production/student paths or search executables."""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import traceback

from entrypoint_reconciliation_fixture import (BUNDLE, build,
    compare_outputs, expected_evidence, sha, verify_command, verify_original, write)

ROOT = Path(__file__).resolve().parents[1]


def verify_package():
    manifest = json.loads((ROOT / 'ENTRYPOINT_PACKAGE.json').read_text())
    for name, digest in manifest['files'].items():
        path = (ROOT / name).resolve()
        if ROOT not in path.parents or sha(path) != digest:
            raise ValueError('Entry-point test package changed: ' + name)
    return manifest


def verify_previous(args, package):
    original = args.original_run.resolve() / 'viSUM-integrated-reconciliation'
    verify_original(original)
    native_status = args.passed_native_run.resolve() / 'results/validation_status.json'
    if sha(native_status) != package['passed_native_status_sha256']:
        raise ValueError('Not the reviewed passed native result')
    native = json.loads(native_status.read_text())
    if native['status'] != 'PASS' or not native['previous_sources_and_outputs_unchanged']:
        raise ValueError('Earlier native validation did not pass')
    if [s['stage'] for s in native['native_nextflow']['stages']] != ['first', 'unchanged_resume', 'manifest_changed', 'bundle_changed']:
        raise ValueError('Earlier native stages incomplete')
    if any(not s['scientific_outputs_identical'] for s in native['native_nextflow']['stages']):
        raise ValueError('Earlier native scientific comparison failed')
    return original


def reference_output(args, sample):
    folder = 'integrated_current_defaults' if sample == 'BMARK_RNA_VIRAL' else 'integrated'
    return args.original_run.resolve() / 'results/full_replay' / sample / folder


def prepare_runtime(output):
    """Offline runtime copy and isolated Nextflow home/cache; no install/update."""
    installed_home = Path(os.environ.get('NXF_HOME', str(Path.home() / '.nextflow')))
    framework = installed_home / 'framework/26.04.6'
    if not framework.is_dir():
        raise ValueError('Offline Nextflow 26.04.6 framework missing: ' + str(framework))
    shutil.copytree(framework, output / 'nextflow-home/framework/26.04.6')
    env = dict(os.environ)
    for name in ('NXF_CLOUDCACHE_PATH', 'NXF_DISABLE_CACHE', 'NXF_CACHE_DIR', 'NXF_PARAMS_FILE',
                 'NXF_WORK', 'NXF_OPTS', 'NXF_ASSETS', 'NXF_DEFAULT_DSL', 'NXF_EXECUTOR'):
        env.pop(name, None)
    env.update(NXF_HOME=str(output / 'nextflow-home'), NXF_PLUGINS_DIR=str(output / 'plugins'),
        NXF_VER='26.04.6', NXF_OFFLINE='true', NXF_DISABLE_CHECK_LATEST='true', NXF_SYNTAX_PARSER='v2',
        PYTHONDONTWRITEBYTECODE='1', TMPDIR='/tmp', TEMP='/tmp', TMP='/tmp')
    validation_bin = output / 'runtime-bin'
    validation_bin.mkdir()
    python = validation_bin / 'python3'
    write(python, '#!/bin/sh\nexec ' + shlex.quote(sys.executable) + ' "$@"\n')
    python.chmod(0o755)
    return env, validation_bin


def command(nextflow, output, project, stage):
    argv = [str(nextflow), '-log', str(output / (stage + '.nextflow.log')),
        'run', str(project / 'visum_nextflow.nf'), '-c', str(output / 'site.config'),
        '-params-file', str(output / ('params-missing.json' if stage == 'missing_manifest' else 'params.json')),
        '-work-dir', str(output / 'work'), '-ansi-log', 'false',
        '-with-trace', str(output / (stage + '.trace.tsv'))]
    if stage == 'unchanged_resume':
        argv.append('-resume')
    return argv


def native(args, output, project, selected):
    if platform.system() != 'Linux':
        raise ValueError('Native entry-point gate requires Linux; local fixtures do not establish it')
    env, runtime_bin = prepare_runtime(output)
    write(output / 'site.config', """// Test-only override; automatic nextflow.config and visum.config still load.
conda.enabled = false
docker.enabled = false
singularity.enabled = false
dag.enabled = false
report.enabled = false
timeline.enabled = false
trace.fields = 'task_id,hash,name,status,exit,workdir'
process.executor = 'local'
process.time = null
process.errorStrategy = 'terminate'
process.maxRetries = 0
process.beforeScript = 'echo "STOP: Only real VIHARMONY may execute in this fixture" >&2; exit 1'
process { withName: 'VIHARMONY' { beforeScript = '' } }
env.PATH = '""" + str(runtime_bin) + ":' + System.getenv('PATH')\n")
    stages, first_hashes = [], {}
    for stage in ('first', 'unchanged_resume', 'missing_manifest'):
        if stage == 'missing_manifest':
            invalid = output / 'missing-reference' / Path(BUNDLE).name
            invalid.parent.mkdir()
            shutil.copyfile(project / BUNDLE, invalid)  # Intentionally no manifest.
            payload = json.loads((output / 'params.json').read_text())
            payload['ictv_reconciliation'] = str(invalid)
            write(output / 'params-missing.json', json.dumps(payload, indent=2) + '\n')
        argv = command(args.nextflow.resolve(), output, project, stage)
        write(output / (stage + '.command.json'), json.dumps(argv, indent=2) + '\n')
        log = output / (stage + '.console.log')
        with log.open('w', encoding='utf-8') as handle:
            process = subprocess.run(argv, cwd=output, env=env, stdout=handle, stderr=subprocess.STDOUT)
        trace = output / (stage + '.trace.tsv')
        rows = []
        if trace.is_file():
            with trace.open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
        if stage == 'missing_manifest':
            if process.returncode == 0 or rows or 'Pinned ICTV reconciliation bundle or manifest was not found' not in log.read_text():
                raise AssertionError('Missing manifest did not fail safely before tasks: ' + str(log))
            stages.append(dict(stage=stage, status='PASS', expected_failure=True, tasks_started=0))
            print('PASS main entry point: missing manifest rejected before any tasks', flush=True)
            continue
        if process.returncode != 0:
            raise AssertionError('Main entry point failed: ' + str(log))
        if len(rows) != 2 or {r['name'].split(' (')[0].split(':')[-1] for r in rows} != {'VIHARMONY'}:
            raise AssertionError('Expected ONLY two actual VIHARMONY tasks: ' + str(rows))
        comparisons, commands = [], []
        for row in rows:
            if row['exit'] != '0' or row['status'] != ('COMPLETED' if stage == 'first' else 'CACHED'):
                raise AssertionError('Unexpected task/cache outcome: ' + str(row))
            sample = next((s for s in selected if '(' + s['sample'] + ')' in row['name']), None)
            if sample is None:
                raise AssertionError('Unexpected task sample')
            sid = sample['sample']
            work = Path(row['workdir'])
            if stage == 'first':
                first_hashes[sid] = row['hash']
            elif first_hashes[sid] != row['hash']:
                raise AssertionError('Unchanged resume task hash changed')
            commands.append(verify_command(work / '.command.sh', sample, project / BUNDLE))
            task_copy = output / 'task-commands' / stage / (sid + '.command.sh')
            task_copy.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(work / '.command.sh', task_copy)
            for source in [sample[key] for key in ('normalized', 'header', 'gate', 'refined', 'regions', 'groups')] + expected_evidence(sample):
                if sha(work / Path(source).name) != sha(Path(source)):
                    raise AssertionError('Wrong/missing staged saved-input bytes: ' + source)
            for name in (BUNDLE, BUNDLE + '.manifest.json', 'assets/ICTV_VMR_MSL41.csv'):
                if sha(work / Path(name).name) != sha(project / name):
                    raise AssertionError('Wrong/missing staged reference bytes: ' + name)
            actual = output / 'published' / (sid + '_results/viharmony')
            comparisons.append(compare_outputs(reference_output(args, sid), actual, sid))
        stages.append(dict(stage=stage, status='PASS', tasks=rows, command_checks=commands,
                           scientific_comparisons=comparisons))
        print('PASS main entry point: ' + stage + '; exact science + reference/evidence staging verified', flush=True)
    return stages


def local_replay(args, output, project, selected):
    """Real saved-input harmonizer checks on this host; explicitly NOT Nextflow."""
    results = []
    for spec in selected:
        for stage, evidence in [('saved_order', expected_evidence(spec)), ('reversed_order', list(reversed(expected_evidence(spec))))]:
            folder = output / 'local-replay' / spec['sample'] / stage
            folder.mkdir(parents=True)
            argv = [sys.executable, '-B', str(project / 'bin/run_viharmony.py'),
                '--sample-id', spec['sample'], '--input-type', spec['type'],
                '--audit-mode', 'full', '--workers', '1',
                '--vcontact3-min-taxonomy-length', str(spec['minimum_length']),
                '--ictv-msl', str(project / 'assets/ICTV_VMR_MSL41.csv'),
                '--ictv-reconciliation', str(project / BUNDLE)]
            for flag, key in [('normalized-fasta','normalized'), ('header-map','header'),
                              ('discovery-gate','gate'), ('refined-fasta','refined'), ('region-map','regions')]:
                argv.extend(['--'+flag, spec[key]])
            argv.extend(['--evidence', *evidence, '--vcontact3-groups', spec['groups'],
                         '--output-prefix', str(folder / spec['sample'])])
            if spec['type'] == 'rna' and spec['pair_floor']:
                argv.append('--rna-pair-homology-floor')
            with (folder / 'replay.log').open('w') as handle:
                process = subprocess.run(argv, stdout=handle, stderr=subprocess.STDOUT)
            if process.returncode:
                raise AssertionError('Local saved-input replay failed: ' + str(folder))
            results.append(dict(stage=stage, **compare_outputs(reference_output(args, spec['sample']), folder, spec['sample'])))
            print('PASS local saved-input check (NOT Nextflow): ' + spec['sample'] + '/' + stage, flush=True)
    return results


def run(args):
    output = args.output_dir.resolve()
    sources = (args.original_run.resolve(), args.passed_native_run.resolve(), ROOT.resolve())
    if output.exists() or any(output == p or p in output.parents for p in sources):
        raise ValueError('Use a NEW output directory outside all source/test trees')
    output.mkdir(parents=True)
    status = dict(status='RUNNING', scope='Unchanged main viSUM workflow + real VIHARMONY + inert saved-file adapters',
        platform=platform.platform(), python=sys.version, expensive_searches_launched=False,
        database_setup_executed=False, production_or_student_runs_changed=False,
        native_entrypoint_executed=False, upstream_tools_and_preflight_tested=False, deployed=False)
    try:
        package = verify_package()
        original = verify_previous(args, package)
        project, selected, protected = build(original, output)
        if args.local_replay:
            status['local_saved_input_checks'] = local_replay(args, output, project, selected)
            status['status'] = 'LOCAL_CHECKS_PASS_NATIVE_PENDING'
        else:
            if args.nextflow is None or not args.nextflow.is_file():
                raise ValueError('Existing Nextflow launcher required; no download attempted')
            status['stages'] = native(args, output, project, selected)
            status.update(status='PASS', native_entrypoint_executed=True)
        if any(sha(project / name) != digest for name, digest in protected.items()):
            raise AssertionError('Protected main/config/scientific file changed during validation')
        verify_previous(args, package)
        verify_package()
        status['original_package_and_protected_files_unchanged'] = True
    except Exception:
        status.update(status='FAIL', error=traceback.format_exc())
        raise
    finally:
        write(output / 'validation_status.json', json.dumps(status, indent=2) + '\n')
    print('Final status: ' + status['status'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-run', type=Path, required=True)
    parser.add_argument('--passed-native-run', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--nextflow', type=Path)
    parser.add_argument('--local-replay', action='store_true', help='Local saved-input check ONLY; does not establish main Nextflow execution')
    run(parser.parse_args())
