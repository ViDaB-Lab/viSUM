#!/usr/bin/env python3
"""Self-contained Linux validation; no setup, searches, student paths or deployment.

The optional native check imports the actual VIHARMONY module and a saved-input
copy process ONLY. It verifies module staging/caching, not the full viSUM DAG.
"""
from __future__ import annotations
import argparse
from contextlib import redirect_stdout, redirect_stderr
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import time
import traceback
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from replay_integrated_ictv_reconciliation import run as replay, compare_worker_outputs, content

TESTS = ('test_ictv_reconciliation.py', 'test_ictv_reference.py', 'test_viharmony.py',
         'test_viharmony_workflow_contracts.py', 'test_viharmony_reconciliation.py')
BUNDLE = 'assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_package(manifest, project_root=ROOT):
    for name, expected in manifest['files'].items():
        path = (project_root / name).resolve()
        if project_root not in path.parents or not path.is_file() or sha(path) != expected:
            raise AssertionError(f'Package input/source changed or missing: {name}')


def native_specs(project_root=ROOT):
    specs = []
    for sample in ('BMARK_DNA_CELLULAR', 'BMARK_RNA_VIRAL'):
        source = project_root / 'saved/inputs' / (sample + '_results')
        fasta = project_root / 'saved/fasta' / (sample + '_results')
        published = project_root / 'saved/published' / sample
        manifest = json.loads((published / 'updated' / (sample + '.harmonizer_manifest.json')).read_text())
        files = dict(normalized=fasta / 'prep' / (sample + '.normalized.fasta'),
            refined=fasta / 'refinement' / (sample + '.refined_candidates.fasta'),
            header=source / 'prep' / (sample + '.header_map.tsv'),
            gate=source / 'discovery_gate' / (sample + '.discovery_gate.tsv'),
            regions=source / 'refinement' / (sample + '.provirus_region_map.tsv'),
            groups=source / 'vcontact3' / (sample + '.vcontact3_group_membership.tsv'))
        evidence = []
        for tool in ('genomad', 'virsorter2', 'cenotetaker3', 'deep6', 'deepmicroclass2',
                     'virbot', 'gianthunter', 'vicat', 'checkv', 'tesorter', 'vitap', 'vcontact3'):
            path = (published / (sample + '.virsorter2_evidence.tsv')) if tool == 'virsorter2' else source / tool / (sample + '.' + tool + '_evidence.tsv')
            if path.is_file():
                evidence.append(path)
        evidence.append(source / 'vicat' / (sample + '.vicat_refined_evidence.tsv'))
        if manifest['input_type'] == 'rna':
            evidence.append(source / 'vicat' / (sample + '.vicat_orf_evidence.tsv'))
        paths = [*files.values(), *evidence]
        if len({p.name for p in paths}) != len(paths) or not all(p.is_file() for p in paths):
            raise AssertionError('Native staged inputs missing or basename collision')
        specs.append(dict(sample=sample, type=manifest['input_type'],
            **{key: path.name for key, path in files.items()},
            evidence=[p.name for p in evidence], files=[str(p) for p in paths],
            minimum_length=manifest['policy']['vcontact3_min_taxonomy_length'],
            pair_floor=manifest['policy']['rna_pair_homology_floor']))
    return specs


def compare_reference_only_outputs(first, latest, sample):
    """Reference byte identity changed deliberately; biological outputs may not."""
    manifests = []
    for folder in (first, latest):
        manifest = json.loads((folder / (sample + '.harmonizer_manifest.json')).read_text())
        for record in manifest['outputs'].values():
            if sha(folder / record['file']) != record['sha256']:
                raise AssertionError('Native output checksum mismatch')
        manifests.append(manifest)
    left, right = manifests
    for label, record in left['outputs'].items():
        if content(first / record['file']) != content(latest / right['outputs'][label]['file']):
            raise AssertionError(f'Biological outputs changed after equivalent reference encoding: {sample}/{label}')
        record.pop('sha256')
        right['outputs'][label].pop('sha256')
    for manifest in manifests:
        for key in ('sha256', 'manifest_sha256'):
            manifest['taxonomy_reconciliation'].pop(key)
    if left != right:
        raise AssertionError('Unexpected non-reference manifest change')


def native_command(nextflow, native, config, parameters, stage, project_root=ROOT):
    """-C and -log are GLOBAL options: they must precede the run command."""
    command = [str(nextflow), '-log', str(native / (stage + '.nextflow.log')),
        '-C', str(config), 'run', str(project_root / 'native_validation.nf'),
        '-params-file', str(parameters), '-work-dir', str(native / 'work'),
        '-ansi-log', 'false', '-with-trace', str(native / (stage + '.trace.tsv'))]
    if stage != 'first':
        command.append('-resume')
    return command


def native_check(output, nextflow, replay_output, workers, project_root=ROOT):
    native = output / 'native'
    native.mkdir()
    reference = native / 'reference'
    reference.mkdir()
    bundle = reference / Path(BUNDLE).name
    manifest = bundle.with_name(bundle.name + '.manifest.json')
    shutil.copyfile(project_root / BUNDLE, bundle)
    shutil.copyfile(project_root / (BUNDLE + '.manifest.json'), manifest)
    specs = native_specs(project_root)
    payload = dict(saved_specs=specs, harmonizer_cpus=workers, max_cpus=workers,
        harmonizer_max_rank='species',  # freeze reporting depth for reconciliation-only comparisons
        harmonizer_memory='32 GB', harmonizer_time=None,
        rna_pair_homology_floor=next(s['pair_floor'] for s in specs if s['type'] == 'rna'),
        ictv_csv=str(project_root / 'assets/ICTV_VMR_MSL41.csv'),
        ictv_reconciliation=str(bundle), outdir=str(native / 'results'))
    parameters = native / 'params.json'
    parameters.write_text(json.dumps(payload, indent=2) + '\n')
    # The process uses the EXACT Python selected by this runner, not an unrelated
    # PATH environment or a newly created Conda installation.
    validation_bin = native / 'bin'
    validation_bin.mkdir()
    wrapper = validation_bin / 'python3'
    wrapper.write_text('#!/bin/sh\nexec ' + shlex.quote(sys.executable) + ' "$@"\n', newline='\n')
    wrapper.chmod(0o755)
    config = native / 'native.config'
    config.write_text("process.executor = 'local'\nprocess.time = null\nconda.enabled = false\n"
        "docker.enabled = false\nsingularity.enabled = false\nexecutor.cpus = " + str(workers) + "\n"
        "executor.memory = '64 GB'\nprocess.errorStrategy = 'terminate'\n"
        "env.PATH = '" + str(validation_bin) + ":' + System.getenv('PATH')\n")
    env = dict(os.environ)
    for name in ('NXF_CLOUDCACHE_PATH', 'NXF_DISABLE_CACHE', 'NXF_CACHE_DIR', 'NXF_PARAMS_FILE', 'NXF_WORK', 'NXF_OPTS'):
        env.pop(name, None)
    env.update(NXF_HOME=str(native / 'nextflow-home'), NXF_PLUGINS_DIR=str(native / 'plugins'),
        NXF_OFFLINE='true', NXF_DISABLE_CHECK_LATEST='true', NXF_SYNTAX_PARSER='v2',
        PYTHONDONTWRITEBYTECODE='1', TMPDIR='/tmp', TEMP='/tmp', TMP='/tmp')
    # No downloads: reuse the installed local Nextflow binary/framework.
    env['NXF_VER'] = '26.04.6'
    installed_home = Path(os.environ.get('NXF_HOME', str(Path.home() / '.nextflow')))
    framework = installed_home / 'framework/26.04.6'
    if not framework.is_dir():
        raise ValueError(f'Offline Nextflow 26.04.6 framework missing: {framework}; no download attempted')
    shutil.copytree(framework, native / 'nextflow-home/framework/26.04.6')
    stages = []
    for stage in ('first', 'unchanged_resume', 'manifest_changed', 'bundle_changed'):
        if stage == 'manifest_changed':
            manifest.write_bytes(manifest.read_bytes() + b'\n  \n')
        if stage == 'bundle_changed':
            bundle.write_bytes(gzip.compress(gzip.decompress(bundle.read_bytes()), compresslevel=1, mtime=0))
            value = json.loads(manifest.read_text())
            value['bundle_sha256'] = sha(bundle)
            manifest.write_text(json.dumps(value, indent=2) + '\n')
        trace = native / (stage + '.trace.tsv')
        command = native_command(nextflow, native, config, parameters, stage, project_root)
        with (native / (stage + '.console.log')).open('w') as log:
            process = subprocess.run(command, cwd=native, env=env, stdout=log, stderr=subprocess.STDOUT)
        if process.returncode:
            raise AssertionError(f'Native {stage} failed: {native / (stage + ".console.log")}')
        with trace.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        if len(rows) != 4:
            raise AssertionError(f'Expected exactly two staging and two VIHARMONY tasks: {stage}/{rows}')
        for row in rows:
            name = row['name'].split(' ')[0].split(':')[-1]
            if name not in {'STAGE_SAVED_INPUTS', 'VIHARMONY'}:
                raise AssertionError(f'Unexpected process scheduled: {row}')
            expected = ('COMPLETED' if stage == 'first' else 'CACHED') if name == 'STAGE_SAVED_INPUTS' else ('CACHED' if stage == 'unchanged_resume' else 'COMPLETED')
            if row['status'] != expected:
                raise AssertionError(f'Unexpected native cache outcome: {stage}/{row}; expected {expected}')
        snapshot = native / (stage + '-outputs')
        shutil.copytree(native / 'results', snapshot)
        for spec in specs:
            sample = spec['sample']
            actual = snapshot / (sample + '_results/viharmony')
            reference_output = replay_output / sample / ('integrated_current_defaults' if spec['type'] == 'rna' else 'integrated')
            if stage in ('first', 'unchanged_resume'):
                compare_worker_outputs(reference_output, actual, sample)
            else:
                original = native / 'first-outputs' / (sample + '_results/viharmony')
                compare_reference_only_outputs(original, actual, sample)
        stages.append(dict(stage=stage, tasks=rows, scientific_outputs_identical=True))
        print(f'PASS native Nextflow: {stage}', flush=True)
    return dict(status='PASS', scope='actual VIHARMONY module + saved-input copy process; no discovery/search modules',
                full_visum_entrypoint_executed=False, stages=stages)


def run(args):
    output = args.output_dir.resolve()
    if output.exists() or ROOT == output or ROOT in output.parents:
        raise ValueError('Use a new output directory OUTSIDE the extracted package')
    output.mkdir(parents=True)
    status = dict(status='RUNNING', python=sys.version, platform=platform.platform(),
                  no_search_programs_launched=True, deployed=False)
    try:
        manifest = json.loads((ROOT / 'PACKAGE_MANIFEST.json').read_text())
        verify_package(manifest)
        status['package_manifest_sha256'] = sha(ROOT / 'PACKAGE_MANIFEST.json')
        print('PASS package/source/saved-input checksums', flush=True)
        suite = unittest.TestSuite()
        for pattern in TESTS:
            suite.addTests(unittest.defaultTestLoader.discover(str(ROOT / 'tests'), pattern=pattern))
        with (output / 'regression.log').open('w') as log, redirect_stdout(log), redirect_stderr(log):
            result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
        status['regression'] = dict(tests_run=result.testsRun, failures=len(result.failures),
                                  errors=len(result.errors), skipped=len(result.skipped))
        if not result.wasSuccessful() or result.skipped or result.testsRun != manifest['unittest_count']:
            raise AssertionError(f'Regression count/failure/skip: {status["regression"]}; see regression.log')
        print(f'PASS {result.testsRun} dependency-free regression tests', flush=True)
        replay_output = output / 'full_replay'
        replay(argparse.Namespace(saved_results=ROOT / 'saved/published', input_root=ROOT / 'saved/inputs',
            fasta_root=ROOT / 'saved/fasta', frozen_script=ROOT / 'baseline/bin/run_viharmony.py',
            bundle=ROOT / BUNDLE, expected_taxonomy=ROOT / 'expected/standalone',
            output=replay_output, parallel_workers=args.workers))
        report = json.loads((replay_output / 'report.json').read_text())
        if report['candidate_records'] != manifest['candidate_records'] or not report['parallel_outputs_identical']:
            raise AssertionError('Complete parallel replay did not run every expected candidate')
        expected_changes = (ROOT / 'expected/full_taxonomy_changes.jsonl').read_text().splitlines()
        actual_changes = (replay_output / 'taxonomy_changes.jsonl').read_text().splitlines()
        if [json.loads(s) for s in actual_changes] != [json.loads(s) for s in expected_changes]:
            raise AssertionError('Detailed taxonomy changes differ from the verified Windows replay')
        status['full_replay'] = dict(status=report['status'], candidate_records=report['candidate_records'],
            taxonomy_changed=report['expected_changed_taxonomy_records_verified'], workers=args.workers,
            identical_viral_retention=True, identical_worker_outputs=True,
            exact_verified_taxonomy_changes=True)
        if args.nextflow:
            if platform.system() != 'Linux':
                raise ValueError('Native Nextflow gate requires Linux; do not claim Windows execution proves it')
            status['native_nextflow'] = native_check(output, args.nextflow.resolve(), replay_output, args.workers)
        else:
            status['native_nextflow'] = dict(status='NOT_RUN', reason='No --nextflow supplied')
        verify_package(manifest)
        status['status'] = 'PASS' if args.nextflow else 'REPLAY_PASS_NATIVE_PENDING'
        status['package_and_saved_inputs_unchanged'] = True
    except Exception:
        status.update(status='FAIL', error=traceback.format_exc())
        raise
    finally:
        (output / 'validation_status.json').write_text(json.dumps(status, indent=2) + '\n')
    print(f'Final status: {status["status"]}; {output / "validation_status.json"}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--nextflow', type=Path)
    args = parser.parse_args()
    if args.workers < 2:
        parser.error('--workers must be at least 2 for the parallel-equivalence gate')
    run(args)
