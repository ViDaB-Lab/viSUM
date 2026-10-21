#!/usr/bin/env python3
"""Continue ONLY the native gate, reading the prior passed Linux replay unchanged."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import traceback

PATCH_ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_patch():
    manifest = json.loads((PATCH_ROOT / 'PATCH_MANIFEST.json').read_text())
    for name, expected in manifest['files'].items():
        path = (PATCH_ROOT / name).resolve()
        if PATCH_ROOT not in path.parents or not path.is_file() or sha(path) != expected:
            raise ValueError('Native continuation package changed: ' + name)
    return manifest


def load_runner(project):
    # Scientific code and replay helpers come from the original checksum-verified
    # package. Only the test runner's command construction is supplied by this fix.
    sys.path.insert(0, str(project / 'tests'))
    source = PATCH_ROOT / 'tests/run_integrated_reconciliation_validation.py'
    spec = importlib.util.spec_from_file_location('ictv_native_corrected_runner', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_previous(previous, expected, runner):
    project = previous / 'viSUM-integrated-reconciliation'
    results = previous / 'results'
    for name, digest in expected['previous_review_files'].items():
        if sha(previous / name) != digest:
            raise ValueError('Prior reviewed result changed or missing: ' + name)
    manifest_path = project / 'PACKAGE_MANIFEST.json'
    if sha(manifest_path) != expected['original_package_manifest_sha256']:
        raise ValueError('Not the original reviewed integrated validation package')
    manifest = json.loads(manifest_path.read_text())
    runner.verify_package(manifest, project)
    status = json.loads((results / 'validation_status.json').read_text())
    if status['regression'] != dict(tests_run=86, failures=0, errors=0, skipped=0):
        raise ValueError('Prior regression checks did not pass')
    replay = json.loads((results / 'full_replay/report.json').read_text())
    required = dict(status='PASS', candidate_records=3068, parallel_workers=8,
        expected_changed_taxonomy_records_verified=213, retention_decisions_identical=True,
        retained_and_provisional_sequences_identical=True, saved_inputs_unchanged=True,
        parallel_outputs_identical=True)
    for key, value in required.items():
        if replay.get(key) != value:
            raise ValueError('Prior full replay did not establish: ' + key)
    for sample, folder in (('BMARK_DNA_CELLULAR', 'integrated'),
                           ('BMARK_RNA_VIRAL', 'integrated_current_defaults')):
        saved = results / 'full_replay' / sample / folder
        output_manifest = json.loads((saved / (sample + '.harmonizer_manifest.json')).read_text())
        for record in output_manifest['outputs'].values():
            name = record['file']
            if Path(name).name != name or not name.startswith(sample + '.') or sha(saved / name) != record['sha256']:
                raise ValueError('Previously generated scientific output changed: ' + name)
    return project, results, manifest, status


def run(args):
    previous = args.previous_run.resolve()
    output = args.output_dir.resolve()
    if output.exists() or previous == output or previous in output.parents:
        raise ValueError('Use a NEW output directory outside the previous run; nothing is overwritten')
    output.mkdir(parents=True)
    status = dict(status='RUNNING', scope='Native VIHARMONY module/cache gate ONLY',
        python=sys.version, platform=platform.platform(), searches_launched=False,
        regression_rerun=False, full_replay_rerun=False, deployed=False)
    try:
        patch = verify_patch()
        status['patch_manifest_sha256'] = sha(PATCH_ROOT / 'PATCH_MANIFEST.json')
        runner = load_runner(previous / 'viSUM-integrated-reconciliation')
        project, old_results, original_manifest, old_status = verify_previous(previous, patch, runner)
        print('PASS: original package, passed 86 tests, passed 3,068-candidate replay, and saved output checksums', flush=True)
        status['reused_passed_regression'] = old_status['regression']
        status['reused_passed_full_replay'] = old_status['full_replay']
        if platform.system() != 'Linux':
            raise ValueError('Execute the native Nextflow gate on Linux; Windows checks do not establish it')
        if not args.nextflow.is_file():
            raise ValueError('Existing Nextflow launcher missing; no download attempted')
        status['native_nextflow'] = runner.native_check(output, args.nextflow.resolve(),
            old_results / 'full_replay', 8, project_root=project)
        # Recheck every original package file and the prior reports after native
        # execution. The original directory is a read-only source to this runner.
        verify_previous(previous, patch, runner)
        verify_patch()
        status.update(status='PASS', previous_sources_and_outputs_unchanged=True)
    except Exception:
        status.update(status='FAIL', error=traceback.format_exc())
        raise
    finally:
        (output / 'validation_status.json').write_text(json.dumps(status, indent=2) + '\n')
    print('PASS: all four native/cache stages; previous full replay reused without rerunning', flush=True)
    print('Report: ' + str(output / 'validation_status.json'), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous-run', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--nextflow', type=Path, required=True)
    run(parser.parse_args())
