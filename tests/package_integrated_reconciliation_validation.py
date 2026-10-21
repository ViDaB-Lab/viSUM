#!/usr/bin/env python3
"""Package ONLY verified saved benchmarks and explicit isolated source files."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = (
    'bin/run_viharmony.py', 'bin/evidence_schema.py', 'bin/ictv_reconciliation.py',
    'bin/build_ictv_reconciliation_bundle.py',
    'bin/reconcile_benchmark_truth.py', 'bin/replay_viharmony_taxonomy.py',
    'bin/standardize_vicat.py', 'bin/standardize_vitap.py', 'bin/standardize_vcontact3.py',
    'assets/ICTV_VMR_MSL41.csv',
    'assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz',
    'assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz.manifest.json',
    'assets/ictv_reconciliation/LICENSE.txt', 'assets/ictv_reconciliation/README.md',
    'visum.config', 'visum_nextflow.nf', 'nextflow.config', 'modules/local/viharmony.nf',
    'tests/replay_integrated_ictv_reconciliation.py',
    'tests/validate_ictv_reconciliation_replay.py',
    'tests/audit_ictv_unassigned_ranks.py',
    'tests/run_integrated_reconciliation_validation.py',
    'tests/test_ictv_reconciliation.py', 'tests/test_ictv_reference.py',
    'tests/test_viharmony.py', 'tests/test_viharmony_workflow_contracts.py',
    'tests/test_viharmony_reconciliation.py',
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def package(verification, archive):
    report_path = verification / 'full_replay_final/report.json'
    summary = json.loads((verification / 'verification_summary.json').read_text())
    report = json.loads(report_path.read_text())
    if sha(report_path) != summary['full_replay_report_sha256'] or report['status'] != 'PASS':
        raise ValueError('Full replay verification snapshot is not the validated PASS report')
    if report['candidate_records'] != 3068 or report['expected_changed_taxonomy_records_verified'] != 213:
        raise ValueError('Unexpected verified benchmark population')
    for name, expected in summary['source_sha256'].items():
        if sha(ROOT / name.replace('\\', '/')) != expected:
            raise ValueError(f'Integrated scientific source changed after verification: {name}')
    files = {name: ROOT / name for name in SOURCE_FILES}
    files['native_validation.nf'] = ROOT / 'tests/integrated_reconciliation_native.nf'
    files['baseline/bin/run_viharmony.py'] = verification / 'baseline_source/bin/run_viharmony.py'
    files['baseline/bin/evidence_schema.py'] = verification / 'baseline_source/bin/evidence_schema.py'
    files['expected/full_replay_report.json'] = report_path
    files['expected/full_taxonomy_changes.jsonl'] = verification / 'full_replay_final/taxonomy_changes.jsonl'
    files['expected/verification_summary.json'] = verification / 'verification_summary.json'
    standalone = ROOT / 'outputs/ictv_reconciliation_prototype/unassigned_rank_correction_verified/benchmark_replay'
    files['expected/standalone/taxonomy_changes.jsonl'] = standalone / 'taxonomy_changes.jsonl'
    counts = {'inputs': 0, 'fasta': 0, 'published': 0}
    for name, expected in report['input_sha256'].items():
        path = Path(name)
        if sha(path) != expected:
            raise ValueError(f'Saved input changed: {path}')
        text = path.as_posix()
        if '/outputs/virsorter2_review_only_final/' in text:
            target = 'saved/published/' + text.split('/outputs/virsorter2_review_only_final/', 1)[1]
            category = 'published'
        elif '/outputs/consolidated_review/inputs/rna_pair_floor_20260909_113555/' in text:
            target = 'saved/inputs/' + text.split('/inputs/rna_pair_floor_20260909_113555/', 1)[1]
            category = 'inputs'
        elif '/outputs/virsorter2_hallmark_replay_v1/fasta/rna_pair_floor_20260909_113555/' in text:
            target = 'saved/fasta/' + text.split('/fasta/rna_pair_floor_20260909_113555/', 1)[1]
            category = 'fasta'
        elif path.resolve() in {files[key].resolve() for key in files if key.startswith('assets/')}:
            continue
        else:
            raise ValueError(f'Unrecognized verified input, refusing implicit packaging: {path}')
        if target in files:
            raise ValueError(f'Duplicate package identity: {target}')
        files[target] = path
        counts[category] += 1
    manifest = dict(schema='integrated-reconciliation-linux-v1', candidate_records=3068,
        changed_taxonomy_records=213, unittest_count=86,
        scope='Saved complete-harmonizer replay plus actual VIHARMONY module cache test; NOT a full pipeline search run',
        input_counts=counts, files={name: sha(path) for name, path in sorted(files.items())})
    archive = archive.resolve()
    sidecar = archive.with_name(archive.name + '.sha256')
    if archive.exists() or sidecar.exists():
        raise ValueError('Use a new archive filename')
    archive.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, 'w:gz', compresslevel=6) as tar:
        for name, path in sorted(files.items()):
            info = tar.gettarinfo(str(path), 'viSUM-integrated-reconciliation/' + name)
            info.mode = 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ''
            with path.open('rb') as handle:
                tar.addfile(info, handle)
        import io
        data = (json.dumps(manifest, indent=2) + '\n').encode()
        info = tarfile.TarInfo('viSUM-integrated-reconciliation/PACKAGE_MANIFEST.json')
        info.size, info.mode = len(data), 0o644
        tar.addfile(info, io.BytesIO(data))
    with sidecar.open('w', encoding='ascii', newline='\n') as handle:
        handle.write(sha(archive) + '  ' + archive.name + '\n')
    print(json.dumps(dict(archive=str(archive), sha256=sha(archive), bytes=archive.stat().st_size,
                         files=len(files), saved_input_counts=counts), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verification', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    args = parser.parse_args()
    package(args.verification, args.archive)
