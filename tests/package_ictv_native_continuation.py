#!/usr/bin/env python3
"""Small test-runner fix; retains the entire original package and passed replay."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[1]
FILES = ('tests/run_integrated_reconciliation_validation.py',
         'tests/resume_integrated_reconciliation_native.py', 'tests/test_ictv_native_command.py')


def package(review, archive):
    with tarfile.open(review) as tar:
        reviewed = {name: tar.extractfile(name).read() for name in (
            'results/validation_status.json', 'results/full_replay/report.json',
            'results/full_replay/taxonomy_changes.jsonl', 'results/native/first.console.log')}
    status = json.loads(reviewed['results/validation_status.json'])
    if status['full_replay']['status'] != 'PASS' or reviewed['results/native/first.console.log'].decode().strip() != "Unknown option: -C -- Check the available commands and options and syntax with 'help'":
        raise ValueError('Not the reviewed passed-science / failed-global-option run')
    data = {name: (ROOT / name).read_bytes() for name in FILES}
    manifest = dict(scope='Native gate only; scientific source unchanged; no previous package files edited',
        original_package_manifest_sha256=status['package_manifest_sha256'],
        previous_review_files={name: hashlib.sha256(value).hexdigest() for name, value in reviewed.items()},
        files={name: hashlib.sha256(value).hexdigest() for name, value in data.items()})
    data['PATCH_MANIFEST.json'] = (json.dumps(manifest, indent=2)+'\n').encode()
    if archive.exists():
        raise ValueError('Use a new archive filename')
    archive.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, 'w:gz') as tar:
        for name, value in data.items():
            info = tarfile.TarInfo('ictv-native-fix/' + name)
            info.size, info.mode = len(value), 0o644
            tar.addfile(info, io.BytesIO(value))
    print(json.dumps(dict(archive=str(archive), bytes=archive.stat().st_size,
                         sha256=hashlib.sha256(archive.read_bytes()).hexdigest()), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    args = parser.parse_args()
    package(args.review, args.archive)
