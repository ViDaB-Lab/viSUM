#!/usr/bin/env python3
"""Verify the frozen local candidate identity. Does not install or execute tools."""
import hashlib
import json
from pathlib import Path, PurePosixPath


def verify(root):
    root = Path(root).resolve()
    manifest = json.loads((root / 'SOURCE_MANIFEST.json').read_text())
    for name, record in manifest['files'].items():
        relative = PurePosixPath(name)
        path = (root / name).resolve()
        if relative.is_absolute() or '..' in relative.parts or root not in path.parents:
            raise ValueError('Unsafe candidate manifest path')
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
            raise ValueError('Candidate file missing/changed: ' + name)
    checks = {}
    for line in (root / 'checksums.sha256').read_text().splitlines():
        digest, name = line.split('  ', 1)
        if name in checks or (root / name).resolve().parent != root and name not in manifest['files']:
            raise ValueError('Unexpected/duplicate checksum path')
        checks[name] = digest
    if set(checks) != set(manifest['files']) | {'SOURCE_MANIFEST.json'}:
        raise ValueError('Checksum inventory does not match candidate manifest')
    for name, digest in checks.items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Checksum mismatch: ' + name)
    return dict(status='PASS', release_id=manifest['release_id'], files=len(manifest['files']),
                manifest_sha256=hashlib.sha256((root / 'SOURCE_MANIFEST.json').read_bytes()).hexdigest(),
                external_database_and_environment_runtime_identities='NOT_VERIFIED_BY_SOURCE_CHECK',
                pipeline_launched=False)


if __name__ == '__main__':
    print(json.dumps(verify(Path(__file__).resolve().parents[1]), indent=2))
