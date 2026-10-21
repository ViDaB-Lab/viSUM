#!/usr/bin/env python3
"""Deterministic private native-validation package, not an alpha release."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

ROOT_NAME = 'viSUM-nonviral-check'
TESTS = ('run_vicat_nonviral_linux_validation.py', 'test_vicat_nonviral_linux_validation.py',
         'test_vicat_nonviral_download.py', 'test_vicat_database_preparation.py',
         'test_setup_workflow_contracts.py', 'test_setup_nextflow_integration.py',
         'package_vicat_nonviral_linux_validation.py')


def package(project, output):
    project, output = Path(project).resolve(), Path(output).resolve()
    if output.exists() or output.with_name(output.name + '.sha256').exists():
        raise FileExistsError('Choose a new archive name; existing packages are not overwritten')
    paths = [project / name for name in ('visum_nextflow.nf', 'visum.config', 'nextflow.config')]
    paths += [p for folder in ('bin', 'modules', 'workflows', 'assets', 'envs')
              for p in (project / folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    paths += [project / 'tests' / name for name in TESTS]
    paths += [project / 'docs/vicat-database-preparation.md', project / 'docs/database-setup.md']
    data, checkout = {}, {}
    for path in sorted(paths):
        if path.is_symlink() or project not in path.resolve().parents:
            raise ValueError('Unsafe source path: ' + str(path))
        name = path.relative_to(project).as_posix()
        original = path.read_bytes()
        data[name] = original.replace(b'\r\n', b'\n') if path.suffix == '.sh' else original
        checkout[name] = hashlib.sha256(original).hexdigest()
    data['README.txt'] = b'''PRIVATE NONVIRAL NATIVE SETUP/REUSE CHECK -- NOT AN ALPHA INSTALLER
Requires Linux, Python 3.10+, Java, Conda, Mamba, and installed Nextflow 26.04.6.
Only the installed Nextflow framework is copied read-only into private storage.
Build environments and all caches/work/output directories are new and private.

1. Portable regressions + actual main setup workflow routing with fake installers.
2. One real listed NCBI bacterial assembly acquisition and blocked-network reuse;
   the real classifier runs, but no real NCBI database is clustered/built.
3. Actual main viSUM --setup/nonviral module/classifier/MMseqs2/DIAMOND/validators
   operate on seven synthetic proteins spanning all six reference classes.
   Only this copied fixture project's viral prepare process, one-entry source
   manifest and downloader host (localhost fixture server) are substituted.
   The production nonviral module, classifier, builder, environment and workflow
   remain byte-for-byte those in the verified package.
4. Repeat setup with explicit original -resume UUID after the fixture server is
   shut down. Both tasks must execute (not CACHED); the real nonviral validation
   must pass and reuse the database without changing its hashes or mtimes.

2 CPU / 8 GB scheduling budget; only a tiny database is built. The build Conda
environment is downloaded into the private test. No student paths, full 480-
assembly acquisition, viral database construction, analyses or benchmarks.
Native success must be reviewed before a new alpha candidate is distributed.
No packages are installed into existing/student environments. No cleanup is run.
Keep the private test directory and upload the printed review tarball on failure
as well as success. The review excludes environments and large source files.
'''
    digest = lambda value: hashlib.sha256(value).hexdigest()
    manifest = dict(schema='visum-nonviral-native-validation-v1',
        files={name: digest(value) for name, value in sorted(data.items())},
        checkout_sha256=checkout, full_installation_tested=False,
        scope='private small Linux setup/acquisition/reuse gate; not deployment or benchmark')
    data['NONVIRAL_PACKAGE.json'] = (json.dumps(manifest, indent=2, sort_keys=True) + '\n').encode()
    data['CHECKSUMS.sha256'] = ''.join(digest(value) + '  ' + name + '\n'
                                     for name, value in sorted(data.items())).encode()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('wb') as stream, gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode='w') as archive:
            for name, value in sorted(data.items()):
                item = tarfile.TarInfo(ROOT_NAME + '/' + name)
                item.size = len(value)
                item.mtime = 0
                item.mode = 0o755 if name.endswith('.sh') else 0o644
                archive.addfile(item, io.BytesIO(value))
    checksum = digest(output.read_bytes())
    output.with_name(output.name + '.sha256').write_bytes((checksum + '  ' + output.name + '\n').encode())
    return dict(archive=str(output), sha256=checksum, files=len(data), bytes=output.stat().st_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True, type=Path)
    parser.add_argument('--archive', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(package(args.project, args.archive), indent=2))
