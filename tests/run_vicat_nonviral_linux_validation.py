#!/usr/bin/env python3
"""Native nonviral setup/reuse gate. No student paths, full DBs or analyses.

The real setup workflow/nonviral module/classifier/builder run on tiny synthetic
sources, served locally. Only the copied fixture's download host, manifest and
viral-preparation process are substituted. Real NCBI acquisition is checked
separately for one listed bacterial assembly, without building its database.
"""
import argparse
from contextlib import redirect_stdout
import csv
import functools
import gzip
import hashlib
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import io
import json
import os
from pathlib import Path
import platform
import random
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import threading
import traceback
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import download_vicat_nonviral_sources as acquisition

TESTS = ('test_vicat_nonviral_download', 'test_vicat_database_preparation',
         'test_setup_workflow_contracts', 'test_setup_nextflow_integration',
         'test_vicat_nonviral_linux_validation')
FRAMEWORK_VERSION = '26.04.6'


def sha(path):
    return acquisition.digest(path)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def verify_package(root=ROOT):
    manifest = json.loads((root / 'NONVIRAL_PACKAGE.json').read_text())
    for name, expected in manifest['files'].items():
        path = (root / name).resolve()
        if root.resolve() not in path.parents or not path.is_file() or sha(path) != expected:
            raise ValueError('Package file missing/changed: ' + name)
    return manifest


def private_environment(output, installed_home, nextflow):
    framework = Path(installed_home).resolve() / 'framework' / FRAMEWORK_VERSION
    if not (framework / f'nextflow-{FRAMEWORK_VERSION}-one.jar').is_file():
        raise ValueError(f'Installed Nextflow {FRAMEWORK_VERSION} framework not found: {framework}')
    home = output / 'nextflow-home'
    shutil.copytree(framework, home / 'framework' / FRAMEWORK_VERSION)
    environment = dict(os.environ)
    for name in ('NXF_HOME', 'NXF_PLUGINS_DIR', 'NXF_WORK', 'NXF_CACHE_DIR', 'NXF_CLOUDCACHE_PATH',
                 'NXF_DISABLE_CACHE', 'NXF_PARAMS_FILE', 'NXF_OPTS', 'NXF_VER', 'PYTHONPATH',
                 'PYTHONHOME', 'CONDA_ENVS_PATH'):
        environment.pop(name, None)
    environment.update(
        NXF_HOME=str(home), NXF_PLUGINS_DIR=str(home / 'plugins'), NXF_OFFLINE='true',
        NXF_VER=FRAMEWORK_VERSION, NXF_SYNTAX_PARSER='v2', NXF_DISABLE_CHECK_LATEST='true',
        CONDA_PKGS_DIRS=str(output / 'conda-pkgs'), CONDA_ENVS_PATH=str(output / 'conda-envs'),
        MAMBA_ROOT_PREFIX=str(output / 'mamba-root'), XDG_CACHE_HOME=str(output / 'cache'),
        PYTHONUSERBASE=str(output / 'python-user'), PYTHONDONTWRITEBYTECODE='1',
        TMPDIR='/tmp', TEMP='/tmp', TMP='/tmp',
        PATH=str(Path(nextflow).resolve().parent) + os.pathsep + environment.get('PATH', ''))
    return environment


def make_remote_fixture(remote):
    accession = 'GCF_000000001.1'
    basename = accession + '_viSUM_fixture'
    directory = remote / 'GCF/000/000/001' / basename
    directory.mkdir(parents=True)
    # Enough length/complexity for the real MMseqs2 linclust, not a tool mock.
    generator = random.Random(16019)
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    fasta = ''.join(f'>P{i}\nM' + ''.join(generator.choice(alphabet) for _ in range(199)) + '\n'
                    for i in range(7))
    proteins = directory / (basename + '_protein.faa.gz')
    proteins.write_bytes(gzip.compress(fasta.encode(), mtime=0))
    feature = '# feature\tseq_type\tgenomic_accession\tproduct_accession\n'
    types = ('chromosome', 'unplaced scaffold', 'plasmid', 'plastid', 'mitochondrion')
    feature += ''.join(f'CDS\t{kind}\tNC_{i}\tP{i}\n' for i, kind in enumerate(types))
    feature += 'CDS\tchromosome\tNC_C\tP5\nCDS\tplasmid\tNC_P\tP5\n'
    features = directory / (basename + '_feature_table.txt.gz')
    features.write_bytes(gzip.compress(feature.encode(), mtime=0))
    (directory / 'md5checksums.txt').write_text(''.join(
        acquisition.digest(path, 'md5') + '  ./' + path.name + '\n' for path in (proteins, features)))
    return 'assembly_accession\tcellular_group\torganism_name\ttaxid\n' + accession + '\tbacteria\tSynthetic fixture\t1\n'


def substitute_viral_prepare(original):
    inputs = original.split('    input:', 1)[1].split('    output:', 1)[0]
    outputs = original.split('    output:', 1)[1].split('    script:', 1)[0]
    return ('process PREPARE_VICAT_DATABASE {\n    cache false\n'
            '    input:' + inputs + '\n    output:' + outputs + '''
    script:
    """
    mkdir vicat_database
    printf 'synthetic viral dependency only\\n' > vicat_database/fixture.txt
    printf 'fixture\\n' > vicat_database_setup_metadata.tsv
    """
}
''')


def prepare_fixture_project(project, url, manifest, root=ROOT):
    for folder in ('bin', 'envs', 'modules', 'workflows', 'assets'):
        shutil.copytree(root / folder, project / folder, ignore=shutil.ignore_patterns('__pycache__'))
    for name in ('visum_nextflow.nf', 'visum.config', 'nextflow.config'):
        shutil.copyfile(root / name, project / name)
    source = project / 'bin/download_vicat_nonviral_sources.py'
    text = source.read_text()
    needle = "NCBI_ROOT = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'"
    if text.count(needle) != 1:
        raise ValueError('Downloader host substitution no longer matches exactly')
    source.write_text(text.replace(needle, 'NCBI_ROOT = ' + repr(url)))
    path = project / 'assets/vicat/nonviral_classaware_v1_20260904/assembly_manifest.tsv'
    path.write_text(manifest)
    viral = project / 'modules/local/vicat_database.nf'
    viral.write_text(substitute_viral_prepare(viral.read_text()))
    # Prove every other runtime byte was copied as packaged, including biology.
    changed = {'bin/download_vicat_nonviral_sources.py',
               'modules/local/vicat_database.nf',
               'assets/vicat/nonviral_classaware_v1_20260904/assembly_manifest.tsv'}
    for path in project.rglob('*'):
        if path.is_file() and path.relative_to(project).as_posix() not in changed:
            assert path.read_bytes() == (root / path.relative_to(project)).read_bytes(), str(path)
    return sorted(changed)


def quote(path):
    return "'" + str(path).replace('\\', '/').replace("'", "\\'") + "'"


def fixture_config(project, output):
    tools = ('genomad', 'virsorter2', 'cenotetaker3', 'deep6', 'deepmicroclass2',
             'virbot', 'gianthunter', 'checkv', 'vitap', 'vcontact3', 'tesorter')
    config = 'includeConfig ' + quote(project / 'visum.config') + '\n'
    config += '\n'.join('params.run_' + tool + ' = false' for tool in tools) + '\n'
    config += '''params.run_vicat = true
params.vicat_nonviral_auto_download = true
params.max_cpus = 2
params.max_memory = '8 GB'
params.vicat_nonviral_build_cpus = 2
params.vicat_nonviral_build_memory = '2 GB'
process.executor = 'local'
process.errorStrategy = 'terminate'
process.maxRetries = 0
process.time = null
docker.enabled = false
singularity.enabled = false
conda.enabled = true
conda.useMamba = true
report.enabled = false
timeline.enabled = false
dag.enabled = false
trace.fields = 'task_id,hash,name,status,exit,workdir'
'''
    config += 'params.vicat_db = ' + quote(output / 'unused-viral-fixture') + '\n'
    config += 'params.dbdir = ' + quote(output / 'databases') + '\n'
    config += 'params.tooldir = ' + quote(output / 'tools') + '\n'
    config += 'params.outdir = ' + quote(output / 'published') + '\n'
    config += 'conda.cacheDir = ' + quote(output / 'conda') + '\n'
    # Fail closed before any unrelated executable can start, even if routing regresses.
    allowed = {'PREPARE_VICAT_DATABASE', 'PREPARE_VICAT_NONVIRAL_DATABASE'}
    processes = set()
    for path in (project / 'modules').rglob('*.nf'):
        processes.update(re.findall(r'process\s+(\w+)\s*\{', path.read_text()))
    guarded = '|'.join(sorted(processes - allowed))
    config += "process { withName: '" + guarded + "' {\n"
    config += "  conda = null\n  beforeScript = 'echo STOP: unrelated process forbidden in nonviral validation >&2; exit 1'\n} }\n"
    return config


def native_command(nextflow, project, native, config, stage, session=None):
    command = [str(nextflow), '-log', str(native / (stage + '.nextflow.log')), '-C', str(config),
               'run', str(project / 'visum_nextflow.nf'), '--setup',
               '-work-dir', str(native / 'work'), '-ansi-log', 'false',
               '-with-trace', str(native / (stage + '.trace.tsv'))]
    if session:
        if not re.fullmatch(r'[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}', session):
            raise ValueError('Invalid explicit resume session')
        command.extend(['-resume', session])
    return command


def snapshot_database(database):
    # Include mtimes: a byte-identical rebuild must not masquerade as reuse.
    result = {}
    for path in database.rglob('*'):
        if path.is_symlink():
            record = dict(kind='symlink', target=os.readlink(path), mtime_ns=path.lstat().st_mtime_ns)
        elif path.is_file():
            record = dict(kind='file', sha256=sha(path), mtime_ns=path.stat().st_mtime_ns)
        else:
            record = dict(kind='directory', mtime_ns=path.stat().st_mtime_ns)
        result[path.relative_to(database).as_posix()] = record
    return result


def check_trace(path, work):
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    expected = {'PREPARE_VICAT_DATABASE', 'PREPARE_VICAT_NONVIRAL_DATABASE'}
    observed = [row['name'].split(' (')[0].split(':')[-1] for row in rows]
    if len(rows) != 2 or set(observed) != expected:
        raise AssertionError('Unexpected/missing task (no analyses allowed): ' + repr(rows))
    for row in rows:
        if row['status'] != 'COMPLETED' or row['exit'] != '0':
            raise AssertionError('Setup must revalidate, not restore a cached task: ' + repr(row))
        if work.resolve() not in Path(row['workdir']).resolve().parents:
            raise AssertionError('Task escaped private work directory')
    return rows


def check_installation_action(path, expected):
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != 1 or rows[0]['installation_action'] != expected or rows[0]['validation'] != 'passed':
        raise AssertionError('Unexpected setup provenance: ' + repr(rows))
    return rows[0]


def native_check(output, nextflow, environment):
    native = output / 'native'
    native.mkdir()
    manifest = make_remote_fixture(native / 'remote')
    requests = []

    class Handler(SimpleHTTPRequestHandler):
        def do_GET(self):
            requests.append(self.path)
            super().do_GET()

        def log_message(self, *_):
            pass

    server = ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(Handler, directory=str(native / 'remote')))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    url = f'http://127.0.0.1:{server.server_port}/'
    project = native / 'fixture-project'
    project.mkdir()
    changes = prepare_fixture_project(project, url, manifest)
    config = native / 'native.config'
    config.write_text(fixture_config(project, native))
    database = native / 'databases/vicat/nonviral_classaware_v1'
    metadata = native / 'published/database_setup/vicat_nonviral_database_setup_metadata.tsv'
    session, first, stages = None, None, []
    try:
        for stage in ('first_setup', 'resume_reuse'):
            if stage == 'resume_reuse':
                # A repeated acquisition now fails: the source server is offline.
                server.shutdown()
                server.server_close()
                thread.join(timeout=5)
            command = native_command(nextflow, project, native, config, stage, session)
            (native / (stage + '.command.txt')).write_text(shlex.join(command) + '\n')
            print('START native ' + stage + ' (logs: ' + str(native / (stage + '.console.log')) + ')', flush=True)
            with (native / (stage + '.console.log')).open('w') as log:
                result = subprocess.run(command, cwd=native, env=environment, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise AssertionError('Nextflow failed; see ' + str(native / (stage + '.console.log')))
            tasks = check_trace(native / (stage + '.trace.tsv'), native / 'work')
            if session is None:
                matches = set(re.findall(r'Session UUID:\s*([0-9a-f-]{36})', (native / (stage + '.nextflow.log')).read_text()))
                if len(matches) != 1:
                    raise AssertionError('Cannot identify original cache session')
                session = matches.pop()
                if not (native / '.nextflow/cache' / session / 'db').is_dir():
                    raise AssertionError('Original cache is outside private launch directory')
            expected = 'downloaded-and-built-from-assembly-manifest' if stage == 'first_setup' else 'reused-existing-database'
            record = check_installation_action(metadata, expected)
            nonviral_task = next(row for row in tasks if row['name'].split(' (')[0].endswith('PREPARE_VICAT_NONVIRAL_DATABASE'))
            if metadata.read_bytes() != (Path(nonviral_task['workdir']) / metadata.name).read_bytes():
                raise AssertionError('Published setup metadata differs from current task')
            copied = native / (stage + '.setup_metadata.tsv')
            shutil.copyfile(metadata, copied)
            current = snapshot_database(database)
            if first is None:
                first = current
                if len(requests) != 4:
                    raise AssertionError('Expected one index, checksum list and two source downloads: ' + repr(requests))
                summary = database / 'vicat_nonviral_classification_summary.tsv'
                with summary.open(newline='') as handle:
                    values = list(csv.DictReader(handle, delimiter='\t'))
                count = next(int(row['value']) for row in values
                             if row['metric'] == 'included_protein_count' and row['category'] == 'all')
                if count != 7:
                    raise AssertionError('Tiny fixture did not classify all seven proteins')
            elif current != first or len(requests) != 4:
                raise AssertionError('Reuse modified database files or contacted source server')
            stages.append(dict(stage=stage, tasks=tasks, metadata=record))
            print('PASS native ' + stage, flush=True)
    finally:
        if thread.is_alive():
            server.shutdown()
            server.server_close()
            thread.join(timeout=5)
    write_json(native / 'fixture_http_requests.json', requests)
    write_json(native / 'database_snapshot.json', first)
    provenance = json.loads((database / 'vicat_nonviral_source_provenance.json').read_text())
    assert provenance['download_sources']['assembly_count'] == 1
    # Copy only small build/reference identity records into the review bundle.
    receipts = output / 'receipts'
    receipts.mkdir()
    for name in ('SHA256SUMS', 'vicat_nonviral_build_config.tsv', 'vicat_nonviral_database_metadata.tsv',
                 'vicat_nonviral_cluster_summary.tsv', 'vicat_nonviral_classification_summary.tsv',
                 'vicat_nonviral_source_manifest.tsv', 'vicat_nonviral_source_provenance.json'):
        shutil.copyfile(database / name, receipts / name)
    prefixes = list((native / 'conda').glob('env-*'))
    if len(prefixes) != 1:
        raise AssertionError('Expected exactly one private build environment: ' + repr(prefixes))
    result = subprocess.run(['conda', 'list', '--prefix', str(prefixes[0]), '--explicit'], env=environment,
                            capture_output=True, text=True, check=True)
    (output / 'environment.explicit.txt').write_text(result.stdout)
    return dict(status='PASS', original_session=session, stages=stages, fixture_only_substitutions=changes,
                setup_revalidated_instead_of_cached=True, database_bytes_and_mtimes_unchanged_on_reuse=True,
                private_environment=str(prefixes[0]), http_requests=len(requests))


def live_check(output, environment=None):
    inventory = ROOT / 'assets/vicat/nonviral_classaware_v1_20260904/assembly_manifest.tsv'
    row = next(row for row in acquisition.read_inventory(inventory) if row['assembly_accession'] == 'GCF_000006765.1')
    live = output / 'live'
    live.mkdir()
    manifest = live / 'one_assembly.tsv'
    with manifest.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerow(row)
    with (live / 'download.console.log').open('w') as log, redirect_stdout(log):
        proteins, features, provenance = acquisition.download_sources(manifest, live / 'sources', workers=1)
    before = provenance.read_bytes()
    with patch.object(acquisition, 'fetch', side_effect=AssertionError('Network attempted on source reuse')):
        with (live / 'reuse.console.log').open('w') as log, redirect_stdout(log):
            acquisition.download_sources(manifest, live / 'sources', workers=1)
    if provenance.read_bytes() != before:
        raise AssertionError('Live provenance changed on offline reuse')
    # Same real classifier; no clustering/DIAMOND build of this reference.
    with (live / 'classification.console.log').open('w') as log:
        subprocess.run([sys.executable, str(ROOT / 'bin/classify_vicat_nonviral_references.py'),
            '--manifest', str(manifest), '--package-root', str(proteins), '--metadata-root', str(features),
            '--output-dir', str(live / 'classified')], env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
    return dict(status='PASS', accession=row['assembly_accession'], assemblies_downloaded=1,
                provenance_sha256=sha(provenance), network_blocked_during_reuse=True,
                full_database_built=False)


def run_tests(output, environment):
    # Child processes inherit this private environment as well.
    previous = dict(os.environ)
    os.environ.clear()
    os.environ.update(environment)
    sys.path.insert(0, str(ROOT / 'tests'))
    try:
        suite = unittest.defaultTestLoader.loadTestsFromNames(TESTS)
        with (output / 'unit.console.log').open('w') as log:
            result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
        if not result.wasSuccessful() or result.skipped:
            raise AssertionError('Regression checks failed or were skipped; see unit.console.log')
        return dict(status='PASS', tests=result.testsRun, skipped=len(result.skipped))
    finally:
        os.environ.clear()
        os.environ.update(previous)


def make_review(output):
    destination = output.parent / 'nonviral_setup_linux_review.tar.gz'
    if destination.exists():
        raise FileExistsError('Review archive already exists')
    selected = []
    # Never traverse the Conda environment/package caches to assemble diagnostics.
    candidates = list(output.glob('*')) + list((output / 'native').glob('*'))
    candidates += list((output / 'receipts').glob('*'))
    candidates += list((output / 'live').rglob('*'))
    candidates += list((output / 'native/work').rglob('*'))
    for path in candidates:
        name = path.relative_to(output).as_posix()
        if not path.is_file() or path.is_symlink():
            continue
        if (name in ('report.json', 'unit.console.log', 'environment.explicit.txt') or name.startswith('receipts/')
            or (name.startswith('native/') and path.parent == output / 'native' and path.suffix in ('.log', '.tsv', '.json', '.txt', '.config'))
            or (name.startswith('live/') and path.name in ('one_assembly.tsv', 'download_provenance.json',
                  'download_summary.json', 'vicat_nonviral_classification_summary.tsv'))
            or (name.startswith('live/') and path.suffix == '.log')
            or (name.startswith('native/work/') and path.name in ('.command.sh', '.command.run', '.command.log', '.command.err', '.exitcode'))):
            selected.append(path)
    with tarfile.open(destination, 'w:gz') as archive:
        for path in sorted(selected):
            archive.add(path, arcname='validation_results/' + path.relative_to(output).as_posix(), recursive=False)
    return destination


def run(args):
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError('Use a new validation output directory')
    output.mkdir(parents=True)
    report = dict(status='FAIL', platform=platform.platform(), expensive_searches_launched=False,
                  student_runs_accessed=False, full_480_assembly_installation_tested=False)
    try:
        if platform.system() != 'Linux' or sys.version_info < (3, 10):
            raise ValueError('Linux and Python 3.10+ are required for the native gate')
        manifest = verify_package()
        report['package_manifest_sha256'] = sha(ROOT / 'NONVIRAL_PACKAGE.json')
        report['packaged_files_verified'] = len(manifest['files'])
        for name in ('java', 'conda', 'mamba', 'bash', 'flock', 'sha256sum'):
            if not shutil.which(name):
                raise ValueError('Required prerequisite missing from PATH: ' + name)
        nextflow = args.nextflow.resolve()
        if not nextflow.is_file():
            raise ValueError('Nextflow launcher not found: ' + str(nextflow))
        environment = private_environment(output, args.nextflow_home, nextflow)
        print('START regression/default-selection checks', flush=True)
        report['regressions'] = run_tests(output, environment)
        print('PASS regression/default-selection checks', flush=True)
        report['live_ncbi'] = live_check(output, environment)
        print('PASS one real NCBI assembly download/classification and offline reuse', flush=True)
        report['native_setup'] = native_check(output, nextflow, environment)
        report['status'] = 'PASS'
    except Exception:
        report['error'] = traceback.format_exc()
        print(report['error'], file=sys.stderr)
    finally:
        write_json(output / 'report.json', report)
        archive = make_review(output)
        print('Final status: ' + report['status'], flush=True)
        print('Upload this review file: ' + str(archive), flush=True)
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--nextflow', required=True, type=Path)
    parser.add_argument('--nextflow-home', required=True, type=Path)
    raise SystemExit(run(parser.parse_args()))
