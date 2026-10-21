"""Disposable native-Nextflow fixture. No real databases, installs or server IO."""

from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys

SOURCE = Path(__file__).resolve().parent
ROOT = SOURCE.parents[1]
ANALYSIS = ('RUN_GENOMAD', 'RUN_VIRSORTER2', 'RUN_VICAT_DIAMOND', 'RUN_VICAT_NONVIRAL_DIAMOND')
MODULES = ('normalize_fasta', 'run_genomad', 'run_virsorter2', 'run_vicat_diamond',
           'run_vicat_nonviral_diamond', 'viharmony')
PREPARES = {
    'genomad_database': ('genomad_db', 'genomad_database_metadata.tsv'),
    'vicat_database': ('vicat_database', 'vicat_database_setup_metadata.tsv'),
    'vicat_nonviral_database': ('vicat_nonviral_database', 'vicat_nonviral_database_setup_metadata.tsv'),
}


def write(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents, encoding='utf-8', newline='\n')


def prepare_body(original, link, metadata):
    """Replace only a scratch copy's installer body; fail closed on format drift."""
    if original.count('\n    script:') != 1 or not original.rstrip().endswith('}'):
        raise ValueError('Prepare module format changed; review the fixture adapter')
    declarations = original.split('\n    script:', 1)[0]
    return declarations + f'''
    script:
    """
    ln -s "${{database_path}}" {link}
    printf 'fixture\\n' > {metadata}
    """
}}
'''


class CacheFixture:
    def __init__(self, root, nextflow=None):
        self.root = Path(root).resolve()
        self.project = self.root / 'project'
        self.home = self.root / 'run-home'
        self.work = self.root / 'work'
        self.nextflow = nextflow or shutil.which('nextflow')
        self.invocations = 0
        self.home.mkdir(parents=True)
        (self.root / 'events').mkdir()
        for folder in ('bin', 'envs'):
            shutil.copytree(ROOT / folder, self.project / folder)
        shutil.copy2(ROOT / 'visum.config', self.project / 'visum.config')
        for name in (*MODULES, *PREPARES):
            destination = self.project / 'modules/local' / f'{name}.nf'
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / 'modules/local' / f'{name}.nf', destination)
            if name in PREPARES:
                write(destination, prepare_body(destination.read_text(), *PREPARES[name]))
        shutil.copy2(SOURCE / 'cache_reuse.nf', self.project / 'cache_reuse.nf')
        shutil.copy2(SOURCE / 'fake_tool.py', self.project / 'fake_tool.py')
        for tool in ('genomad', 'virsorter', 'diamond', 'mmseqs', 'python3', 'python'):
            wrapper = self.project / 'fake-bin' / tool
            write(wrapper, f'#!/bin/sh\nexec {shlex.quote(sys.executable)} '
                  f'{shlex.quote(str(self.project / "fake_tool.py"))} {tool} "$@"\n')
            wrapper.chmod(0o755)
        data = self.project / 'data'
        for kind in ('dna', 'rna'):
            write(data / f'{kind}.fasta', '>contig1 fixture\n' + 'ACGT' * 300 + '\n')
            write(data / f'sample_{kind}.orfs.faa', '>orf1\nMKKK\n')
        for name in ('orf_map.tsv', 'gate.tsv', 'regions.tsv', 'groups.tsv', 'msl.csv', 'helper.txt'):
            write(data / name, 'fixture\n')
        write(data / 'refined.fasta', '')
        assets = self.project / 'databases'
        for name, marker in (('genomad', 'fixture.db'), ('virsorter2', 'fixture.db'),
                             ('vicat', 'IMGVR5_UViG_representatives.dmnd'), ('nonviral', 'vicat_nonviral.dmnd')):
            write(assets / name / marker, 'fixture\n')
        for name in ('viral-metadata.tsv', 'nonviral-metadata.tsv'):
            write(assets / name, 'fixture\n')
        write(self.project / 'test.config', '''
conda.enabled = false
docker.enabled = false
singularity.enabled = false
dag.enabled = false
report.enabled = false
timeline.enabled = false
trace.enabled = true
trace.fields = 'task_id,hash,name,status,exit,workdir'
executor.cpus = 2
executor.memory = '256 MB'
process {
    executor = 'local'
    cpus = 1
    memory = '64 MB'
    time = '90s'
    errorStrategy = 'terminate'
    maxRetries = 0
    withName: 'NORMALIZE_FASTA|RUN_.*|PREPARE_.*|VIHARMONY' {
        cpus = 1
        memory = '64 MB'
        time = '90s'
        errorStrategy = 'terminate'
        maxRetries = 0
    }
    withName: 'VIHARMONY' {
        cpus = { params.harmonizer_cpus as int }
    }
}
''')

    def environment(self):
        env = dict(os.environ)
        for key in ('NXF_WORK', 'NXF_ASSETS', 'NXF_OPTS', 'NXF_DEFAULT_DSL', 'NXF_EXECUTOR',
                    'NXF_CACHE_DIR', 'NXF_DISABLE_CACHE', 'NXF_PARAMS_FILE', 'JAVA_TOOL_OPTIONS'):
            env.pop(key, None)
        env['VISUM_CACHE_TEST_ROOT'] = str(self.root)
        env['PATH'] = str(self.project / 'fake-bin') + os.pathsep + env['PATH']
        env['NXF_ANSI_LOG'] = 'false'
        env['NXF_OFFLINE'] = 'true'
        env['NXF_OPTS'] = '-Xms128m -Xmx1g -XX:ActiveProcessorCount=2'
        return env

    def commandline(self, *, resume=False, workers=1, inline=False):
        if type(inline) is not bool:
            raise ValueError('Fixture inline flag must be a Boolean, not a string')
        if type(workers) is not int or not 1 <= workers <= 2:
            raise ValueError('Fixture workers must be an integer between one and two')
        # Fixed cwd is the launch-context contract for the future launcher.
        # This is NOT a test of today's production visum launcher.
        self.invocations += 1
        label = f'run-{self.invocations}'
        self.trace = self.root / f'{label}.trace.tsv'
        self.log = self.root / f'{label}.nextflow.log'
        parameters = self.root / f'{label}.params.json'
        write(parameters, json.dumps({
            'outdir': str(self.root / 'results'),
            'max_cpus': 2,
            'harmonizer_cpus': workers,
            'fixture_inline_prepare': inline,
        }, indent=2) + '\n')
        configs = ','.join(str(self.project / name) for name in ('visum.config', 'test.config'))
        command = [self.nextflow, '-C', configs, '-log', str(self.log),
                   'run', str(self.project / getattr(self, 'workflow', 'cache_reuse.nf')),
                   '-work-dir', str(self.work), '-with-trace', str(self.trace), '-ansi-log', 'false',
                   '-params-file', str(parameters)]
        if resume:
            command += ['-resume']
        self.console = self.root / f'{label}.console.log'
        return command

    def start(self, **options):
        if os.name != 'posix' or not self.nextflow:
            raise RuntimeError('Native cache tests require Linux/Unix with Nextflow on PATH')
        command = self.commandline(**options)
        with self.console.open('w') as stream:
            return subprocess.Popen(command, cwd=self.home, env=self.environment(),
                                    stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)

    def wait(self, process, expected=0):
        try:
            code = process.wait(timeout=180)
        except subprocess.TimeoutExpired:
            self.stop(process)
            raise AssertionError(f'Fixture timed out. Logs: {self.console}')
        if code != expected:
            raise AssertionError(f'Nextflow exit {code}, expected {expected}:\n{self.console.read_text()}')
        return self.rows()

    @staticmethod
    def stop(process):
        if process.poll() is not None:
            return
        os.killpg(process.pid, signal.SIGINT)
        try:
            process.wait(timeout=20)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=5)

    def rows(self):
        with self.trace.open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def counts(self):
        return Counter((event['tool'], event['sample']) for event in
                       (json.loads(path.read_text()) for path in (self.root / 'events').glob('*.json')))

    def outputs(self):
        return {str(path.relative_to(self.root / 'results')): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (self.root / 'results').rglob('*') if path.is_file()}
