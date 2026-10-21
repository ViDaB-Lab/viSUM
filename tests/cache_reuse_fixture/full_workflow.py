"""Exercise the actual viSUM entry point; only external classifier binaries are fakes."""

import json
from pathlib import Path
import shutil
import sys

from .harness import CacheFixture, ROOT, write

sys.path.insert(0, str(ROOT / 'tests'))
import test_reference_database_validation as reference_fixtures


class FullWorkflowFixture(CacheFixture):
    def __init__(self, root, nextflow=None):
        super().__init__(root, nextflow)
        for folder in ('modules', 'workflows', 'assets'):
            shutil.copytree(ROOT / folder, self.project / folder, dirs_exist_ok=True)
        for name in ('visum_nextflow.nf', 'nextflow.config'):
            shutil.copy2(ROOT / name, self.project / name)
        self.workflow = 'visum_nextflow.nf'
        builder = reference_fixtures.ReferenceValidationTests()
        builder.root = self.project / 'databases/genomad'
        builder.genomad_fixture()
        write(self.project / 'data/samples.csv', 'prefix,type,fasta\nsample_dna,dna,dna.fasta\nsample_rna,rna,rna.fasta\n')
        config = (self.project / 'test.config').read_text()
        config = config.replace("'NORMALIZE_FASTA|RUN_.*|PREPARE_.*|VIHARMONY'", "'.*'")
        # Test the shipped resource-budget wiring, not a fixture-only budget.
        config = config.replace('executor.cpus = 2\n', '').replace("executor.memory = '256 MB'\n", '')
        write(self.project / 'test.config', config)

    def environment(self):
        env = super().environment()
        env['VISUM_TEST_REAL_HARMONY'] = '1'
        return env

    def commandline(self, *, setup=False, **options):
        command = super().commandline(**options)
        parameters = Path(command[command.index('-params-file') + 1])
        payload = json.loads(parameters.read_text())
        for tool in ('genomad', 'virsorter2', 'cenotetaker3', 'deep6', 'deepmicroclass2',
                     'virbot', 'gianthunter', 'vicat', 'checkv', 'tesorter', 'vitap', 'vcontact3'):
            payload['run_' + tool] = tool == 'genomad'
        payload.update(setup=setup, max_memory='256 MB', genomad_db=str(self.project / 'databases/genomad'))
        if setup:
            payload['outdir'] = str(self.root / 'setup-results')
        else:
            payload.update(prefix_many=str(self.project / 'data/samples.csv'), indir=str(self.project / 'data'))
        write(parameters, json.dumps(payload, indent=2) + '\n')
        # Native automatic nextflow.config loading, plus one test/site override.
        index = command.index('-C')
        command[index:index + 2] = ['-c', str(self.project / 'test.config')]
        return command

    def start_setup(self):
        home, work = self.home, self.work
        self.home = self.root / 'setup-home'
        self.work = self.root / 'setup-work'
        self.home.mkdir(exist_ok=True)
        try:
            return self.start(setup=True)
        finally:
            self.home, self.work = home, work
