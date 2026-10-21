"""Native launch/setup/resume checks against the actual workflow, with tiny inputs."""

import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cache_reuse_fixture.full_workflow import FullWorkflowFixture

NEXTFLOW = shutil.which('nextflow')


@unittest.skipUnless(os.name == 'posix' and NEXTFLOW, 'Requires Linux/Unix and Nextflow on PATH')
class NativeWorkflowTests(unittest.TestCase):
    def setUp(self):
        parent = os.environ.get('VISUM_CACHE_TEST_ARTIFACTS')
        if parent:
            Path(parent).mkdir(parents=True, exist_ok=True)
            root = Path(tempfile.mkdtemp(prefix=self._testMethodName + '-', dir=parent))
        else:
            temp = tempfile.TemporaryDirectory(prefix='visum-native-workflow-')
            self.addCleanup(temp.cleanup)
            root = Path(temp.name)
        self.fixture = FullWorkflowFixture(root, NEXTFLOW)
        print(f'\nActual-workflow fixture: {root}', flush=True)

    def run_flow(self, **options):
        process = self.fixture.start(**options)
        self.addCleanup(self.fixture.stop, process)
        rows = self.fixture.wait(process)
        self.assertFalse(any(row['name'].startswith('PREPARE_') for row in rows), rows)
        return rows

    def test_real_workflow_resume_survives_separate_repeated_setup(self):
        first = self.run_flow()
        self.assertEqual(len(first), 13, first)
        self.assertTrue(all(row['status'] == 'COMPLETED' for row in first), first)
        self.assertIn('cpus=2; memory=256 MB', self.fixture.log.read_text())
        events = self.fixture.counts()
        outputs = self.fixture.outputs()
        self.assertEqual(len(events), 4)  # Two fake geNomad runs and two real viHARMONY runs.
        for _ in range(2):
            process = self.fixture.start_setup()
            self.addCleanup(self.fixture.stop, process)
            setup_rows = self.fixture.wait(process)
            self.assertEqual([r['name'].split(' (')[0].split(':')[-1] for r in setup_rows], ['PREPARE_GENOMAD_DATABASE'])
            self.assertEqual(self.fixture.counts(), events)
        caller = self.fixture.root / 'different-caller'
        caller.mkdir()
        from contextlib import chdir
        with chdir(caller):
            resumed = self.run_flow(resume=True)
        self.assertEqual(len(resumed), 13, resumed)
        for row in resumed:
            self.assertEqual(row['status'], 'COMPLETED' if row['name'].startswith('VALIDATE_PREPARED_ASSETS') else 'CACHED', row)
        self.assertEqual(self.fixture.counts(), events)
        self.assertEqual(self.fixture.outputs(), outputs)

    def assert_preflight_failure_without_analyses(self, text):
        process = self.fixture.start()
        self.addCleanup(self.fixture.stop, process)
        self.assertNotEqual(process.wait(timeout=180), 0)
        self.assertIn(text, self.fixture.console.read_text())
        self.assertFalse(self.fixture.counts())
        rows = self.fixture.rows()
        self.assertEqual(len(rows), 1, rows)
        self.assertTrue(rows[0]['name'].startswith('VALIDATE_PREPARED_ASSETS'), rows)

    def test_missing_database_fails_before_scientific_tasks(self):
        target = self.fixture.project / 'databases/genomad/genomad_db_h.index'
        target.unlink()  # Only an asset created by this fixture.
        self.assert_preflight_failure_without_analyses('Run Nextflow with --setup')

    def test_duplicate_sample_prefix_fails_before_scientific_tasks(self):
        (self.fixture.project / 'data/samples.csv').write_text('prefix,type,fasta\nsame,dna,dna.fasta\nsame,rna,rna.fasta\n')
        self.assert_preflight_failure_without_analyses('Duplicate sample prefix')

    def test_installing_database_fails_before_scientific_tasks(self):
        marker = self.fixture.project / 'databases/genomad/.visum_installing'
        marker.touch()
        self.assert_preflight_failure_without_analyses('Installation is incomplete or in progress')
        self.assertTrue(marker.exists())


if __name__ == '__main__':
    if os.name != 'posix' or not NEXTFLOW:
        raise SystemExit('ERROR: Actual-workflow native tests were NOT RUN. Use Linux/Unix with Nextflow on PATH.')
    unittest.main(verbosity=2)
