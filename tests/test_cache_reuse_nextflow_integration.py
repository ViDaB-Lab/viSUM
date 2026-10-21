"""Real Nextflow cache tests with tiny executables, not real biological analyses.

Direct execution fails if Nextflow is unavailable; discovery reports skips.
Set VISUM_CACHE_TEST_ARTIFACTS to an empty/new directory to retain all logs.
"""

from collections import Counter
from contextlib import chdir
import os
import json
from pathlib import Path
import shutil
import sys
import tempfile
import time
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cache_reuse_fixture.harness import ANALYSIS, CacheFixture

NEXTFLOW = shutil.which('nextflow')
SAMPLES = ('sample_dna', 'sample_rna')
EXPECTED = Counter({(tool, sample): 1 for tool in (*ANALYSIS, 'VIHARMONY') for sample in SAMPLES})


def process_name(row):
    return row['name'].split(' (')[0].split(':')[-1]


@unittest.skipUnless(os.name == 'posix' and NEXTFLOW, 'Requires Linux/Unix and Nextflow on PATH')
class NativeCacheReuseTests(unittest.TestCase):
    def setUp(self):
        artifacts = os.environ.get('VISUM_CACHE_TEST_ARTIFACTS')
        if artifacts:
            parent = Path(artifacts).resolve()
            parent.mkdir(parents=True, exist_ok=True)
            root = Path(tempfile.mkdtemp(prefix=self._testMethodName + '-', dir=parent))
        else:
            temporary = tempfile.TemporaryDirectory(prefix='visum-cache-test-')
            self.addCleanup(temporary.cleanup)
            root = Path(temporary.name)
        self.fixture = CacheFixture(root, NEXTFLOW)
        print(f'\nCache fixture: {root}', flush=True)

    def run_flow(self, **options):
        inline = options.get('inline', False)
        caller = options.pop('caller', Path.cwd())
        with chdir(caller):
            process = self.fixture.start(**options)
        self.addCleanup(self.fixture.stop, process)
        rows = self.fixture.wait(process)
        prepares = [process_name(row) for row in rows if process_name(row).startswith('PREPARE_')]
        expected = ['PREPARE_GENOMAD_DATABASE', 'PREPARE_VICAT_DATABASE', 'PREPARE_VICAT_NONVIRAL_DATABASE'] if inline else []
        self.assertCountEqual(prepares, expected, 'Fixture selected the wrong database-preparation branch')
        self.assertIn(f'CACHE_FIXTURE inline_prepare={str(inline).lower()}', self.fixture.log.read_text())
        return rows

    def assert_cached_searches(self, rows):
        searches = [row for row in rows if process_name(row) in ANALYSIS]
        self.assertEqual(len(searches), len(ANALYSIS) * len(SAMPLES), rows)
        self.assertTrue(all(row['status'] == 'CACHED' for row in searches), searches)

    def test_completed_resume_from_different_caller_reuses_every_task(self):
        caller_a, caller_b = self.fixture.root / 'caller-a', self.fixture.root / 'caller-b'
        caller_a.mkdir()
        caller_b.mkdir()
        first = self.run_flow(caller=caller_a)
        self.assertEqual(len(first), 12, first)
        self.assertTrue(all(row['status'] == 'COMPLETED' for row in first), first)
        self.assertEqual(self.fixture.counts(), EXPECTED)
        outputs = self.fixture.outputs()
        self.assertTrue(outputs)
        # Simulate a read-only setup check writing its own report, not search inputs.
        for marker in (self.fixture.project / 'databases').rglob('*'):
            if marker.is_file():
                self.assertTrue(marker.read_bytes())
        (self.fixture.root / 'setup-report.txt').write_text('existing assets validated\n')
        resumed = self.run_flow(resume=True, caller=caller_b)
        self.assertEqual(len(resumed), 12, resumed)
        self.assertTrue(all(row['status'] == 'CACHED' for row in resumed), resumed)
        self.assertEqual(self.fixture.counts(), EXPECTED)
        self.assertEqual(self.fixture.outputs(), outputs)

    def test_interrupt_harmony_then_resume_does_not_repeat_searches(self):
        hold = self.fixture.root / 'hold-harmony'
        hold.touch()
        process = self.fixture.start()
        self.addCleanup(self.fixture.stop, process)
        deadline = time.monotonic() + 90
        while not all(self.fixture.counts()[('VIHARMONY', sample)] == 1 for sample in SAMPLES):
            self.assertIsNone(process.poll(), self.fixture.console.read_text())
            self.assertLess(time.monotonic(), deadline, self.fixture.console.read_text())
            time.sleep(0.1)
        self.fixture.stop(process)
        self.assertNotEqual(process.returncode, 0)
        self.assertEqual(self.fixture.counts(), EXPECTED)
        hold.unlink()
        resumed = self.run_flow(resume=True, workers=2)
        self.assert_cached_searches(resumed)
        counts = self.fixture.counts()
        for tool in ANALYSIS:
            for sample in SAMPLES:
                self.assertEqual(counts[(tool, sample)], 1)
        for sample in SAMPLES:
            self.assertEqual(counts[('VIHARMONY', sample)], 2)
        harmony = [row for row in resumed if process_name(row) == 'VIHARMONY']
        self.assertEqual(len(harmony), 2)
        self.assertTrue(all(row['status'] == 'COMPLETED' for row in harmony), harmony)
        events = [json.loads(path.read_text()) for path in (self.fixture.root / 'events').glob('*.json')]
        for sample in SAMPLES:
            self.assertCountEqual([e['workers'] for e in events if e['tool'] == 'VIHARMONY' and e['sample'] == sample], [1, 2])

    def test_changing_workers_preserves_completed_cache(self):
        self.run_flow()
        outputs = self.fixture.outputs()
        resumed = self.run_flow(resume=True, workers=2)
        self.assert_cached_searches(resumed)
        self.assertEqual(self.fixture.counts(), EXPECTED)
        harmony = [row for row in resumed if process_name(row) == 'VIHARMONY']
        self.assertEqual(len(harmony), 2)
        self.assertTrue(all(row['status'] == 'CACHED' for row in harmony), harmony)
        self.assertEqual(self.fixture.outputs(), outputs)

    def test_missing_cached_output_reexecutes_affected_search(self):
        first = self.run_flow()
        row = next(r for r in first if process_name(r) == 'RUN_GENOMAD' and '(sample_dna)' in r['name'])
        target = Path(row['workdir']) / 'sample_dna_virus_summary.tsv'
        target.resolve().relative_to(self.fixture.work.resolve())  # Fixture-owned target only.
        target.unlink()
        resumed = self.run_flow(resume=True)
        counts = self.fixture.counts()
        self.assertEqual(counts[('RUN_GENOMAD', 'sample_dna')], 2)
        for tool in ANALYSIS:
            for sample in SAMPLES:
                if (tool, sample) != ('RUN_GENOMAD', 'sample_dna'):
                    self.assertEqual(counts[(tool, sample)], 1)
        rerun = [r for r in resumed if process_name(r) == 'RUN_GENOMAD' and '(sample_dna)' in r['name']]
        self.assertEqual(len(rerun), 1)
        self.assertEqual(rerun[0]['status'], 'COMPLETED')

    def test_inline_noncached_preparation_is_a_cache_invalidation_control(self):
        first = self.run_flow(inline=True)
        resumed = self.run_flow(resume=True, inline=True)
        # Characterization of the old wiring, NOT the desired alpha behavior.
        # Fresh preparation work directories change staged database paths.
        for tool in ('RUN_VICAT_DIAMOND', 'RUN_VICAT_NONVIRAL_DIAMOND'):
            rows = [r for r in resumed if process_name(r) == tool]
            self.assertEqual(len(rows), 2)
            self.assertTrue(all(r['status'] == 'COMPLETED' for r in rows), rows)
            for sample in SAMPLES:
                self.assertEqual(self.fixture.counts()[(tool, sample)], 2)
        old = {process_name(r): r['workdir'] for r in first if process_name(r).startswith('PREPARE_')}
        new = {process_name(r): r['workdir'] for r in resumed if process_name(r).startswith('PREPARE_')}
        self.assertEqual(len(old), 3)
        self.assertEqual(old.keys(), new.keys())
        self.assertTrue(all(old[name] != new[name] for name in old), (old, new))
        for sample in SAMPLES:
            self.assertEqual(self.fixture.counts()[('RUN_VIRSORTER2', sample)], 1)


if __name__ == '__main__':
    if os.name != 'posix' or not NEXTFLOW:
        raise SystemExit('ERROR: Native cache tests were NOT RUN. Use Linux/Unix with Nextflow on PATH.')
    unittest.main(verbosity=2)
