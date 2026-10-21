"""Portable checks of the offline fixture; these are NOT native cache tests."""

from contextlib import chdir
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cache_reuse_fixture import fake_tool
from cache_reuse_fixture.harness import CacheFixture, MODULES, PREPARES, ROOT, prepare_body


class CacheFixtureTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='visum-cache-fixture-')
        self.addCleanup(temporary.cleanup)
        self.fixture = CacheFixture(Path(temporary.name))

    def test_analysis_modules_and_decision_script_are_copied_unchanged(self):
        for name in MODULES:
            path = Path('modules/local') / f'{name}.nf'
            self.assertEqual((ROOT / path).read_bytes(), (self.fixture.project / path).read_bytes())
        for name in ('run_viharmony.py', 'evidence_schema.py', 'normalize_fasta.py'):
            self.assertEqual((ROOT / 'bin' / name).read_bytes(), (self.fixture.project / 'bin' / name).read_bytes())

    def test_prepare_adapters_preserve_declarations_including_cache_false(self):
        for name in PREPARES:
            original = (ROOT / 'modules/local' / f'{name}.nf').read_text()
            adapted = (self.fixture.project / 'modules/local' / f'{name}.nf').read_text()
            self.assertIn('cache false', adapted)
            self.assertEqual(original.split('\n    script:')[0], adapted.split('\n    script:')[0])
            self.assertNotIn('download', adapted.split('\n    script:')[1])

    def test_prepare_adapter_rejects_unknown_format(self):
        for original in ('process X {}', 'process X {\n    script:\nfoo',
                         'process X {\n    script:\nfoo\n    script:\nbar\n}'):
            with self.assertRaises(ValueError):
                prepare_body(original, 'db', 'metadata')

    def test_fake_harmony_declares_exactly_required_production_outputs(self):
        module = (ROOT / 'modules/local/viharmony.nf').read_text()
        self.assertEqual(len(fake_tool.HARMONY_SUFFIXES), 21)
        for pattern in fake_tool.HARMONY_PATHS:
            self.assertIn('path("' + pattern.replace('{sample}', '${prefix}') + '")', module)
        self.assertNotIn('review_candidates.fasta', module)

    def execute(self, tool, *arguments):
        with chdir(self.fixture.root), patch.dict(os.environ, VISUM_CACHE_TEST_ROOT=str(self.fixture.root)):
            return fake_tool.main([tool, *arguments])

    def test_genomad_version_is_not_an_analysis_execution(self):
        self.assertEqual(self.execute('genomad', '--version'), 0)
        self.assertFalse(self.fixture.counts())
        self.assertEqual(self.execute('genomad', 'end-to-end', '--threads', '1',
                                     'sample_dna.fasta', 'sample_dna.genomad', 'db'), 0)
        self.assertEqual(self.fixture.counts()[('RUN_GENOMAD', 'sample_dna')], 1)
        summary = self.fixture.root / 'sample_dna.genomad/sample_dna_summary'
        self.assertEqual(len(list(summary.iterdir())), 5)

    def test_fake_genomad_tables_pass_real_standardizer_for_dna_and_rna(self):
        for kind in ('dna', 'rna'):
            with self.subTest(input_type=kind):
                prefix = 'schema_' + kind
                fasta = self.fixture.root / (prefix + '.fasta')
                header_map = self.fixture.root / (prefix + '.header_map.tsv')
                result = subprocess.run([sys.executable, '-B', str(ROOT / 'bin/normalize_fasta.py'),
                    '--input', str(self.fixture.project / 'data' / (kind + '.fasta')),
                    '--sample-id', prefix, '--output-fasta', str(fasta), '--output-map', str(header_map)],
                    capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(self.execute('genomad', 'end-to-end', str(fasta),
                    str(self.fixture.root / (prefix + '.genomad')), 'fixture-db'), 0)
                summary = self.fixture.root / (prefix + '.genomad') / (prefix + '_summary')
                metadata = self.fixture.root / (prefix + '.genomad_run_metadata.tsv')
                # Same zero-call metadata contract written by RUN_GENOMAD.
                columns = ['sample_id', 'input_type', 'input_fasta', 'input_sequence_count',
                    'genomad_version', 'threads', 'splits', 'cleanup', 'score_calibration_requested',
                    'score_calibration_applied', 'run_status', 'virus_call_count', 'plasmid_call_count']
                with metadata.open('w', newline='') as handle:
                    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                    writer.writerow(columns)
                    writer.writerow([prefix, kind, fasta.name, 1, 'genomad cache-test-fixture', 1, 1,
                        'true', 'true', 'false', 'completed_no_viruses_detected', 0, 0])
                output = self.fixture.root / (prefix + '.genomad_evidence.tsv')
                command = [sys.executable, '-B', str(ROOT / 'bin/standardize_genomad.py'),
                    '--sample-id', prefix, '--input-type', kind, '--header-map', str(header_map),
                    '--virus-summary', str(summary / (prefix + '_virus_summary.tsv')),
                    '--virus-genes', str(summary / (prefix + '_virus_genes.tsv')),
                    '--plasmid-summary', str(summary / (prefix + '_plasmid_summary.tsv')),
                    '--run-metadata', str(metadata), '--output', str(output)]
                result = subprocess.run(command, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                with output.open() as handle:
                    reader = csv.DictReader(handle, delimiter='\t')
                    self.assertIn('classification', reader.fieldnames)
                    self.assertEqual(list(reader), [])
                # Each of the three tables must retain its schema validation.
                for suffix in fake_tool.GENOMAD_TABLE_COLUMNS:
                    with self.subTest(input_type=kind, malformed_table=suffix):
                        table = summary / (prefix + '_' + suffix)
                        original = table.read_bytes()
                        try:
                            table.write_text('seq_name\n')
                            invalid = subprocess.run(command, capture_output=True, text=True, timeout=30)
                            self.assertNotEqual(invalid.returncode, 0)
                            self.assertIn('Missing columns', invalid.stderr)
                        finally:
                            table.write_bytes(original)

    def test_virsorter_and_both_diamond_executions_are_counted_separately(self):
        self.execute('virsorter', 'run', '-i', 'sample_rna.fasta', '-w', 'sample_rna.virsorter2')
        self.execute('diamond', 'blastp', '--query', 'sample_rna.orfs.faa', '--out', 'sample_rna.vicat_diamond.raw.tsv')
        self.execute('diamond', 'blastp', '--query', 'sample_rna.orfs.faa', '--out', 'sample_rna.vicat_nonviral_diamond.raw.tsv')
        for tool in ('RUN_VIRSORTER2', 'RUN_VICAT_DIAMOND', 'RUN_VICAT_NONVIRAL_DIAMOND'):
            self.assertEqual(self.fixture.counts()[(tool, 'sample_rna')], 1)
        self.assertFalse(list((self.fixture.root / 'events').glob('*.pending')))

    def test_harmony_substitute_is_deterministic_and_records_worker_rerun(self):
        for workers in ('1', '2'):
            self.execute('python3', '/fixture/run_viharmony.py', '--sample-id', 'sample_dna', '--workers', workers)
            for pattern in fake_tool.HARMONY_PATHS:
                self.assertEqual((self.fixture.root / pattern.format(sample='sample_dna')).read_text(), 'cache-test-only:sample_dna\n')
        self.assertEqual(self.fixture.counts()[('VIHARMONY', 'sample_dna')], 2)

    def test_python_wrapper_passes_real_normalization_through(self):
        data = self.fixture.project / 'data'
        self.assertEqual(self.execute('python3', str(self.fixture.project / 'bin/normalize_fasta.py'),
            '--input', str(data / 'dna.fasta'), '--sample-id', 'sample_dna',
            '--output-fasta', 'normalized.fasta', '--output-map', 'map.tsv'), 0)
        self.assertIn('>sample_dna__c000001', (self.fixture.root / 'normalized.fasta').read_text())
        self.assertFalse(self.fixture.counts())

    def test_fixture_environment_does_not_inherit_other_run_cache_settings(self):
        with patch.dict(os.environ, NXF_WORK='/do-not-use', NXF_CACHE_DIR='/do-not-use',
                        NXF_EXECUTOR='awsbatch', JAVA_TOOL_OPTIONS='-wrong'):
            env = self.fixture.environment()
        for name in ('NXF_WORK', 'NXF_CACHE_DIR', 'NXF_EXECUTOR', 'JAVA_TOOL_OPTIONS'):
            self.assertNotIn(name, env)
        self.assertEqual(env['NXF_OFFLINE'], 'true')
        self.assertEqual(env['NXF_OPTS'], '-Xms128m -Xmx1g -XX:ActiveProcessorCount=2')
        self.assertEqual(Path(env['VISUM_CACHE_TEST_ROOT']), self.fixture.root)

    def test_parameters_round_trip_as_booleans_and_integers_not_cli_strings(self):
        for inline in (False, True):
            with self.subTest(inline=inline):
                command = self.fixture.commandline(inline=inline, workers=2)
                self.assertNotIn('--fixture_inline_prepare', command)
                parameters = Path(command[command.index('-params-file') + 1])
                payload = json.loads(parameters.read_text())
                self.assertIs(payload['fixture_inline_prepare'], inline)
                self.assertIs(type(payload['max_cpus']), int)
                self.assertIs(type(payload['harmonizer_cpus']), int)
                self.assertEqual(payload['harmonizer_cpus'], 2)
                self.assertEqual(payload['outdir'], str(self.fixture.root / 'results'))

    def test_commandline_rejects_truthy_false_strings_and_invalid_workers(self):
        for inline in ('false', 'true', None, 0, 1):
            with self.assertRaises(ValueError):
                self.fixture.commandline(inline=inline)
        for workers in ('2', 0, 3, True, 1.5):
            with self.assertRaises(ValueError):
                self.fixture.commandline(workers=workers)

    def test_native_fixture_rejects_untyped_boolean_before_scheduling(self):
        fixture = (ROOT / 'tests/cache_reuse_fixture/cache_reuse.nf').read_text()
        check = fixture.index('params.fixture_inline_prepare instanceof Boolean')
        self.assertLess(check, fixture.index('NORMALIZE_FASTA(Channel.of('))
        self.assertIn('CACHE_FIXTURE inline_prepare=', fixture)

    def test_unknown_tool_is_rejected(self):
        with self.assertRaises(ValueError):
            self.execute('unexpected-tool')

    @unittest.skipIf(os.name == 'posix' and shutil.which('nextflow'), 'Native runtime available')
    def test_direct_native_runner_fails_instead_of_reporting_skips_as_passes(self):
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'tests/test_cache_reuse_nextflow_integration.py')],
                                capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('NOT RUN', result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
