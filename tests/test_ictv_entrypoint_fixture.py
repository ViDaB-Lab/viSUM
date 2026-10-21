"""Local contracts for main-workflow saved-input adapter and strict comparisons.

These tests DO NOT claim to execute Nextflow. The separate native runner does.
"""
from collections import Counter
import gzip
import json
import os
from pathlib import Path
import shlex
import tempfile
import unittest

import entrypoint_reconciliation_fixture as fixture
import run_ictv_entrypoint_validation as runner


class EntrypointFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        original = os.environ.get('ICTV_ENTRYPOINT_ORIGINAL')
        if not original:
            raise ValueError('Set ICTV_ENTRYPOINT_ORIGINAL to the verified original extracted package')
        cls.original = Path(original).resolve()
        fixture.verify_original(cls.original)
        cls.temporary = tempfile.TemporaryDirectory(prefix='visum-entrypoint-contracts-')
        cls.root = Path(cls.temporary.name)
        cls.project, cls.specs, cls.protected = fixture.build(cls.original, cls.root)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_main_entrypoint_config_and_scientific_sources_are_byte_identical(self):
        for name, digest in self.protected.items():
            self.assertEqual(digest, fixture.sha(self.original / name), name)
            self.assertEqual(fixture.sha(self.project / name), digest, name)

    def test_all_imports_resolved_only_harmonizer_is_real_process(self):
        for name, path in fixture.imports((self.project / 'visum_nextflow.nf').read_text()):
            module = (self.project / (path.removeprefix('./') + '.nf')).read_text()
            if name == 'VIHARMONY':
                self.assertIn('process VIHARMONY {', module)
            else:
                self.assertIn('workflow ' + name + ' {', module)
                self.assertNotIn('process ', module)
                for keyword in ('script:', 'shell:', 'exec:', 'stub:'):
                    self.assertNotIn(keyword, module)

    def test_adapters_have_no_command_execution_apis(self):
        source = '\n'.join(fixture.adapters().values())
        for forbidden in ('.execute(', 'ProcessBuilder', 'Runtime.getRuntime', 'bash ', 'python ', 'python3 ', 'curl ', 'wget '):
            self.assertNotIn(forbidden, source)

    def test_setup_fails_closed_instead_of_preparing(self):
        source = (self.project / 'workflows/database_setup.nf').read_text()
        self.assertIn("error 'STOP: Database setup is forbidden", source)
        self.assertNotIn('include {', source)

    def test_all_tool_branches_enabled_and_pinned_defaults_not_overridden(self):
        params = json.loads((self.root / 'params.json').read_text())
        for tool in fixture.TOOLS:
            self.assertTrue(params['run_' + tool])
        self.assertFalse(params['setup'])
        self.assertEqual(params['harmonizer_cpus'], 8)
        self.assertEqual(params['harmonizer_max_rank'], 'species')
        self.assertEqual(params['max_memory'], '64 GB')
        self.assertIsNone(params['harmonizer_time'])
        self.assertNotIn('ictv_csv', params)
        self.assertNotIn('ictv_reconciliation', params)

    def test_modern_missing_context_is_explicitly_empty_not_invented_evidence(self):
        for spec in self.specs:
            self.assertEqual(len(Path(spec['context']).read_text().splitlines()), 1)
            if spec['type'] == 'dna':
                self.assertNotIn(spec['loci'], fixture.expected_evidence(spec))
            else:
                self.assertIn(spec['loci'], fixture.expected_evidence(spec))

    def test_all_saved_evidence_and_normalized_inputs_exist(self):
        for spec in self.specs:
            paths = fixture.expected_evidence(spec)
            self.assertEqual(len(paths), len({Path(p).name for p in paths}))
            self.assertEqual(len(paths), 13 if spec['type'] == 'rna' else 12)
            self.assertTrue(all(Path(p).is_file() for p in paths))
            self.assertEqual(Path(spec['normalized']).name, spec['sample'] + '.normalized.fasta')

    def test_multioutput_adapters_cover_vicat_context_loci_projection_and_groups(self):
        adapters = fixture.adapters()
        for name in ('context', 'loci', 'provirus', 'evidence'):
            self.assertIn(name + ' = out_' + name, adapters['STANDARDIZE_VICAT'])
        self.assertIn('projection = out_projection', adapters['PROJECT_VICAT_REFINED'])
        self.assertIn('groups = out_groups', adapters['STANDARDIZE_VCONTACT3'])
        self.assertIn('completed = out_completed', adapters['RUN_TESORTER'])

    def test_adapter_interfaces_match_real_module_contracts(self):
        manifest = Path(fixture.__file__).resolve().parents[1] / 'ENTRYPOINT_PACKAGE.json'
        if manifest.is_file():
            contracts = json.loads(manifest.read_text())['adapter_interfaces_verified']
        else:
            contracts = fixture.interface_contracts(Path(fixture.__file__).resolve().parents[1])
        self.assertEqual(set(contracts), set(fixture.adapters()))
        self.assertEqual(contracts['RUN_CHECKV']['used_outputs']['results'], 9)
        self.assertEqual(contracts['RUN_VITAP']['inputs'], 3)

    def test_nextflow_global_log_precedes_run_and_automatic_config_loading_is_preserved(self):
        for stage in ('first', 'unchanged_resume', 'missing_manifest'):
            argv = runner.command(Path('/nextflow'), self.root, self.project, stage)
            self.assertLess(argv.index('-log'), argv.index('run'))
            self.assertNotIn('-C', argv)
            self.assertEqual(argv[argv.index('run')+1], str(self.project / 'visum_nextflow.nf'))
            self.assertEqual('-resume' in argv, stage == 'unchanged_resume')

    def fake_command(self, spec):
        argv = ['python3', '/project/bin/run_viharmony.py', '--sample-id', spec['sample'], '--input-type', spec['type']]
        for flag, key in [('normalized-fasta','normalized'), ('header-map','header'), ('discovery-gate','gate'),
                          ('refined-fasta','refined'), ('region-map','regions')]:
            argv.extend(['--'+flag, Path(spec[key]).name])
        argv.extend(['--evidence', *[Path(p).name for p in fixture.expected_evidence(spec)],
            '--vcontact3-groups', Path(spec['groups']).name,
            '--ictv-reconciliation', Path(fixture.BUNDLE).name, '--workers', '8',
            '--ictv-msl', 'ICTV_VMR_MSL41.csv', '--audit-mode', 'full',
            '--vcontact3-min-taxonomy-length', str(spec['minimum_length'])])
        if spec['type'] == 'rna' and spec['pair_floor']:
            argv.append('--rna-pair-homology-floor')
        return argv

    def test_actual_command_verifier_checks_all_input_artifacts(self):
        for spec in self.specs:
            path = self.root / (spec['sample'] + '.command.sh')
            argv = self.fake_command(spec)
            fixture.write(path, shlex.join(argv) + '\n')
            result = fixture.verify_command(path, spec, self.project / fixture.BUNDLE)
            self.assertEqual(result['evidence_files_verified'], len(fixture.expected_evidence(spec)))
            argv.pop(argv.index('--evidence')+1)
            fixture.write(path, shlex.join(argv) + '\n')
            with self.assertRaisesRegex(AssertionError, 'dropped/duplicated/mixed evidence'):
                fixture.verify_command(path, spec, self.project / fixture.BUNDLE)

    def test_project_cannot_be_overwritten(self):
        with self.assertRaisesRegex(ValueError, 'already exists'):
            fixture.build(self.original, self.root)


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def outputs(self, label, audit, fasta=b'>a\nACGT\n'):
        root = self.root / label
        root.mkdir()
        (root/'sample.evidence_audit.tsv.gz').write_bytes(gzip.compress(audit.encode()))
        (root/'sample.final.normalized.fasta').write_bytes(fasta)
        manifest = dict(policy=dict(example='unchanged'), outputs={
            'audit':dict(file='sample.evidence_audit.tsv.gz', sha256=fixture.sha(root/'sample.evidence_audit.tsv.gz')),
            'fasta':dict(file='sample.final.normalized.fasta', sha256=fixture.sha(root/'sample.final.normalized.fasta'))})
        fixture.write(root/'sample.harmonizer_manifest.json', json.dumps(manifest))
        return root

    def test_only_audit_row_column_order_is_normalized(self):
        first=self.outputs('a', 'id\ttool\n1\ta\n2\tb\n')
        second=self.outputs('b', 'tool\tid\nb\t2\na\t1\n')
        self.assertTrue(fixture.compare_outputs(first,second,'sample')['scientific_outputs_identical'])

    def test_duplicate_audit_row_is_not_hidden(self):
        first=self.outputs('a', 'id\ttool\n1\ta\n')
        second=self.outputs('b', 'id\ttool\n1\ta\n1\ta\n')
        with self.assertRaisesRegex(AssertionError,'Scientific output differs'):
            fixture.compare_outputs(first,second,'sample')

    def test_changed_fasta_is_not_hidden(self):
        first=self.outputs('a', 'id\n1\n')
        second=self.outputs('b', 'id\n1\n', b'>a\nACGA\n')
        with self.assertRaisesRegex(AssertionError,'Scientific output differs'):
            fixture.compare_outputs(first,second,'sample')

    def test_changed_policy_is_not_hidden(self):
        first=self.outputs('a', 'id\n1\n')
        second=self.outputs('b', 'id\n1\n')
        manifest=second/'sample.harmonizer_manifest.json'
        value=json.loads(manifest.read_text());value['policy']['example']='changed'
        fixture.write(manifest,json.dumps(value))
        with self.assertRaisesRegex(AssertionError,'Policy/reference'):
            fixture.compare_outputs(first,second,'sample')

    def test_bad_output_checksum_is_not_hidden(self):
        first=self.outputs('a', 'id\n1\n')
        second=self.outputs('b', 'id\n1\n')
        (second/'sample.final.normalized.fasta').write_bytes(b'>bad\nN\n')
        with self.assertRaisesRegex(AssertionError,'checksum'):
            fixture.compare_outputs(first,second,'sample')


if __name__ == '__main__':
    unittest.main(verbosity=2)
