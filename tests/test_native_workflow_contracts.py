"""Contracts for the actual entry point and its offline native test fixture."""

import json
from pathlib import Path
import re
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from cache_reuse_fixture.full_workflow import FullWorkflowFixture


class NativeWorkflowContracts(unittest.TestCase):
    def test_strict_syntax_uses_qualified_serializer_without_groovy_imports(self):
        workflow = (ROOT / 'visum_nextflow.nf').read_text()
        self.assertNotRegex(workflow, r'(?m)^\s*import\s')
        self.assertIn('groovy.json.JsonOutput.toJson(spec)', workflow)
        self.assertIn('java.util.Base64.getEncoder()', workflow)

    def test_analysis_imports_no_global_preparation_modules(self):
        workflow = (ROOT / 'visum_nextflow.nf').read_text()
        setup = (ROOT / 'workflows/database_setup.nf').read_text()
        preparations = re.findall(r'include \{ (PREPARE_\w+) \}', setup)
        self.assertEqual(len(preparations), 12)
        for process in preparations:
            self.assertNotIn('include { ' + process + ' }', workflow)
            self.assertNotIn(process + '(', workflow)

    def test_preflight_is_a_read_only_noncacheable_barrier(self):
        module = (ROOT / 'modules/local/validate_prepared_assets.nf').read_text()
        workflow = (ROOT / 'visum_nextflow.nf').read_text()
        self.assertIn('cache false', module)
        self.assertIn('val samples, emit: samples', module)
        self.assertIn('ch_samples.toList()', workflow)
        self.assertIn('NORMALIZE_FASTA(ch_ready_samples)', workflow)
        self.assertNotIn('asset_validation.json', workflow)
        for forbidden in ('wget', 'curl', 'download', 'install.lock'):
            self.assertNotIn(forbidden, module)

    def test_analysis_inputs_exclude_preparation_reports(self):
        for name in ('run_deep6', 'run_deepmicroclass2', 'run_virbot', 'run_gianthunter',
                     'run_checkv', 'run_vitap', 'run_vcontact3', 'run_vicat_diamond',
                     'run_vicat_nonviral_diamond', 'prepare_vicat_reference_subset', 'prepare_vicat_hits',
                     'standardize_vitap'):
            module = (ROOT / 'modules/local' / (name + '.nf')).read_text()
            inputs = module.split('    input:')[1].split('    output:')[0]
            self.assertNotRegex(inputs, r'path\(\w+(?:database|installation)_metadata\)')

    def test_all_imported_analysis_inputs_exclude_preparation_metadata(self):
        workflow = (ROOT / 'visum_nextflow.nf').read_text()
        modules = re.findall(r"from './modules/local/(\w+)'", workflow)
        self.assertGreater(len(modules), 20)
        for name in modules:
            with self.subTest(module=name):
                module = (ROOT / 'modules/local' / (name + '.nf')).read_text()
                inputs = module.split('    input:')[1].split('    output:')[0]
                self.assertNotRegex(inputs, r'\w+(?:database|installation|setup)_metadata')

    def test_shipped_config_is_autodiscovered_and_wrapper_is_only_forwarding(self):
        self.assertIn("includeConfig 'visum.config'", (ROOT / 'nextflow.config').read_text())
        launcher = (ROOT / 'visum').read_text()
        self.assertIn('exec "$nextflow_bin" run "$ROOT_DIR/visum_nextflow.nf" "$@"', launcher)
        for forbidden in ('taskset', 'NXF_OPTS', 'VISUM_LOCAL_CPU_CAP', '-resume', 'MAX_CPUS'):
            self.assertNotIn(forbidden, launcher)

    def test_full_fixture_uses_actual_workflow_and_typed_native_parameters(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = FullWorkflowFixture(Path(temporary))
            command = fixture.commandline()
            self.assertEqual(command[command.index('run') + 1], str(fixture.project / 'visum_nextflow.nf'))
            self.assertNotIn('-C', command)
            payload = json.loads(Path(command[command.index('-params-file') + 1]).read_text())
            self.assertIs(payload['run_genomad'], True)
            self.assertIs(payload['setup'], False)
            self.assertTrue(all(payload['run_' + tool] is False for tool in ('deep6', 'vicat', 'checkv', 'virsorter2')))
            self.assertIn('prefix_many', payload)
            setup = fixture.commandline(setup=True)
            payload = json.loads(Path(setup[setup.index('-params-file') + 1]).read_text())
            self.assertIs(payload['setup'], True)
            self.assertNotIn('prefix_many', payload)
            for folder in ('modules', 'workflows', 'assets'):
                self.assertTrue((fixture.project / folder).is_dir())

    def test_modified_analysis_shell_templates_parse(self):
        bash = str(Path('C:/Program Files/Git/bin/bash.exe')) if os.name == 'nt' else shutil.which('bash')
        if not bash or not Path(bash).exists():
            self.skipTest('Bash unavailable')
        for name in ('run_deep6', 'run_deepmicroclass2', 'run_virbot', 'run_gianthunter', 'run_checkv',
                     'run_vitap', 'run_vcontact3', 'run_vicat_diamond', 'run_vicat_nonviral_diamond',
                     'prepare_vicat_hits', 'prepare_vicat_reference_subset', 'validate_prepared_assets',
                     'standardize_vitap'):
            with self.subTest(module=name):
                module = (ROOT / 'modules/local' / (name + '.nf')).read_text()
                script = module.split('"""', 1)[1].rsplit('"""', 1)[0]
                script = re.sub(r'(?<!\\)\$\{[^}]+\}', 'fixture', script).replace('\\$', '$')
                result = subprocess.run([bash, '-n'], input=script, text=True, capture_output=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_asset_root_defaults_are_resolved_after_overrides(self):
        setup = (ROOT / 'workflows/database_setup.nf').read_text()
        workflow = (ROOT / 'visum_nextflow.nf').read_text()
        config = (ROOT / 'visum.config').read_text()
        for option in ('vitap_dir', 'vcontact3_dir', 'vicat_managed_db', 'vicat_nonviral_db'):
            self.assertRegex(config, option + r'\s*= null')
            self.assertIn('params.' + option + ' ?:', workflow)
            self.assertIn('params.' + option + ' ?:', setup)


if __name__ == '__main__':
    unittest.main(verbosity=2)
