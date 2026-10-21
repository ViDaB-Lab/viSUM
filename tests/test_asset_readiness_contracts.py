"""Shared setup/preflight contracts; no downloads, databases or inference."""

import ast
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
sys.path.insert(0, str(ROOT / 'tests'))
import validate_reference_database as contracts
import validate_prepared_assets as preflight
import test_prepared_assets as fixtures


class AssetReadinessContracts(unittest.TestCase):
    def setUp(self):
        self.fixtures = fixtures.PreparedAssetsTests()
        self.fixtures.setUp()
        self.addCleanup(self.fixtures.doCleanups)

    def cli(self, tool, spec, *extra):
        command = [sys.executable, '-B', str(ROOT / 'bin/validate_reference_database.py'), tool, spec['path']]
        for key, option in (('models', '--models'), ('hmm_version', '--hmm-version'), ('domain', '--domain')):
            if key in spec:
                command += [option, spec[key]]
        return subprocess.run(command + list(extra), capture_output=True, text=True, timeout=30)

    def test_preflight_imports_exact_setup_contract(self):
        self.assertIs(preflight.validate_asset, contracts.validate_asset)
        ast.parse((ROOT / 'bin/validate_reference_database.py').read_text(), feature_version=(3, 6))

    def test_all_twelve_cli_contracts_accept_complete_and_reject_empty_assets(self):
        for tool in contracts.TOOLS:
            with self.subTest(tool=tool):
                spec = self.fixtures.fixture(tool)
                result = self.cli(tool, spec)
                self.assertEqual(result.returncode, 0, result.stderr)
                target = sorted(p for p in Path(spec['path']).rglob('*') if p.is_file())[0]
                target.write_bytes(b'')
                with self.assertRaises((ValueError, OSError, KeyError)):
                    preflight.validate_asset(tool, spec)
                result = self.cli(tool, spec)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('validation failed', result.stderr)

    def test_installing_markers_rejected_for_every_tool_without_mutation(self):
        for tool in contracts.TOOLS:
            with self.subTest(tool=tool):
                spec = self.fixtures.fixture(tool)
                marker = Path(spec['path']) / '.visum_installing'
                marker.touch()
                before = marker.stat().st_mtime_ns
                with self.assertRaisesRegex(ValueError, 'incomplete or in progress'):
                    preflight.validate_asset(tool, spec)
                result = self.cli(tool, spec)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(marker.stat().st_mtime_ns, before)
                self.assertEqual(marker.read_bytes(), b'')
                # Installer-only allowance does not bypass the structural checks.
                contracts.validate_asset(tool, spec, allow_installing=True)
                target = next(p for p in Path(spec['path']).rglob('*') if p.is_file() and p != marker)
                target.write_bytes(b'')
                with self.assertRaises((ValueError, OSError, KeyError)):
                    contracts.validate_asset(tool, spec, allow_installing=True)

    def test_external_deep6_models_installing_marker_rejected(self):
        spec = self.fixtures.fixture('deep6')
        (Path(spec['models']) / '.visum_installing').touch()
        with self.assertRaisesRegex(ValueError, 'incomplete or in progress'):
            preflight.validate_asset('deep6', spec)

    def test_interrupted_managed_installations_cannot_take_reuse_shortcut(self):
        guard = 'if [[ ! -e "\\$INSTALLING_MARKER" && ! -L "\\$INSTALLING_MARKER" ]] && validate_database "\\$DB_DEST"; then'
        for name in ('virsorter2_database', 'cenotetaker3_database', 'gianthunter_database'):
            source = (ROOT / 'modules/local' / (name + '.nf')).read_text()
            self.assertIn(guard, source)
            self.assertLess(source.index('flock 9'), source.index(guard))

    def test_checkv_junk_file_is_not_an_hmm_in_either_interface(self):
        spec = self.fixtures.fixture('checkv')
        hmm = Path(spec['path']) / 'hmm_db/checkv_hmms/test.hmm'
        hmm.rename(hmm.with_suffix('.txt'))
        with self.assertRaisesRegex(ValueError, 'Missing CheckV HMM'):
            preflight.validate_asset('checkv', spec)
        self.assertNotEqual(self.cli('checkv', spec).returncode, 0)

    def test_checkv_quick_check_does_not_open_all_hmms(self):
        spec = self.fixtures.fixture('checkv')
        for number in range(20):
            self.fixtures.asset('checkv/hmm_db/checkv_hmms/{}.hmm'.format(number))
        with patch.object(contracts, 'nonempty', wraps=contracts.nonempty) as read:
            preflight.validate_asset('checkv', spec)
        hmm_calls = [c for c in read.call_args_list if str(c.args[0]).endswith('.hmm')]
        self.assertEqual(len(hmm_calls), 1)

    def test_vs2_empty_environment_history_rejected(self):
        spec = self.fixtures.fixture('virsorter2')
        history = next(Path(spec['path']).glob('conda_envs/**/conda-meta/history'))
        history.write_bytes(b'')
        with self.assertRaisesRegex(ValueError, 'history'):
            preflight.validate_asset('virsorter2', spec)
        self.assertNotEqual(self.cli('virsorter2', spec).returncode, 0)

    def test_all_twelve_asset_checks_leave_files_and_mtimes_unchanged(self):
        specs = {tool: self.fixtures.fixture(tool) for tool in contracts.TOOLS}
        def inventory():
            return {str(p): (p.read_bytes(), p.stat().st_mtime_ns)
                    for p in self.fixtures.root.rglob('*') if p.is_file()}
        before = inventory()
        for tool, spec in specs.items():
            preflight.validate_asset(tool, spec)
        self.assertEqual(inventory(), before)
        report = preflight.validate(self.fixtures.spec(specs))
        self.assertEqual(report['runtime_validation'], 'not-performed')

    def test_all_prepares_delegate_structural_validation_and_shells_parse(self):
        setup = (ROOT / 'workflows/database_setup.nf').read_text()
        modules = re.findall(r"from '../modules/local/(\w+)'", setup)
        self.assertEqual(len(modules), 12)
        bash = 'C:/Program Files/Git/bin/bash.exe' if os.name == 'nt' else shutil.which('bash')
        for name in modules:
            with self.subTest(module=name):
                module = (ROOT / 'modules/local' / (name + '.nf')).read_text()
                self.assertIn('validate_reference_database.py', module)
                if name == 'vcontact3_database':
                    shell = module.split("    '''", 1)[1].rsplit("    '''", 1)[0]
                    shell = re.sub(r'!\{[^}]+\}', 'fixture', shell)
                else:
                    shell = module.split('    """', 1)[1].rsplit('    """', 1)[0]
                    shell = re.sub(r'(?<!\\)\$\{[^}]+\}', 'fixture', shell).replace('\\$', '$')
                if bash and Path(bash).exists():
                    result = subprocess.run([bash, '-n'], input=shell, text=True, capture_output=True, timeout=30)
                    self.assertEqual(result.returncode, 0, result.stderr)

    def test_default_memory_budget_covers_non_build_requests(self):
        config = (ROOT / 'visum.config').read_text()
        default = int(re.search(r"max_memory\s*=\s*'(\d+) GB'", config)[1])
        self.assertEqual(default, 64)
        safe = config.split('  local_safe {', 1)[1].split('\n  }', 1)[0]
        self.assertIn("params.max_memory = '64 GB'", safe)
        for name, amount in re.findall(r"(?m)^\s*(\w+_memory)\s*=\s*'(\d+) GB'", config):
            if '_build_memory' not in name:
                self.assertLessEqual(int(amount), default, name)
        self.assertEqual(json.loads((ROOT / 'docs/example-analysis.json').read_text())['max_memory'], '64 GB')

    def test_deep6_setup_uses_analysis_environment_and_import_smoke(self):
        setup = (ROOT / 'modules/local/deep6_database.nf').read_text()
        analysis = (ROOT / 'modules/local/run_deep6.nf').read_text()
        declaration = 'conda "${projectDir}/envs/deep6.yml"'
        self.assertIn(declaration, setup)
        self.assertIn(declaration, analysis)
        self.assertIn('import tensorflow, keras, h5py, numpy, pandas, sklearn, Bio', setup)
        self.assertIn('imports-passed-not-inference-tested', setup)
        self.assertLess(setup.index('import tensorflow'), setup.index('git clone'))
        env = (ROOT / 'envs/deep6.yml').read_text()
        for pin in ('python=3.6.15', 'tensorflow=2.4.1', 'keras=2.4.3', 'h5py=2.10.0'):
            self.assertIn(pin, env)
        self.assertIn('  - git\n', env)


if __name__ == '__main__':
    unittest.main(verbosity=2)
