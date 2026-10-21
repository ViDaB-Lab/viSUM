"""Regression for the exact Linux CLI failure and native-only continuation."""
from pathlib import Path
import os
import importlib.util
import subprocess
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_integrated_reconciliation_validation as runner
import resume_integrated_reconciliation_native as resume

WARNING = "Unknown option: -C -- Check the available commands and options and syntax with 'help'\n"


def global_option_boundary(argv):
    """Minimal parser-boundary reproducer, NOT a real Nextflow execution."""
    return WARNING if '-C' in argv[argv.index('run') + 1:] else None


class NativeCommandTests(unittest.TestCase):
    def command(self, stage='first'):
        return runner.native_command(Path('/nextflow'), Path('/native'),
            Path('/native/native.config'), Path('/native/params.json'), stage, Path('/project'))

    def test_original_order_reproduces_exact_reported_warning(self):
        self.assertEqual(global_option_boundary(['/nextflow','run','main.nf','-C','native.config']), WARNING)

    def test_fixed_order_places_both_global_options_before_run(self):
        command = self.command()
        self.assertIsNone(global_option_boundary(command))
        self.assertLess(command.index('-C'), command.index('run'))
        self.assertLess(command.index('-log'), command.index('run'))
        self.assertEqual(command.count('-C'), 1)
        self.assertNotIn('-resume', command)

    def test_only_followup_stages_resume_same_work_directory(self):
        for stage in ('unchanged_resume','manifest_changed','bundle_changed'):
            command = self.command(stage)
            self.assertEqual(command[-1], '-resume')
            self.assertEqual(command[command.index('-work-dir') + 1], str(Path('/native/work')))
            self.assertIsNone(global_option_boundary(command))

    def test_main_workflow_is_not_the_native_module_fixture(self):
        command = self.command()
        self.assertEqual(command[command.index('run') + 1], str(Path('/project/native_validation.nf')))
        self.assertNotIn('visum_nextflow.nf', ' '.join(command))

    def test_continuation_never_invokes_full_replay_or_regression(self):
        source = Path(resume.__file__).read_text()
        self.assertNotIn('runner.run(', source)
        self.assertNotIn('runner.replay(', source)
        self.assertNotIn('unittest.', source)
        self.assertIn('runner.native_check(', source)

    def test_only_module_and_copy_process_are_imported_by_fixture(self):
        root = Path(os.environ.get('ICTV_VALIDATION_ORIGINAL_PACKAGE', str(runner.ROOT)))
        fixture = root / 'native_validation.nf'
        if not fixture.is_file():
            fixture = root / 'tests/integrated_reconciliation_native.nf'
        source = fixture.read_text()
        self.assertEqual(source.count('include {'), 1)
        self.assertIn('include { VIHARMONY }', source)
        for name in ('RUN_GENOMAD', 'RUN_VICAT_DIAMOND', 'RUN_VIRSORTER2', 'PREPARE_GENOMAD_DATABASE'):
            self.assertNotIn(name, source)


if __name__ == '__main__':
    unittest.main(verbosity=2)
