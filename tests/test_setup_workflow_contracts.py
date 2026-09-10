from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SetupWorkflowTests(unittest.TestCase):
    def test_setup_branches_before_input_and_returns(self):
        text = (ROOT / 'visum_nextflow.nf').read_text()
        start = text.index("if( parseBooleanParameter(params.setup, '--setup') )")
        end = text.index('    def singleMode', start)
        branch = text[start:end]
        self.assertIn('DATABASE_SETUP()', branch)
        self.assertIn('return', branch)
        self.assertIn('does not accept analysis input parameters', branch)
        self.assertLess(end, text.index('    NORMALIZE_FASTA(ch_samples)'))

    def test_only_global_prepare_processes_imported(self):
        text = (ROOT / 'workflows/database_setup.nf').read_text()
        names = re.findall(r'include \{ (\w+) \}', text)
        self.assertEqual(len(names), 12)
        self.assertTrue(all(name.startswith('PREPARE_') for name in names))
        for forbidden in ('NORMALIZE_FASTA', 'PREPARE_VICAT_REFERENCE_SUBSET',
                          'PREPARE_VICAT_HITS', 'RUN_', 'STANDARDIZE_', 'DISCOVERY_GATE', 'PREDICT_'):
            self.assertNotIn(forbidden, text)
        self.assertIn('at least one enabled database-bearing tool', text)
        self.assertIn('TEsorter has no standalone', text)

    def test_tool_selection_and_serial_dependencies(self):
        text = (ROOT / 'workflows/database_setup.nf').read_text()
        self.assertEqual(text.count('ready = PREPARE_'), 11)
        for tool in ('genomad', 'virsorter2', 'cenotetaker3', 'deep6', 'deepmicroclass2',
                     'virbot', 'gianthunter', 'vicat', 'checkv', 'vitap', 'vcontact3'):
            self.assertIn(f'if( flags.{tool} )', text)
        self.assertNotIn('Channel.of(', text)
        self.assertIn('ch_vicat_nonviral_database_request = PREPARE_VICAT_DATABASE.out.database', text)

    def test_validation_is_not_skipped_on_resume(self):
        config = (ROOT / 'visum.config').read_text()
        self.assertIn("withName: 'DATABASE_SETUP:.*'", config)
        self.assertIn('cache = false', config)
        for name in ('genomad_database', 'vcontact3_database'):
            text = (ROOT / f'modules/local/{name}.nf').read_text()
            self.assertIn('validate_reference_database.py', text)
            self.assertIn('--check-mmseqs', text)
            self.assertIn('cache false', text)


if __name__ == '__main__':
    unittest.main()
