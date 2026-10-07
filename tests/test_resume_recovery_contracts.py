from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ResumeRecoveryContracts(unittest.TestCase):
    def test_recovery_is_opt_in_and_preserves_completed_prepare_paths(self):
        workflow = (ROOT / 'visum_nextflow.nf').read_text(encoding='utf-8')
        config = (ROOT / 'visum.config').read_text(encoding='utf-8')
        for parameter, link, metadata in (
            ('resume_genomad_prepare_workdir', 'genomad_db', 'genomad_database_metadata.tsv'),
            ('resume_vicat_prepare_workdir', 'vicat_database', 'vicat_database_setup_metadata.tsv'),
            ('resume_vicat_nonviral_prepare_workdir', 'vicat_nonviral_database', 'vicat_nonviral_database_setup_metadata.tsv'),
        ):
            self.assertIn(f'{parameter} = null', config)
            self.assertIn(f'params.{parameter}', workflow)
            self.assertIn(f"'{link}', '{metadata}'", workflow)
        self.assertIn("java.nio.file.Files.readString(exitCode).trim() != '0'", workflow)
        self.assertIn('database.toRealPath() != file(expectedDatabase).toRealPath()', workflow)
        self.assertIn('ch_genomad_database = Channel.value(recovered[0])', workflow)
        self.assertIn('ch_vicat_database = Channel.value(recoveredViral)', workflow)
        self.assertIn('ch_vicat_nonviral_database = Channel.value(recoveredNonviral)', workflow)
        self.assertIn('PREPARE_GENOMAD_DATABASE(ch_genomad_database_request)', workflow)
        self.assertIn('PREPARE_VICAT_DATABASE(ch_vicat_database_request)', workflow)
        self.assertIn('PREPARE_VICAT_NONVIRAL_DATABASE(ch_vicat_nonviral_database_request)', workflow)


if __name__ == '__main__':
    unittest.main()
