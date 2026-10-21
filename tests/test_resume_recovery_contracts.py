from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


class NativeResumeContracts(unittest.TestCase):
    def test_resume_does_not_depend_on_prepare_work_directories(self):
        workflow = (ROOT / 'visum_nextflow.nf').read_text()
        config = (ROOT / 'visum.config').read_text()
        self.assertNotIn('def recoveredDatabaseOutput', workflow)
        for option in ('resume_genomad_prepare_workdir', 'resume_vicat_prepare_workdir', 'resume_vicat_nonviral_prepare_workdir'):
            self.assertNotIn(option, config)
            self.assertIn(option, workflow)  # Rejected explicitly, never used to recover inputs.
        self.assertIn('is obsolete in this build', workflow)
        self.assertIn('do not migrate an active old-build run', workflow)
        for tool in ('genomad', 'vicat', 'vicat_nonviral'):
            self.assertIn('Channel.value(file(prepared.' + tool + '.path))', workflow)


if __name__ == '__main__':
    unittest.main()
