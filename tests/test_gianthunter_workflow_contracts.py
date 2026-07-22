import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class GiantHunterWorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = (ROOT / "visum.config").read_text(encoding="utf-8")
        cls.workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        cls.database = (
            ROOT / "modules" / "local" / "gianthunter_database.nf"
        ).read_text(encoding="utf-8")
        cls.analysis = (
            ROOT / "modules" / "local" / "run_gianthunter.nf"
        ).read_text(encoding="utf-8")
        cls.environment = (ROOT / "envs" / "gianthunter.yml").read_text(
            encoding="utf-8"
        )

    def test_switch_is_wired_and_analysis_is_dna_only(self):
        self.assertIn("run_gianthunter            = true", self.config)
        self.assertIn("params.run_gianthunter", self.workflow)
        self.assertIn("type == 'dna'", self.workflow)
        self.assertIn("if( runGianthunter )", self.workflow)

    def test_code_and_database_are_pinned(self):
        self.assertIn(
            "d53c466687b64062f23c52a9c0b9c35c60914b98",
            self.environment,
        )
        self.assertIn(
            "499e8efa89c9f6c0ffab1692b94c590f32b34c5aae1eaf44c36d7f22074e1c63",
            self.config,
        )
        self.assertIn("sha256sum --check", self.database)

    def test_setup_and_analysis_share_one_environment(self):
        conda_line = 'conda "${projectDir}/envs/gianthunter.yml"'
        self.assertIn(conda_line, self.database)
        self.assertIn(conda_line, self.analysis)

    def test_analysis_uses_the_pinned_cli_and_handles_zero_calls(self):
        self.assertIn("--outpth", self.analysis)
        self.assertNotIn('--out "', self.analysis)
        self.assertIn("skipped_no_sequences_meeting_minimum_length", self.analysis)
        self.assertIn("completed_no_reference_protein_hits", self.analysis)
        self.assertIn("completed_no_giant_virus_calls", self.analysis)

    def test_database_validation_covers_runtime_assets(self):
        for required_name in (
            "database.dmnd",
            "RefVirus.dmnd",
            "RefVirus.faa",
            "RefVirus_anno.pkl",
            "transformer.pth",
            "names.csv",
            "nodes.csv",
            "taxid.csv",
        ):
            self.assertIn(required_name, self.database)


if __name__ == "__main__":
    unittest.main()
