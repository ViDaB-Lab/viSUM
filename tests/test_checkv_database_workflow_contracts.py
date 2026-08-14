import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class CheckVDatabaseWorkflowContractTests(unittest.TestCase):
    def test_workflow_wires_optional_checkv_database_setup(self) -> None:
        workflow = (REPOSITORY_ROOT / "visum_nextflow.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("PREPARE_CHECKV_DATABASE", workflow)
        self.assertIn("params.run_checkv", workflow)
        self.assertIn("params.checkv_db", workflow)
        self.assertIn("params.checkv_auto_download", workflow)

    def test_database_module_validates_both_database_components(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/checkv_database.nf").read_text(
            encoding="utf-8"
        )
        for required_path in (
            "genome_db/checkv_reps.faa",
            "genome_db/checkv_reps.fna",
            "genome_db/checkv_reps.tsv",
            "genome_db/checkv_reps.dmnd",
            "hmm_db/checkv_hmms.tsv",
            "hmm_db/genome_lengths.tsv",
            "hmm_db/checkv_hmms",
        ):
            self.assertIn(required_path, module)
        self.assertIn("diamond dbinfo", module)

    def test_database_install_is_staged_locked_and_reusable(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/checkv_database.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("flock 9", module)
        self.assertIn("checkv download_database", module)
        self.assertIn("skipped-existing-managed-database", module)
        self.assertIn("skipped-user-supplied-database", module)


if __name__ == "__main__":
    unittest.main()
