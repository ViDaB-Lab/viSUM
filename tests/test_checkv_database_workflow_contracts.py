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

    def test_analysis_runs_only_on_discovery_candidates(self) -> None:
        workflow = (REPOSITORY_ROOT / "visum_nextflow.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("include { RUN_CHECKV }", workflow)
        self.assertIn("DISCOVERY_GATE.out.candidates", workflow)
        self.assertIn("RUN_CHECKV(", workflow)

    def test_analysis_preserves_primary_reports_and_zero_candidate_runs(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/run_checkv.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("checkv end_to_end", module)
        self.assertIn("skipped_no_discovery_candidates", module)
        self.assertIn("quality-summary coverage mismatch", module)
        for expected_output in (
            "checkv_quality_summary.tsv",
            "checkv_completeness.tsv",
            "checkv_contamination.tsv",
            "checkv_complete_genomes.tsv",
            "checkv_proviruses.fna",
            "checkv_run_metadata.tsv",
        ):
            self.assertIn(expected_output, module)


if __name__ == "__main__":
    unittest.main()
