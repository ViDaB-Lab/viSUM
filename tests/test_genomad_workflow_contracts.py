import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class GenomadWorkflowContractTests(unittest.TestCase):
    def test_score_calibration_is_enabled_by_default_and_audited(self) -> None:
        config = (REPOSITORY_ROOT / "visum.config").read_text(encoding="utf-8")
        module = (REPOSITORY_ROOT / "modules/local/run_genomad.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("genomad_score_calibration = true", config)
        self.assertIn("'--genomad_score_calibration'", (
            REPOSITORY_ROOT / "visum_nextflow.nf"
        ).read_text(encoding="utf-8"))
        self.assertIn("--enable-score-calibration", module)
        self.assertIn("score_calibration_requested", module)
        self.assertIn("score_calibration_applied", module)
        self.assertIn('"\\$INPUT_SEQUENCE_COUNT" -ge 1000', module)

    def test_virus_gene_table_is_wired_into_standardization(self) -> None:
        workflow = (REPOSITORY_ROOT / "visum_nextflow.nf").read_text(
            encoding="utf-8"
        )
        module = (
            REPOSITORY_ROOT / "modules/local/standardize_genomad.nf"
        ).read_text(encoding="utf-8")

        self.assertIn(
            "tuple(prefix, type, virusSummary, virusGenes, plasmidSummary, metadata)",
            workflow,
        )
        self.assertIn("path(virus_genes)", module)
        self.assertIn('--virus-genes "${virus_genes}"', module)

    def test_standardizer_exports_uscg_and_strength_fields(self) -> None:
        standardizer = (
            REPOSITORY_ROOT / "bin/standardize_genomad.py"
        ).read_text(encoding="utf-8")

        self.assertIn('"n_uscg"', standardizer)
        self.assertIn('"evidence_strength"', standardizer)
        self.assertIn('"strength_basis"', standardizer)
        self.assertIn('"genomad_conservative_preset"', standardizer)


if __name__ == "__main__":
    unittest.main()
