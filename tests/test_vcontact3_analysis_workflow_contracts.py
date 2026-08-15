from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class Vcontact3AnalysisWorkflowContractTests(unittest.TestCase):
    def test_workflow_runs_analysis_after_provirus_refinement(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")

        self.assertIn("include { RUN_VCONTACT3 }", workflow)
        self.assertIn("RUN_VCONTACT3(\n            REFINE_PROVIRAL_REGIONS.out.refined", workflow)
        self.assertLess(
            workflow.index("REFINE_PROVIRAL_REGIONS(ch_provirus_refinement_inputs)"),
            workflow.index("RUN_VCONTACT3(\n            REFINE_PROVIRAL_REGIONS.out.refined"),
        )
        self.assertIn("vcontact3DbDomain", workflow)
        self.assertIn("['both', 'prokaryotes', 'eukaryotes']", workflow)

    def test_module_uses_official_command_and_resource_budget(self) -> None:
        module = (ROOT / "modules/local/run_vcontact3.nf").read_text(encoding="utf-8")

        self.assertIn("vcontact3 run", module)
        self.assertIn('--nucleotide "${refined_fasta}"', module)
        self.assertIn('--db-domain "\\$DOMAIN"', module)
        self.assertIn('--pyrodigal-gv', module)
        self.assertNotIn('--virus-only', module)
        self.assertIn('--threads "${task.cpus}"', module)
        self.assertIn('--no-progress', module)
        self.assertIn('DOMAINS=(prokaryotes eukaryotes)', module)
        self.assertIn('cpus { Math.min(params.vcontact3_cpus as int, params.max_cpus as int) }', module)

    def test_module_preserves_domain_specific_review_outputs(self) -> None:
        module = (ROOT / "modules/local/run_vcontact3.nf").read_text(encoding="utf-8")

        self.assertIn('path("${prefix}.vcontact3_*_final_assignments.csv")', module)
        self.assertIn('path("${prefix}.vcontact3_*_performance_metrics.csv")', module)
        self.assertIn('path("${prefix}.vcontact3_run_metadata.tsv")', module)
        self.assertIn("skipped_no_refined_candidates", module)
        self.assertIn("user_assignment_row_count", module)

    def test_config_documents_domain_mode(self) -> None:
        config = (ROOT / "visum.config").read_text(encoding="utf-8")

        self.assertIn("vcontact3_db_domain       = 'both'", config)
        self.assertIn("prokaryotes", config)
        self.assertIn("eukaryotes", config)


if __name__ == "__main__":
    unittest.main()
