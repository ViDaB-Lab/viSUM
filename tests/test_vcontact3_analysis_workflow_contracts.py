from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class Vcontact3AnalysisWorkflowContractTests(unittest.TestCase):
    def test_workflow_runs_analysis_after_provirus_refinement(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")

        self.assertIn("include { RUN_VCONTACT3 }", workflow)
        self.assertIn("include { STANDARDIZE_VCONTACT3 }", workflow)
        self.assertIn("RUN_VCONTACT3(\n            ch_refined_for_taxonomy", workflow)
        self.assertLess(
            workflow.index("REFINE_PROVIRAL_REGIONS(ch_provirus_refinement_inputs)"),
            workflow.index("RUN_VCONTACT3(\n            ch_refined_for_taxonomy"),
        )
        self.assertIn("vcontact3DbDomain", workflow)
        self.assertIn("['both', 'prokaryotes', 'eukaryotes']", workflow)
        self.assertIn("STANDARDIZE_VCONTACT3(RUN_VCONTACT3.out.results)", workflow)

    def test_database_output_is_reusable_for_every_sample(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        assignment = workflow.index(
            "ch_vcontact3_database = PREPARE_VCONTACT3_DATABASE.out.database"
        )
        invocation = workflow.index("RUN_VCONTACT3(", assignment)
        database_block = workflow[assignment:invocation]

        self.assertIn("tuple(database, metadata)", database_block)
        self.assertIn(".first()", database_block)
        self.assertIn("if( runVcontact3 ) harmonyOnlyTools << 'vcontact3'", workflow)

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
        self.assertIn(r'\$DOMAIN_OUTPUT/exports/final_assignments.csv', module)
        self.assertIn(r'\$DOMAIN_OUTPUT/exports/performance_metrics.csv', module)
        self.assertIn('cpus { Math.min(params.vcontact3_cpus as int, params.max_cpus as int) }', module)

    def test_module_preserves_domain_specific_review_outputs(self) -> None:
        module = (ROOT / "modules/local/run_vcontact3.nf").read_text(encoding="utf-8")

        self.assertIn('path("${prefix}.vcontact3_*_final_assignments.csv")', module)
        self.assertIn('path("${prefix}.vcontact3_*_performance_metrics.csv")', module)
        self.assertIn('path("${prefix}.vcontact3_run_metadata.tsv")', module)
        self.assertIn("skipped_no_refined_candidates", module)
        self.assertIn("candidate_assignment_row_count", module)
        self.assertIn("clustered_candidate_count", module)

    def test_standardizer_publishes_formal_evidence_and_group_membership(self) -> None:
        module = (ROOT / "modules/local/standardize_vcontact3.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("standardize_vcontact3.py", module)
        self.assertIn('path("${prefix}.vcontact3_evidence.tsv")', module)
        self.assertIn('path("${prefix}.vcontact3_group_membership.tsv")', module)
        self.assertIn("--assignments", module)
        self.assertIn("--region-map", module)

    def test_config_documents_domain_mode(self) -> None:
        config = (ROOT / "visum.config").read_text(encoding="utf-8")

        self.assertIn("vcontact3_db_domain       = 'both'", config)
        self.assertIn("prokaryotes", config)
        self.assertIn("eukaryotes", config)
        self.assertIn("vcontact3_standardizer_memory", config)
        self.assertIn("vcontact3_standardizer_time", config)


if __name__ == "__main__":
    unittest.main()
