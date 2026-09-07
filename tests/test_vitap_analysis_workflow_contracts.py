from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class VitapAnalysisWorkflowContractTests(unittest.TestCase):
    def test_workflow_runs_vitap_after_provirus_refinement(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")

        self.assertIn("include { RUN_VITAP }", workflow)
        self.assertIn("include { STANDARDIZE_VITAP }", workflow)
        self.assertIn("RUN_VITAP(\n            ch_refined_for_taxonomy", workflow)
        self.assertLess(
            workflow.index("REFINE_PROVIRAL_REGIONS(\n        ch_provirus_refinement_inputs,"),
            workflow.index("RUN_VITAP(\n            ch_refined_for_taxonomy"),
        )
        self.assertNotIn("run_vitap_assignment.py", workflow)
        self.assertIn("STANDARDIZE_VITAP(\n            RUN_VITAP.out.results", workflow)
        self.assertIn("STANDARDIZE_VITAP.out.evidence", workflow)

    def test_database_output_is_a_reusable_value_for_analysis_and_standardization(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        assignment = workflow.index(
            "ch_vitap_database = PREPARE_VITAP_DATABASE.out.database"
        )
        invocation = workflow.index("RUN_VITAP(", assignment)
        database_block = workflow[assignment:invocation]

        self.assertIn("tuple(database, metadata)", database_block)
        self.assertIn(".first()", database_block)
        self.assertIn("if( runVitap ) harmonyOnlyTools << 'vitap'", workflow)

    def test_module_uses_refined_fasta_and_enforces_task_cpu_budget(self) -> None:
        module = (ROOT / "modules/local/run_vitap.nf").read_text(encoding="utf-8")

        self.assertIn("path(refined_fasta)", module)
        self.assertIn('cpus { Math.min(params.vitap_cpus as int, params.max_cpus as int) }', module)
        self.assertIn("VITAP assignment", module)
        self.assertIn('-p "${task.cpus}"', module)
        self.assertNotIn("run_vitap_assignment.py", module)
        self.assertIn('export OMP_NUM_THREADS="${task.cpus}"', module)
        self.assertIn("skipped_no_refined_candidates", module)
        self.assertIn("Genome_ID\\tlineage\\tlineage_score/participation_index", module)
        self.assertIn("tr -d '\\r'", module)

    def test_module_preserves_review_and_provenance_outputs(self) -> None:
        module = (ROOT / "modules/local/run_vitap.nf").read_text(encoding="utf-8")

        self.assertIn('path("${prefix}.vitap_best_determined_lineages.tsv")', module)
        self.assertIn('path("${prefix}.vitap_all_lineages.tsv")', module)
        self.assertIn('path("${prefix}.vitap_uniref90_fallback.tsv")', module)
        self.assertIn("path(provirus_region_map)", module)
        self.assertIn('path("${prefix}.vitap_raw")', module)
        self.assertIn("pattern: '*.vitap_*.*'", module)
        self.assertNotIn('pattern: "${prefix}.vitap_*.*"', module)
        self.assertIn("vitap_include_low_confidence", (ROOT / "visum.config").read_text())

    def test_standardizer_uses_refinement_map_and_vmr_taxonomy(self) -> None:
        module = (ROOT / "modules/local/standardize_vitap.nf").read_text(encoding="utf-8")

        self.assertIn("path(provirus_region_map)", module)
        self.assertIn("--region-map", module)
        self.assertIn("--best-lineages", module)
        self.assertIn("--all-lineages", module)
        self.assertIn("--uniref-fallback", module)
        self.assertIn("--vmr", module)
        self.assertIn('path("${prefix}.vitap_evidence.tsv")', module)
        self.assertIn('path("${prefix}.vitap_standardization_audit.tsv")', module)


if __name__ == "__main__":
    unittest.main()
