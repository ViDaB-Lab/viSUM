from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ViharmonyWorkflowContractTests(unittest.TestCase):
    def test_workflow_runs_harmonizer_after_refinement_and_taxonomy_tools(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        self.assertIn("include { VIHARMONY }", workflow)
        self.assertGreater(workflow.index("VIHARMONY(ch_harmony_inputs)"), workflow.index("STANDARDIZE_VCONTACT3(RUN_VCONTACT3.out.results)"))
        self.assertIn("ch_harmony_evidence_for_sample", workflow)
        self.assertIn("assets/empty_discovery_evidence.tsv", workflow)

    def test_module_declares_agreed_final_outputs(self) -> None:
        module = (ROOT / "modules" / "local" / "viharmony.nf").read_text(encoding="utf-8")
        for suffix in (
            ".final.normalized.fasta", ".final.original_ids.fasta",
            ".final_metadata.tsv", ".review_queue.tsv",
            ".sequence_disposition.tsv", ".sequence_map.tsv",
            ".harmonizer_manifest.json", ".database_candidates.fasta",
            ".database_candidates.tsv",
        ):
            self.assertIn(suffix, module)
        self.assertIn("--audit-mode", module)
        self.assertIn("--vcontact3-min-taxonomy-length", module)

    def test_config_uses_one_canonical_msl_parameter(self) -> None:
        config = (ROOT / "visum.config").read_text(encoding="utf-8")
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        self.assertIn("harmonizer_audit", config)
        self.assertIn('ictv_csv = "${projectDir}/assets/ICTV_VMR_MSL41.csv"', config)
        self.assertNotIn("harmonizer_ictv_csv", config)
        self.assertIn("def ictvCsv = file(params.ictv_csv)", workflow)
        self.assertNotIn("harmonizerIctvCsv", workflow)
        self.assertIn("vcontact3_min_taxonomy_length", config)


if __name__ == "__main__":
    unittest.main()
