from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class TEsorterAnalysisWorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        cls.module = (ROOT / "modules/local/run_tesorter.nf").read_text(
            encoding="utf-8"
        )
        cls.config = (ROOT / "visum.config").read_text(encoding="utf-8")

    def test_analysis_consumes_refined_candidates_for_both_input_types(self) -> None:
        self.assertIn("include { RUN_TESORTER }", self.workflow)
        self.assertIn(
            "RUN_TESORTER(REFINE_PROVIRAL_REGIONS.out.refined)", self.workflow
        )
        self.assertNotIn("if( type == 'dna'", self.module)
        self.assertNotIn("if( type == 'rna'", self.module)

    def test_official_element_mode_command_uses_bundled_rexdb(self) -> None:
        self.assertIn('TEsorter "${refined_fasta}"', self.module)
        self.assertIn("-db rexdb", self.module)
        self.assertIn('-pre "\\$OUTPUT_PREFIX"', self.module)
        self.assertIn('-p "${task.cpus}"', self.module)
        self.assertNotIn("PREPARE_TESORTER_DATABASE", self.workflow)

    def test_switch_resources_and_outputs_are_declared(self) -> None:
        for declaration in (
            "run_tesorter           = true",
            "tesorter_cpus          = 4",
            "tesorter_memory        = '8 GB'",
            "tesorter_time          = '12h'",
        ):
            self.assertIn(declaration, self.config)
        self.assertIn("tesorter.rexdb.cls.tsv", self.module)
        self.assertIn("tesorter.rexdb.dom.tsv", self.module)
        self.assertIn("tesorter.rexdb.dom.gff3", self.module)
        self.assertIn("tesorter_run_metadata.tsv", self.module)
        self.assertIn("emit: completed", self.module)

    def test_taxonomy_and_harmonizer_wait_for_tesorter_when_enabled(self) -> None:
        self.assertIn(".join(RUN_TESORTER.out.completed)", self.workflow)
        self.assertIn(
            "RUN_VITAP(\n            ch_refined_for_taxonomy,",
            self.workflow,
        )
        self.assertIn(
            "RUN_VCONTACT3(\n            ch_refined_for_taxonomy,",
            self.workflow,
        )
        self.assertIn(
            "ch_harmony_inputs = ch_refined_for_taxonomy",
            self.workflow,
        )


if __name__ == "__main__":
    unittest.main()
