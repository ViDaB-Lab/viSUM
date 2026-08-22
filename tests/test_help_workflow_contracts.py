from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class HelpWorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        cls.config = (ROOT / "visum.config").read_text(encoding="utf-8")

    def test_help_is_declared_and_checked_before_pipeline_validation(self) -> None:
        self.assertIn("help         = false", self.config)
        help_check = self.workflow.index("if( helpRequested )")
        self.assertLess(help_check, self.workflow.index("if( params.threads != null )"))
        self.assertLess(help_check, self.workflow.index("def singleMode = params.input != null"))
        self.assertIn("println visumHelp()\n        return", self.workflow)

    def test_help_documents_primary_interfaces_and_every_program_switch(self) -> None:
        help_start = self.workflow.index("def visumHelp()")
        workflow_start = self.workflow.index("workflow {", help_start)
        help_text = self.workflow[help_start:workflow_start]

        for option in (
            "--help", "--input", "--prefix", "--type", "--prefix_many",
            "--indir", "--outdir", "--max_cpus", "--max_memory",
            "--harmonizer_audit", "--ictv_csv",
        ):
            self.assertIn(option, help_text)

        for switch in (
            "run_genomad", "run_virsorter2", "run_cenotetaker3", "run_deep6",
            "run_deepmicroclass2", "run_virbot", "run_gianthunter", "run_vicat",
            "run_checkv", "run_tesorter", "run_vitap", "run_vcontact3",
        ):
            self.assertIn(f"--{switch}", help_text)


if __name__ == "__main__":
    unittest.main()
