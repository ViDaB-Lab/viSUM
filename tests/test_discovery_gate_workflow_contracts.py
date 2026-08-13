import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class DiscoveryGateWorkflowContractTests(unittest.TestCase):
    def test_workflow_collects_all_standardized_evidence(self) -> None:
        workflow = (REPOSITORY_ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        for tool in (
            "GENOMAD",
            "VIRSORTER2",
            "CENOTETAKER3",
            "DEEP6",
            "DEEPMICROCLASS2",
            "VIRBOT",
            "GIANTHUNTER",
            "VICAT",
        ):
            self.assertIn(f"STANDARDIZE_{tool}.out.evidence", workflow)
        self.assertIn("DISCOVERY_GATE(ch_discovery_gate_inputs)", workflow)
        self.assertIn("groupTuple()", workflow)

    def test_module_publishes_candidate_and_audit_outputs(self) -> None:
        module = (
            REPOSITORY_ROOT / "modules" / "local" / "discovery_gate.nf"
        ).read_text(encoding="utf-8")
        self.assertIn(".discovery_candidates.fasta", module)
        self.assertIn(".discovery_noncandidates.fasta", module)
        self.assertIn(".discovery_gate.tsv", module)
        self.assertIn(".discovery_gate_summary.tsv", module)

    def test_zero_evidence_uses_a_staged_placeholder_without_counting_it(self) -> None:
        workflow = (REPOSITORY_ROOT / "visum_nextflow.nf").read_text(
            encoding="utf-8"
        )
        module = (
            REPOSITORY_ROOT / "modules" / "local" / "discovery_gate.nf"
        ).read_text(encoding="utf-8")
        placeholder = (
            REPOSITORY_ROOT / "assets" / "empty_discovery_evidence.tsv"
        )

        self.assertTrue(placeholder.is_file())
        self.assertIn("empty_discovery_evidence.tsv", workflow)
        self.assertIn("joined.size() < 5 || joined[4] == null", workflow)
        self.assertIn("val(evidence_file_count)", module)
        self.assertIn("evidence_file_count > 0", module)


if __name__ == "__main__":
    unittest.main()
