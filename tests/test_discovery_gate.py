import argparse
import csv
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "bin"))

import discovery_gate


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_ids(path: Path) -> list[str]:
    return [
        line[1:].strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]


class DiscoveryGateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)
        self.fasta = self.directory / "normalized.fasta"
        self.header_map = self.directory / "header_map.tsv"
        self.candidates = self.directory / "candidates.fasta"
        self.noncandidates = self.directory / "noncandidates.fasta"
        self.audit = self.directory / "audit.tsv"
        self.summary = self.directory / "summary.tsv"

        self.fasta.write_text(
            ">sample__c000001\nAAAA\n"
            ">sample__c000002\nCCCC\n"
            ">sample__c000003\nGGGG\n"
            ">sample__c000004\nTTTT\n",
            encoding="utf-8",
        )
        write_tsv(
            self.header_map,
            [
                "sample_id",
                "sequence_id",
                "parent_sequence_id",
                "record_type",
                "length",
            ],
            [
                {
                    "sample_id": "sample",
                    "sequence_id": f"sample__c{index:06d}",
                    "parent_sequence_id": "",
                    "record_type": "input_contig",
                    "length": "4",
                }
                for index in range(1, 5)
            ],
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def arguments(self, evidence: list[Path]) -> argparse.Namespace:
        return argparse.Namespace(
            sample_id="sample",
            input_type="dna",
            fasta=self.fasta,
            header_map=self.header_map,
            evidence=evidence,
            output_candidates=self.candidates,
            output_noncandidates=self.noncandidates,
            output_audit=self.audit,
            output_summary=self.summary,
        )

    def test_routes_any_qualified_viral_evidence_and_preserves_conflicts(self) -> None:
        viral = self.directory / "viral.tsv"
        cellular = self.directory / "cellular.tsv"
        columns = sorted(discovery_gate.REQUIRED_EVIDENCE_COLUMNS)
        write_tsv(
            viral,
            columns,
            [
                {
                    "sample_id": "sample",
                    "sequence_id": "sample__c000001",
                    "parent_sequence_id": "",
                    "tool": "genomad",
                    "classification": "virus",
                },
                {
                    "sample_id": "sample",
                    "sequence_id": "sample__c000001|provirus_2_3",
                    "parent_sequence_id": "sample__c000001",
                    "tool": "virsorter2",
                    "classification": "virus",
                },
                {
                    "sample_id": "sample",
                    "sequence_id": "sample__c000002",
                    "parent_sequence_id": "",
                    "tool": "cenotetaker3",
                    "classification": "virus",
                },
            ],
        )
        write_tsv(
            cellular,
            columns,
            [
                {
                    "sample_id": "sample",
                    "sequence_id": "sample__c000002",
                    "parent_sequence_id": "",
                    "tool": "deepmicroclass2",
                    "classification": "cellular",
                },
                {
                    "sample_id": "sample",
                    "sequence_id": "sample__c000003",
                    "parent_sequence_id": "",
                    "tool": "deepmicroclass2",
                    "classification": "plasmid",
                },
            ],
        )

        discovery_gate.run(self.arguments([viral, cellular]))

        self.assertEqual(
            fasta_ids(self.candidates),
            ["sample__c000001", "sample__c000002"],
        )
        self.assertEqual(
            fasta_ids(self.noncandidates),
            ["sample__c000003", "sample__c000004"],
        )
        rows = {row["sequence_id"]: row for row in read_tsv(self.audit)}
        self.assertEqual(rows["sample__c000001"]["discovery_status"], "viral")
        self.assertEqual(rows["sample__c000001"]["viral_tool_count"], "2")
        self.assertEqual(rows["sample__c000002"]["discovery_status"], "ambiguous")
        self.assertEqual(rows["sample__c000002"]["advance_to_refinement"], "true")
        self.assertEqual(rows["sample__c000003"]["discovery_status"], "likely_nonviral")
        self.assertEqual(rows["sample__c000004"]["discovery_status"], "unresolved")

    def test_no_evidence_writes_empty_candidates_and_audits_every_contig(self) -> None:
        discovery_gate.run(self.arguments([]))

        self.assertEqual(self.candidates.read_text(encoding="utf-8"), "")
        self.assertEqual(len(read_tsv(self.audit)), 4)
        summary = read_tsv(self.summary)[0]
        self.assertEqual(summary["candidate_sequence_count"], "0")
        self.assertEqual(summary["unresolved_count"], "4")

    def test_rejects_evidence_for_an_unknown_parent(self) -> None:
        evidence = self.directory / "bad.tsv"
        write_tsv(
            evidence,
            sorted(discovery_gate.REQUIRED_EVIDENCE_COLUMNS),
            [
                {
                    "sample_id": "sample",
                    "sequence_id": "child",
                    "parent_sequence_id": "missing",
                    "tool": "virsorter2",
                    "classification": "virus",
                }
            ],
        )

        with self.assertRaisesRegex(ValueError, "absent from the header map"):
            discovery_gate.run(self.arguments([evidence]))

    def test_multiple_regions_from_one_tool_count_as_one_tool(self) -> None:
        status, advance, reason = discovery_gate.classify(
            {"virus": ["virsorter2", "virsorter2"]}
        )
        self.assertEqual(status, "likely_viral")
        self.assertTrue(advance)
        self.assertIn("one tool", reason)

    def test_competitive_vicat_internal_cellular_support_is_one_tool_conflict(self) -> None:
        evidence = self.directory / "vicat.tsv"
        columns = sorted(discovery_gate.REQUIRED_EVIDENCE_COLUMNS | {"cellular_supported_loci"})
        write_tsv(evidence, columns, [{
            "sample_id": "sample", "sequence_id": "sample__c000001",
            "parent_sequence_id": "", "tool": "vicat", "classification": "virus",
            "cellular_supported_loci": "3",
        }])
        discovery_gate.run(self.arguments([evidence]))
        row = read_tsv(self.audit)[0]
        self.assertEqual(row["discovery_status"], "ambiguous")
        self.assertEqual(row["viral_tool_count"], "1")
        self.assertEqual(row["cellular_tool_count"], "1")
        self.assertEqual(row["viral_tools"], "vicat")
        self.assertEqual(row["cellular_tools"], "vicat")
        self.assertEqual(row["advance_to_refinement"], "true")


if __name__ == "__main__":
    unittest.main()
