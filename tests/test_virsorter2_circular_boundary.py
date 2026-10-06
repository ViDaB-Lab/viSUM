"""Regression coverage for VirSorter2's unrolled circular full boundaries."""

import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))

import discovery_gate
import run_viharmony
import standardize_virsorter2


SAMPLE = "sample"
PARENT = "sample__c000001"
PREDICTION = f"{PARENT}||full"
CORE = "ACGT" * 480 + "ACG"
SEQUENCE = CORE + CORE[:73]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class CircularFullBoundaryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.header_map = self.root / "header_map.tsv"
        self.score = self.root / "score.tsv"
        self.boundary = self.root / "boundary.tsv"
        self.metadata = self.root / "metadata.tsv"
        self.normalized = self.root / "normalized.fasta"
        self.viral = self.root / "viral.fasta"
        self.output = self.root / "evidence.tsv"
        self.audit = self.root / "boundary_audit.tsv"

        write_tsv(self.header_map, [{
            "sample_id": SAMPLE, "sequence_id": PARENT,
            "parent_sequence_id": "", "record_type": "input_contig",
            "length": "1996",
        }])
        write_tsv(self.score, [{
            "seqname": PREDICTION, "dsDNAphage": "0.967",
            "max_score": "0.967", "max_score_group": "dsDNAphage",
            "length": "1923", "hallmark": "0", "viral": "100.000",
            "cellular": "0.000",
        }])
        write_tsv(self.boundary, [{
            "seqname": PARENT, "trim_orf_index_start": "1",
            "trim_orf_index_end": "5", "trim_bp_start": "12",
            "trim_bp_end": "2017", "partial": "0", "hallmark_cnt": "0",
            "shape": "circular", "seqname_new": PREDICTION,
            "final_max_score": "0.967", "final_max_score_group": "dsDNAphage",
        }])
        write_tsv(self.metadata, [{
            "sample_id": SAMPLE, "input_type": "rna",
            "virsorter2_version": "2.2.4", "classifier_groups": "dsDNAphage",
            "min_length": "1500", "min_score": "0.5",
            "run_status": "completed_with_virus_calls", "virus_call_count": "1",
        }])
        self.normalized.write_text(f">{PARENT}\n{SEQUENCE}\n", encoding="utf-8")
        self.viral.write_text(f">{PREDICTION}\n{SEQUENCE}\n", encoding="utf-8")

    def run_adapter(self) -> None:
        argv = [
            "standardize_virsorter2.py", "--sample-id", SAMPLE,
            "--input-type", "rna", "--header-map", str(self.header_map),
            "--score-table", str(self.score), "--boundary-table", str(self.boundary),
            "--run-metadata", str(self.metadata),
            "--normalized-fasta", str(self.normalized),
            "--viral-fasta", str(self.viral),
            "--boundary-audit", str(self.audit), "--output", str(self.output),
        ]
        with patch.object(sys, "argv", argv):
            standardize_virsorter2.main()

    def test_overlapping_circular_full_call_keeps_vote_without_linear_boundary(self) -> None:
        self.assertEqual(len(SEQUENCE), 1996)
        self.assertEqual(SEQUENCE[:73], SEQUENCE[-73:])
        self.run_adapter()

        evidence = read_tsv(self.output)
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0]["record_type"], "input_contig")
        self.assertEqual(evidence[0]["coordinates"], "")
        self.assertEqual(evidence[0]["score"], "0.967")
        self.assertEqual(evidence[0]["evidence_strength"], "strong")
        audit = read_tsv(self.audit)
        self.assertEqual(len(audit), 1)
        self.assertEqual(audit[0]["native_coordinates"], "12-2017")
        self.assertEqual(audit[0]["boundary_status"], "circular_full_boundary_unresolved")
        self.assertEqual(audit[0]["verified_sequence_sha256"], hashlib.sha256(SEQUENCE.encode()).hexdigest())

        headers = discovery_gate.load_header_map(self.header_map, SAMPLE)
        votes = discovery_gate.load_evidence([self.output], SAMPLE, headers)
        self.assertEqual(discovery_gate.classify(votes[PARENT]),
                         ("likely_viral", True, "qualified viral evidence from one tool"))
        self.assertTrue(run_viharmony.evidence_applies_to_final(
            PREDICTION, PARENT, "input_contig", None,
            PARENT, PARENT, "input_contig", None, len(SEQUENCE), {PARENT},
        ))

    def test_wrong_viral_sequence_does_not_inherit_the_strong_vote(self) -> None:
        bad = SEQUENCE[:-1] + ("A" if SEQUENCE[-1] != "A" else "C")
        self.viral.write_text(f">{PREDICTION}\n{bad}\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "does not match its normalized parent"):
            self.run_adapter()
        self.assertFalse(self.output.exists())

    def test_grossly_out_of_range_circular_boundary_still_fails(self) -> None:
        rows = read_tsv(self.boundary)
        rows[0]["trim_bp_end"] = str(2 * len(SEQUENCE) + 1)
        write_tsv(self.boundary, rows)
        with self.assertRaisesRegex(ValueError, "Invalid boundary coordinates"):
            self.run_adapter()
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
