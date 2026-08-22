from pathlib import Path
import csv
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "standardize_tesorter.py"


def write_tsv(path: Path, columns: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


class StandardizeTEsorterTests(unittest.TestCase):
    def test_tesorter_safe_refined_id_maps_back_to_canonical_region(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            region_map = root / "regions.tsv"
            classifications = root / "classifications.tsv"
            domains = root / "domains.tsv"
            metadata = root / "metadata.tsv"
            evidence = root / "evidence.tsv"
            audit = root / "audit.tsv"
            canonical_id = "sample__c000001|viral_region_101_900"
            reported_id = "sample__c000001_viral_region_101_900"

            write_tsv(
                region_map,
                [
                    "sample_id", "input_type", "sequence_id",
                    "parent_sequence_id", "record_type", "coordinates",
                    "refined_length",
                ],
                [[
                    "sample", "dna", canonical_id, "sample__c000001",
                    "provirus", "101-900", "800",
                ]],
            )
            write_tsv(
                classifications,
                ["#TE", "Order", "Superfamily", "Clade", "Complete", "Strand", "Domains"],
                [[reported_id, "LTR", "Gypsy", "unknown", "no", "+", "RT|Ty3_gypsy"]],
            )
            write_tsv(
                domains,
                ["#id", "length", "evalue", "coverge", "probability", "score"],
                [[f"{reported_id}|domain", "100", "1e-20", "80", "0.9", "0.7"]],
            )
            write_tsv(
                metadata,
                [
                    "sample_id", "input_type", "te_classification_count",
                    "domain_row_count", "run_status",
                ],
                [["sample", "dna", "1", "1", "completed_with_te_classifications"]],
            )

            subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--sample-id", "sample", "--input-type", "dna",
                    "--region-map", str(region_map),
                    "--classifications", str(classifications),
                    "--domains", str(domains),
                    "--run-metadata", str(metadata),
                    "--output-evidence", str(evidence),
                    "--output-audit", str(audit),
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
            )

            with evidence.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["sequence_id"], canonical_id)
            self.assertEqual(rows[0]["parent_sequence_id"], "sample__c000001")
            self.assertEqual(rows[0]["coordinates"], "101-900")

    def test_direct_complete_second_pass_and_viral_like_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            region_map = root / "regions.tsv"
            classifications = root / "classifications.tsv"
            domains = root / "domains.tsv"
            metadata = root / "metadata.tsv"
            evidence = root / "evidence.tsv"
            audit = root / "audit.tsv"

            write_tsv(
                region_map,
                [
                    "sample_id", "input_type", "sequence_id",
                    "parent_sequence_id", "record_type", "coordinates",
                    "refined_length",
                ],
                [
                    ["sample", "dna", "seq1", "", "input_contig", "", "1000"],
                    ["sample", "dna", "seq2", "", "input_contig", "", "2000"],
                    ["sample", "dna", "seq3", "", "input_contig", "", "3000"],
                ],
            )
            write_tsv(
                classifications,
                ["#TE", "Order", "Superfamily", "Clade", "Complete", "Strand", "Domains"],
                [
                    ["seq1", "LTR", "Gypsy", "unknown", "no", "+", "RT|Ty3_gypsy"],
                    ["seq2", "LTR", "Gypsy", "unknown", "yes", "+", "GAG|Gypsy RT|Ty3_gypsy"],
                    ["seq3", "LTR", "Retrovirus", "unknown", "none", "+", "none"],
                ],
            )
            write_tsv(
                domains,
                ["#id", "length", "evalue", "coverge", "probability", "score"],
                [
                    ["seq1|domain", "100", "1e-20", "80", "0.9", "0.7"],
                    ["seq2|domain1", "100", "1e-30", "90", "0.95", "0.8"],
                    ["seq2|domain2", "80", "1e-10", "70", "0.8", "0.5"],
                ],
            )
            write_tsv(
                metadata,
                [
                    "sample_id", "input_type", "te_classification_count",
                    "domain_row_count", "run_status",
                ],
                [["sample", "dna", "3", "3", "completed_with_te_classifications"]],
            )

            subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--sample-id", "sample", "--input-type", "dna",
                    "--region-map", str(region_map),
                    "--classifications", str(classifications),
                    "--domains", str(domains),
                    "--run-metadata", str(metadata),
                    "--output-evidence", str(evidence),
                    "--output-audit", str(audit),
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
            )

            with evidence.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["evidence_strength"] for row in rows], ["qualified", "strong", "weak"])
            self.assertEqual(rows[0]["assignment_method"], "direct_hmm")
            self.assertEqual(rows[0]["classification"], "retroelement")
            self.assertEqual(rows[1]["domain_count"], "2")
            self.assertEqual(rows[2]["classification"], "viral_like_mobile_element")
            self.assertEqual(rows[2]["assignment_method"], "second_pass")

    def test_empty_success_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            region_map = root / "regions.tsv"
            classifications = root / "classifications.tsv"
            domains = root / "domains.tsv"
            metadata = root / "metadata.tsv"
            evidence = root / "evidence.tsv"
            audit = root / "audit.tsv"
            write_tsv(
                region_map,
                [
                    "sample_id", "input_type", "sequence_id",
                    "parent_sequence_id", "record_type", "coordinates",
                    "refined_length",
                ],
                [["sample", "rna", "seq1", "", "input_contig", "", "1000"]],
            )
            classifications.write_text("", encoding="utf-8")
            domains.write_text("", encoding="utf-8")
            write_tsv(
                metadata,
                [
                    "sample_id", "input_type", "te_classification_count",
                    "domain_row_count", "run_status",
                ],
                [["sample", "rna", "0", "0", "completed_no_te_classifications"]],
            )
            subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--sample-id", "sample", "--input-type", "rna",
                    "--region-map", str(region_map),
                    "--classifications", str(classifications),
                    "--domains", str(domains),
                    "--run-metadata", str(metadata),
                    "--output-evidence", str(evidence),
                    "--output-audit", str(audit),
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(len(evidence.read_text(encoding="utf-8").splitlines()), 1)


if __name__ == "__main__":
    unittest.main()
