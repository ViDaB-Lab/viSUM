import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import standardize_checkv


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class StandardizeCheckVTests(unittest.TestCase):
    def build_inputs(self, directory: Path) -> dict[str, Path]:
        paths = {
            name: directory / name
            for name in (
                "candidates.fasta",
                "quality.tsv",
                "completeness.tsv",
                "contamination.tsv",
                "complete.tsv",
                "metadata.tsv",
                "evidence.tsv",
            )
        }
        paths["candidates.fasta"].write_text(
            ">sample__c000001\n" + "A" * 100 + "\n"
            ">sample__c000002\n" + "C" * 80 + "\n"
            ">sample__c000003\n" + "G" * 60 + "\n",
            encoding="utf-8",
        )
        quality_rows = [
            {
                "contig_id": "sample__c000001",
                "contig_length": "100",
                "provirus": "Yes",
                "proviral_length": "60",
                "gene_count": "10",
                "viral_genes": "4",
                "host_genes": "3",
                "checkv_quality": "Low-quality",
                "miuvig_quality": "Genome-fragment",
                "completeness": "25.0",
                "completeness_method": "AAI-based",
                "contamination": "40.0",
                "kmer_freq": "1.0",
                "warnings": "",
            },
            {
                "contig_id": "sample__c000002",
                "contig_length": "80",
                "provirus": "No",
                "proviral_length": "NA",
                "gene_count": "5",
                "viral_genes": "2",
                "host_genes": "0",
                "checkv_quality": "High-quality",
                "miuvig_quality": "High-quality",
                "completeness": "95.0",
                "completeness_method": "AAI-based",
                "contamination": "0.0",
                "kmer_freq": "1.0",
                "warnings": "",
            },
            {
                "contig_id": "sample__c000003",
                "contig_length": "60",
                "provirus": "No",
                "proviral_length": "NA",
                "gene_count": "2",
                "viral_genes": "0",
                "host_genes": "0",
                "checkv_quality": "Not-determined",
                "miuvig_quality": "Genome-fragment",
                "completeness": "NA",
                "completeness_method": "NA",
                "contamination": "0.0",
                "kmer_freq": "1.0",
                "warnings": "",
            },
        ]
        write_tsv(
            paths["quality.tsv"],
            list(quality_rows[0]),
            quality_rows,
        )
        completeness_rows = [
            {
                "contig_id": row["contig_id"],
                "contig_length": row["contig_length"],
                "viral_length": (
                    row["proviral_length"] if row["provirus"] == "Yes" else row["contig_length"]
                ),
                "aai_confidence": (
                    "medium" if row["contig_id"] == "sample__c000002" else "low"
                ),
                "aai_id": "75.0" if row["contig_id"] == "sample__c000002" else "35.0",
                "aai_af": "0.80" if row["contig_id"] == "sample__c000002" else "0.20",
                "aai_num_hits": "5" if row["contig_id"] == "sample__c000002" else "1",
                "hmm_completeness_lower": "NA",
                "hmm_completeness_upper": "NA",
                "hmm_num_hits": "0",
            }
            for row in quality_rows
        ]
        write_tsv(
            paths["completeness.tsv"],
            sorted(standardize_checkv.COMPLETENESS_REQUIRED),
            completeness_rows,
        )
        contamination_rows = [
            {
                "contig_id": "sample__c000001",
                "contig_length": "100",
                "provirus": "Yes",
                "proviral_length": "60",
                "region_types": "host,viral,host",
                "region_coords_bp": "1-20,21-80,81-100",
            },
            {
                "contig_id": "sample__c000002",
                "contig_length": "80",
                "provirus": "No",
                "proviral_length": "NA",
                "region_types": "NA",
                "region_coords_bp": "NA",
            },
            {
                "contig_id": "sample__c000003",
                "contig_length": "60",
                "provirus": "No",
                "proviral_length": "NA",
                "region_types": "NA",
                "region_coords_bp": "NA",
            },
        ]
        write_tsv(
            paths["contamination.tsv"],
            sorted(standardize_checkv.CONTAMINATION_REQUIRED),
            contamination_rows,
        )
        write_tsv(
            paths["complete.tsv"],
            sorted(standardize_checkv.COMPLETE_GENOMES_REQUIRED),
            [
                {
                    "contig_id": "sample__c000003",
                    "contig_length": "60",
                    "prediction_type": "DTR",
                    "confidence_level": "not-determined",
                }
            ],
        )
        write_tsv(
            paths["metadata.tsv"],
            sorted(standardize_checkv.METADATA_REQUIRED),
            [
                {
                    "sample_id": "sample",
                    "input_type": "dna",
                    "candidate_sequence_count": "3",
                    "quality_summary_row_count": "3",
                    "determined_quality_count": "2",
                    "provirus_count": "1",
                    "checkv_version": "CheckV, version 1.1.1",
                    "run_status": "completed",
                }
            ],
        )
        return paths

    def run_standardizer(self, paths: dict[str, Path]) -> None:
        argv = [
            "standardize_checkv.py",
            "--sample-id",
            "sample",
            "--input-type",
            "dna",
            "--candidate-fasta",
            str(paths["candidates.fasta"]),
            "--quality-summary",
            str(paths["quality.tsv"]),
            "--completeness",
            str(paths["completeness.tsv"]),
            "--contamination",
            str(paths["contamination.tsv"]),
            "--complete-genomes",
            str(paths["complete.tsv"]),
            "--run-metadata",
            str(paths["metadata.tsv"]),
            "--output",
            str(paths["evidence.tsv"]),
        ]
        with patch.object(sys, "argv", argv):
            standardize_checkv.main()

    def test_emits_gene_context_without_determined_quality(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            self.run_standardizer(paths)
            rows = read_tsv(paths["evidence.tsv"])

            self.assertEqual(len(rows), 3)
            provirus = next(row for row in rows if row["record_type"] == "provirus")
            self.assertEqual(provirus["sequence_id"], "sample__c000001|provirus_21_80")
            self.assertEqual(provirus["parent_sequence_id"], "sample__c000001")
            self.assertEqual(provirus["coordinates"], "21-80")
            self.assertEqual(provirus["length"], "60")
            self.assertEqual(provirus["topology"], "Provirus")
            self.assertEqual(provirus["d__Domain"], "d__Viruses")
            self.assertEqual(provirus["evidence_strength"], "qualified")
            self.assertEqual(
                provirus["strength_basis"],
                "checkv_provirus_boundary_with_viral_genes",
            )

            high_quality = next(
                row for row in rows if row["sequence_id"] == "sample__c000002"
            )
            self.assertEqual(high_quality["checkv_quality"], "High-quality")
            self.assertEqual(high_quality["score"], "95.0")
            self.assertEqual(high_quality["evidence_strength"], "strong")
            self.assertEqual(
                high_quality["strength_basis"],
                "checkv_high_quality_confident_aai",
            )
            self.assertEqual(high_quality["aai_confidence"], "medium")
            undetermined = next(row for row in rows if row["sequence_id"] == "sample__c000003")
            self.assertEqual(undetermined["classification"], "unclassified")
            self.assertEqual(undetermined["evidence_strength"], "weak")

    def test_undetermined_quality_preserves_host_only_and_mixed_content(self) -> None:
        for viral, host, expected in [(0, 9, "cellular"), (2, 9, "virus"), (2, 0, "virus")]:
            with self.subTest(viral=viral, host=host), tempfile.TemporaryDirectory() as temporary_directory:
                paths = self.build_inputs(Path(temporary_directory))
                quality = read_tsv(paths["quality.tsv"])
                row = next(r for r in quality if r["contig_id"] == "sample__c000003")
                row.update(viral_genes=str(viral), host_genes=str(host), warnings="retained context")
                write_tsv(paths["quality.tsv"], list(quality[0]), quality)
                self.run_standardizer(paths)
                result = next(r for r in read_tsv(paths["evidence.tsv"]) if r["sequence_id"] == "sample__c000003")
                self.assertEqual(result["classification"], expected)
                self.assertEqual(result["viral_genes"], str(viral))
                self.assertEqual(result["host_genes"], str(host))
                self.assertEqual(result["warnings"], "retained context")
                self.assertEqual(result["score"], "")

    def test_rejects_proviral_length_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            rows = read_tsv(paths["quality.tsv"])
            rows[0]["proviral_length"] = "59"
            write_tsv(paths["quality.tsv"], list(rows[0]), rows)
            with self.assertRaises(SystemExit) as raised:
                self.run_standardizer(paths)
            self.assertIn("proviral coordinates disagree", str(raised.exception))

    def test_high_quality_hmm_only_call_remains_qualified(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            rows = read_tsv(paths["completeness.tsv"])
            high_quality = next(
                row for row in rows if row["contig_id"] == "sample__c000002"
            )
            high_quality["aai_confidence"] = "low"
            high_quality["hmm_completeness_lower"] = "92.0"
            high_quality["hmm_completeness_upper"] = "100.0"
            high_quality["hmm_num_hits"] = "3"
            write_tsv(paths["completeness.tsv"], list(rows[0]), rows)

            self.run_standardizer(paths)
            evidence = next(
                row
                for row in read_tsv(paths["evidence.tsv"])
                if row["sequence_id"] == "sample__c000002"
            )
            self.assertEqual(evidence["evidence_strength"], "qualified")
            self.assertEqual(
                evidence["strength_basis"],
                "checkv_high_quality_viral_genes_without_confident_aai",
            )
            self.assertEqual(evidence["hmm_completeness_lower"], "92.0")
            self.assertEqual(evidence["hmm_hit_count"], "3")

    def test_medium_quality_viral_genes_without_host_genes_is_qualified(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            quality_rows = read_tsv(paths["quality.tsv"])
            row = next(
                item for item in quality_rows if item["contig_id"] == "sample__c000002"
            )
            row["checkv_quality"] = "Medium-quality"
            write_tsv(paths["quality.tsv"], list(quality_rows[0]), quality_rows)
            completeness_rows = read_tsv(paths["completeness.tsv"])
            detail = next(
                item
                for item in completeness_rows
                if item["contig_id"] == "sample__c000002"
            )
            detail["aai_confidence"] = "low"
            write_tsv(
                paths["completeness.tsv"], list(completeness_rows[0]), completeness_rows
            )

            self.run_standardizer(paths)
            evidence = next(
                item
                for item in read_tsv(paths["evidence.tsv"])
                if item["sequence_id"] == "sample__c000002"
            )
            self.assertEqual(evidence["classification"], "virus")
            self.assertEqual(evidence["evidence_strength"], "qualified")
            self.assertEqual(
                evidence["strength_basis"],
                "checkv_medium_quality_viral_genes_without_host_genes",
            )

    def test_host_genes_without_viral_genes_is_cellular_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            quality_rows = read_tsv(paths["quality.tsv"])
            row = next(
                item for item in quality_rows if item["contig_id"] == "sample__c000002"
            )
            row["viral_genes"] = "0"
            row["host_genes"] = "3"
            write_tsv(paths["quality.tsv"], list(quality_rows[0]), quality_rows)

            self.run_standardizer(paths)
            evidence = next(
                item
                for item in read_tsv(paths["evidence.tsv"])
                if item["sequence_id"] == "sample__c000002"
            )
            self.assertEqual(evidence["classification"], "cellular")
            self.assertEqual(evidence["evidence_strength"], "qualified")
            self.assertEqual(evidence["d__Domain"], "d__unclassified")

    def test_zero_candidates_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            paths = self.build_inputs(directory)
            paths["candidates.fasta"].write_text("", encoding="utf-8")
            write_tsv(paths["quality.tsv"], sorted(standardize_checkv.QUALITY_REQUIRED), [])
            write_tsv(
                paths["completeness.tsv"],
                sorted(standardize_checkv.COMPLETENESS_REQUIRED),
                [],
            )
            write_tsv(
                paths["contamination.tsv"],
                sorted(standardize_checkv.CONTAMINATION_REQUIRED),
                [],
            )
            write_tsv(
                paths["complete.tsv"],
                sorted(standardize_checkv.COMPLETE_GENOMES_REQUIRED),
                [],
            )
            write_tsv(
                paths["metadata.tsv"],
                sorted(standardize_checkv.METADATA_REQUIRED),
                [
                    {
                        "sample_id": "sample",
                        "input_type": "dna",
                        "candidate_sequence_count": "0",
                        "quality_summary_row_count": "0",
                        "determined_quality_count": "0",
                        "provirus_count": "0",
                        "checkv_version": "CheckV, version 1.1.1",
                        "run_status": "skipped_no_discovery_candidates",
                    }
                ],
            )
            self.run_standardizer(paths)
            self.assertEqual(read_tsv(paths["evidence.tsv"]), [])


if __name__ == "__main__":
    unittest.main()
