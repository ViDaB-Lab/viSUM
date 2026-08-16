from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import standardize_vcontact3


def write_table(
    path: Path, columns: list[str], rows: list[dict[str, str]], delimiter: str
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter=delimiter, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def assignment_row(sequence_id: str, **updates: str) -> dict[str, str]:
    row = {column: "" for column in sorted(standardize_vcontact3.ASSIGNMENT_REQUIRED)}
    row.update(
        {
            "Genome": sequence_id,
            "GenomeName": sequence_id,
            "Proteins": "5",
            "Reference": "False",
            "Size_Kb": "1.0",
        }
    )
    row.update(updates)
    return row


class StandardizeVcontact3Tests(unittest.TestCase):
    def build_inputs(self, directory: Path) -> dict[str, Path]:
        paths = {
            "regions": directory / "regions.tsv",
            "metadata": directory / "metadata.tsv",
            "prokaryotes": directory
            / "sample.vcontact3_prokaryotes_final_assignments.csv",
            "eukaryotes": directory
            / "sample.vcontact3_eukaryotes_final_assignments.csv",
            "evidence": directory / "evidence.tsv",
            "groups": directory / "groups.tsv",
        }
        region_columns = [
            "sample_id",
            "input_type",
            "sequence_id",
            "parent_sequence_id",
            "record_type",
            "coordinates",
            "refined_length",
        ]
        region_rows = [
            {
                "sample_id": "sample",
                "input_type": "rna",
                "sequence_id": "sample__c000001",
                "parent_sequence_id": "",
                "record_type": "input_contig",
                "coordinates": "",
                "refined_length": "1000",
            },
            {
                "sample_id": "sample",
                "input_type": "rna",
                "sequence_id": "sample__c000002|provirus_101_700",
                "parent_sequence_id": "sample__c000002",
                "record_type": "provirus",
                "coordinates": "101-700",
                "refined_length": "600",
            },
            {
                "sample_id": "sample",
                "input_type": "rna",
                "sequence_id": "sample__c000003",
                "parent_sequence_id": "",
                "record_type": "input_contig",
                "coordinates": "",
                "refined_length": "800",
            },
            {
                "sample_id": "sample",
                "input_type": "rna",
                "sequence_id": "sample__c000004",
                "parent_sequence_id": "",
                "record_type": "input_contig",
                "coordinates": "",
                "refined_length": "900",
            },
        ]
        write_table(paths["regions"], region_columns, region_rows, "\t")

        metadata_columns = [
            "sample_id",
            "input_type",
            "database_domain",
            "refined_sequence_count",
            "assignment_row_count",
            "candidate_assignment_row_count",
            "clustered_candidate_count",
            "threads",
            "vcontact3_version",
            "database_path",
            "database_version",
            "run_status",
            "input_fasta",
        ]
        metadata_rows = [
            {
                "sample_id": "sample",
                "input_type": "rna",
                "database_domain": domain,
                "refined_sequence_count": "4",
                "assignment_row_count": "10",
                "candidate_assignment_row_count": "4",
                "clustered_candidate_count": "3",
                "threads": "8",
                "vcontact3_version": "3.2.4",
                "database_path": "/db/232",
                "database_version": "232",
                "run_status": "completed",
                "input_fasta": "sample.refined_candidates.fasta",
            }
            for domain in ("prokaryotes", "eukaryotes")
        ]
        write_table(paths["metadata"], metadata_columns, metadata_rows, "\t")

        eukaryotes = [
            assignment_row(
                "sample__c000001",
                realm_prediction="Monodnaviria",
                kingdom_prediction="Shotokuvirae",
                kingdom_evidence="reference",
                kingdom_network_support="1.0",
                phylum_prediction="Cossaviricota",
                phylum_evidence="reference",
                phylum_network_support="1.0",
                class_prediction="Papovaviricetes",
                class_evidence="reference",
                class_network_support="1.0",
                order_prediction="Sepolyvirales",
                order_evidence="reference",
                order_network_support="1.0",
                family_prediction="Polyomaviridae",
                family_evidence="reference",
                family_network_support="1.0",
                genus_prediction="Alphapolyomavirus|Betapolyomavirus",
                genus_evidence="reference",
                genus_network_support="1.0",
                host_domain="vertebrates",
            ),
            assignment_row(
                "sample__c000002|provirus_101_700",
                realm_prediction="default",
                order_prediction="novel_order_12_of_default",
                order_evidence="novel",
                order_network_support="0.0",
            ),
            assignment_row(
                "sample__c000003", Reference="", realm_prediction="singleton"
            ),
            assignment_row(
                "sample__c000004", realm_prediction="Riboviria"
            ),
            assignment_row(
                "REFERENCE_1", Reference="True", realm_prediction="Monodnaviria"
            ),
        ]
        prokaryotes = [
            assignment_row(
                "sample__c000001", Reference="", realm_prediction="singleton"
            ),
            assignment_row(
                "sample__c000002|provirus_101_700",
                realm_prediction="Duplodnaviria",
                class_prediction="Caudoviricetes",
                class_evidence="reference",
                class_network_support="1.0",
                order_prediction="novel_order_105_of_Caudoviricetes",
                order_evidence="novel",
                order_network_support="0.0",
            ),
            assignment_row(
                "sample__c000003", Reference="", realm_prediction="singleton"
            ),
            assignment_row(
                "sample__c000004", Reference="", realm_prediction="singleton"
            ),
        ]
        columns = sorted(standardize_vcontact3.ASSIGNMENT_REQUIRED)
        write_table(paths["eukaryotes"], columns, eukaryotes, ",")
        write_table(paths["prokaryotes"], columns, prokaryotes, ",")
        return paths

    def run_standardizer(self, paths: dict[str, Path]) -> None:
        argv = [
            "standardize_vcontact3.py",
            "--sample-id",
            "sample",
            "--input-type",
            "rna",
            "--region-map",
            str(paths["regions"]),
            "--assignments",
            str(paths["prokaryotes"]),
            "--assignments",
            str(paths["eukaryotes"]),
            "--run-metadata",
            str(paths["metadata"]),
            "--output-evidence",
            str(paths["evidence"]),
            "--output-groups",
            str(paths["groups"]),
        ]
        with patch.object(sys, "argv", argv):
            standardize_vcontact3.main()

    def test_preserves_all_groups_and_writes_only_qualified_taxonomy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            self.run_standardizer(paths)
            evidence = read_tsv(paths["evidence"])
            groups = read_tsv(paths["groups"])

            self.assertEqual(len(groups), 8)
            self.assertEqual(len(evidence), 3)
            self.assertEqual(
                {row["sequence_id"] for row in evidence},
                {
                    "sample__c000001",
                    "sample__c000002|provirus_101_700",
                    "sample__c000004",
                },
            )

            polyoma = next(
                row
                for row in evidence
                if row["sequence_id"] == "sample__c000001"
            )
            self.assertEqual(polyoma["f__Family"], "f__Polyomaviridae")
            self.assertEqual(polyoma["g__Genus"], "g__unclassified")
            self.assertEqual(polyoma["classification_rank"], "family")
            self.assertEqual(polyoma["score_type"], "vcontact3_network_support")
            self.assertEqual(polyoma["vcontact3_database_domain"], "eukaryotes")

            provirus = next(
                row
                for row in evidence
                if row["sequence_id"] == "sample__c000002|provirus_101_700"
            )
            self.assertEqual(provirus["parent_sequence_id"], "sample__c000002")
            self.assertEqual(provirus["coordinates"], "101-700")
            self.assertEqual(provirus["c__Class"], "c__Caudoviricetes")
            self.assertEqual(provirus["o__Order"], "o__unclassified")

            realm_only = next(
                row
                for row in evidence
                if row["sequence_id"] == "sample__c000004"
            )
            self.assertEqual(realm_only["r__Realm"], "r__Riboviria")
            self.assertEqual(realm_only["score"], "")
            self.assertEqual(
                realm_only["vcontact3_assignment_method"], "realm_only"
            )

            novel = next(
                row
                for row in groups
                if row["sequence_id"] == "sample__c000002|provirus_101_700"
                and row["database_domain"] == "prokaryotes"
            )
            self.assertEqual(
                novel["order_prediction"],
                "novel_order_105_of_Caudoviricetes",
            )
            self.assertEqual(novel["order_evidence"], "novel")
            self.assertEqual(novel["novel_group_ranks"], "order")
            self.assertEqual(
                novel["standardization_decision"],
                "reference_supported_taxonomy",
            )

            singleton = next(
                row
                for row in groups
                if row["sequence_id"] == "sample__c000003"
                and row["database_domain"] == "eukaryotes"
            )
            self.assertEqual(singleton["candidate_assignment_status"], "singleton")
            self.assertEqual(singleton["standardization_decision"], "singleton")

    def test_missing_candidate_assignment_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            with paths["prokaryotes"].open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            rows = [row for row in rows if row["Genome"] != "sample__c000003"]
            write_table(
                paths["prokaryotes"],
                sorted(standardize_vcontact3.ASSIGNMENT_REQUIRED),
                rows,
                ",",
            )
            with self.assertRaisesRegex(ValueError, "missing 1 refined candidates"):
                self.run_standardizer(paths)

    def test_zero_candidates_writes_header_only_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            write_table(paths["regions"], sorted(standardize_vcontact3.REGION_REQUIRED), [], "\t")
            metadata_rows = read_tsv(paths["metadata"])
            for row in metadata_rows:
                row["refined_sequence_count"] = "0"
                row["run_status"] = "skipped_no_refined_candidates"
            write_table(paths["metadata"], list(metadata_rows[0]), metadata_rows, "\t")
            self.run_standardizer(paths)
            self.assertEqual(read_tsv(paths["evidence"]), [])
            self.assertEqual(read_tsv(paths["groups"]), [])


if __name__ == "__main__":
    unittest.main()
