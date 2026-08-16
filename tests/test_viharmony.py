from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import run_viharmony


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ViharmonyTests(unittest.TestCase):
    def test_combines_parent_evidence_region_taxonomy_and_original_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            normalized = directory / "sample.normalized.fasta"
            refined = directory / "sample.refined.fasta"
            normalized.write_text(">sample__c000001\n" + "A" * 100 + "\n>sample__c000002\n" + "C" * 100 + "\n", encoding="utf-8")
            refined.write_text(">sample__c000001\n" + "A" * 100 + "\n>sample__c000002|viral_region_11_80\n" + "C" * 70 + "\n", encoding="utf-8")

            header_map = directory / "header.tsv"
            header_columns = ["sample_id", "sequence_id", "parent_sequence_id", "record_type", "original_id", "original_header", "length", "coordinates", "extraction_tool"]
            write_tsv(header_map, header_columns, [
                {"sample_id": "sample", "sequence_id": "sample__c000001", "parent_sequence_id": "", "record_type": "input_contig", "original_id": "original-one", "original_header": "original-one note", "length": "100", "coordinates": "", "extraction_tool": ""},
                {"sample_id": "sample", "sequence_id": "sample__c000002", "parent_sequence_id": "", "record_type": "input_contig", "original_id": "original-two", "original_header": "original-two", "length": "100", "coordinates": "", "extraction_tool": ""},
            ])
            gate = directory / "gate.tsv"
            write_tsv(gate, ["sequence_id", "advance_to_refinement", "discovery_status"], [
                {"sequence_id": "sample__c000001", "advance_to_refinement": "true", "discovery_status": "viral"},
                {"sequence_id": "sample__c000002", "advance_to_refinement": "true", "discovery_status": "likely_viral"},
            ])
            region_map = directory / "regions.tsv"
            region_columns = ["sequence_id", "parent_sequence_id", "record_type", "coordinates", "original_length", "refined_length"]
            write_tsv(region_map, region_columns, [
                {"sequence_id": "sample__c000001", "parent_sequence_id": "", "record_type": "input_contig", "coordinates": "", "original_length": "100", "refined_length": "100"},
                {"sequence_id": "sample__c000002|viral_region_11_80", "parent_sequence_id": "sample__c000002", "record_type": "provirus", "coordinates": "11-80", "original_length": "100", "refined_length": "70"},
            ])

            msl = directory / "MSL41.csv"
            msl.write_text(
                "Realm,Kingdom,Phylum,Class,Order,Family,Genus,Species\n"
                "Adnaviria,Zilligvirae,Taleaviricota,Tokiviricetes,Ligamenvirales,Chiyouviridae,Wargodvirus,Wargodvirus xiongnu\n",
                encoding="utf-8",
            )
            evidence_columns = ["sample_id", "sequence_id", "parent_sequence_id", "tool", "classification", "evidence_strength", "score_type", *run_viharmony.RANK_COLUMNS.values(), "vitap_confidence_level", "vitap_assignment_method", "classification_rank"]
            evidence = directory / "evidence.tsv"
            base_taxonomy = {
                "d__Domain": "d__Viruses", "r__Realm": "r__Adnaviria", "k__Kingdom": "k__Zilligvirae",
                "p__Phylum": "p__Taleaviricota", "c__Class": "c__Tokiviricetes", "o__Order": "o__Ligamenvirales",
                "f__Family": "f__Chiyouviridae", "g__Genus": "g__unclassified", "s__Species": "s__unclassified",
            }
            rows = []
            for sequence_id, tool, strength in [
                ("sample__c000001", "genomad", "strong"),
                ("sample__c000001", "virsorter2", "strong"),
                ("sample__c000002", "genomad", "qualified"),
                ("sample__c000002|viral_region_11_80", "vitap", ""),
            ]:
                row = {column: "" for column in evidence_columns}
                row.update({"sample_id": "sample", "sequence_id": sequence_id, "parent_sequence_id": "", "tool": tool, "classification": "virus", "evidence_strength": strength})
                if tool in {"genomad", "vitap"}:
                    row.update(base_taxonomy)
                    row["classification_rank"] = "family"
                if tool == "vitap":
                    row["vitap_confidence_level"] = "High-confidence"
                    row["vitap_assignment_method"] = "graph"
                rows.append(row)
            write_tsv(evidence, evidence_columns, rows)

            groups = directory / "groups.tsv"
            write_tsv(groups, ["sequence_id", "genus_prediction"], [
                {"sequence_id": "sample__c000002|viral_region_11_80", "genus_prediction": "novel_genus_7_of_Chiyouviridae"}
            ])
            args = Namespace(
                sample_id="sample", input_type="dna", normalized_fasta=normalized,
                header_map=header_map, discovery_gate=gate, refined_fasta=refined,
                region_map=region_map, ictv_msl=msl, evidence=[evidence],
                vcontact3_groups=[groups], audit_mode="full",
                output_prefix=directory / "sample",
            )
            run_viharmony.run(args)

            metadata = read_tsv(directory / "sample.final_metadata.tsv")
            self.assertEqual(metadata[0]["viral_confidence"], "high")
            self.assertEqual(metadata[0]["provirus_coordinates"], "NA")
            self.assertEqual(metadata[1]["provirus_coordinates"], "11-80")
            self.assertEqual(metadata[1]["strict_taxonomy_rank"], "family")
            self.assertIn("g__viharmony_sample_novel_genus_7_of_Chiyouviridae", metadata[1]["analysis_taxonomy"])
            self.assertIn(">original-two|provirus_11_80", (directory / "sample.final.original_ids.fasta").read_text(encoding="utf-8"))
            self.assertTrue((directory / "sample.evidence_audit.tsv.gz").exists())
            manifest = json.loads((directory / "sample.harmonizer_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema_version"], "viharmony-0.1")

    def test_zero_candidates_still_writes_disposition_and_all_final_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            normalized = directory / "sample.normalized.fasta"
            normalized.write_text(">sample__c000001\nAAAA\n", encoding="utf-8")
            refined = directory / "sample.refined.fasta"
            refined.write_text("", encoding="utf-8")
            header = directory / "header.tsv"
            write_tsv(header, ["sample_id", "sequence_id", "original_id", "original_header", "length"], [{"sample_id": "sample", "sequence_id": "sample__c000001", "original_id": "one", "original_header": "one", "length": "4"}])
            gate = directory / "gate.tsv"
            write_tsv(gate, ["sequence_id", "advance_to_refinement", "discovery_status"], [{"sequence_id": "sample__c000001", "advance_to_refinement": "false", "discovery_status": "unresolved"}])
            regions = directory / "regions.tsv"
            write_tsv(regions, ["sequence_id", "parent_sequence_id", "record_type", "coordinates", "original_length", "refined_length"], [])
            msl = directory / "MSL41.csv"
            msl.write_text("Realm,Kingdom,Phylum,Class,Order,Family,Genus,Species\nAdnaviria,Zilligvirae,Taleaviricota,Tokiviricetes,Ligamenvirales,Chiyouviridae,Wargodvirus,Wargodvirus xiongnu\n", encoding="utf-8")
            args = Namespace(sample_id="sample", input_type="dna", normalized_fasta=normalized, header_map=header, discovery_gate=gate, refined_fasta=refined, region_map=regions, ictv_msl=msl, evidence=[], vcontact3_groups=[], audit_mode="compact", output_prefix=directory / "sample")
            run_viharmony.run(args)
            self.assertEqual(read_tsv(directory / "sample.final_metadata.tsv"), [])
            disposition = read_tsv(directory / "sample.sequence_disposition.tsv")
            self.assertEqual(disposition[0]["disposition"], "discovery_noncandidate_insufficient_support")
            self.assertFalse((directory / "sample.evidence_audit.tsv.gz").exists())


if __name__ == "__main__":
    unittest.main()
