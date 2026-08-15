import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import standardize_vitap


def write_table(
    path: Path, columns: list[str], rows: list[dict[str, str]], delimiter: str = "\t"
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


class StandardizeVitapTests(unittest.TestCase):
    def build_inputs(self, directory: Path) -> dict[str, Path]:
        paths = {
            name: directory / name
            for name in (
                "regions.tsv",
                "best.tsv",
                "all.tsv",
                "fallback.tsv",
                "metadata.tsv",
                "vmr.csv",
                "evidence.tsv",
                "audit.tsv",
            )
        }
        map_columns = [
            "sample_id",
            "input_type",
            "sequence_id",
            "parent_sequence_id",
            "record_type",
            "coordinates",
            "original_length",
            "refined_length",
            "boundary_source",
            "supporting_boundary_tools",
            "boundary_status",
        ]
        write_table(
            paths["regions.tsv"],
            map_columns,
            [
                {
                    "sample_id": "sample",
                    "input_type": "rna",
                    "sequence_id": "sample__c000001",
                    "parent_sequence_id": "",
                    "record_type": "input_contig",
                    "coordinates": "",
                    "original_length": "1000",
                    "refined_length": "1000",
                    "boundary_source": "",
                    "supporting_boundary_tools": "",
                    "boundary_status": "unchanged",
                },
                {
                    "sample_id": "sample",
                    "input_type": "rna",
                    "sequence_id": "sample__c000002|provirus_101_700",
                    "parent_sequence_id": "sample__c000002",
                    "record_type": "provirus",
                    "coordinates": "101-700",
                    "original_length": "900",
                    "refined_length": "600",
                    "boundary_source": "checkv",
                    "supporting_boundary_tools": "checkv",
                    "boundary_status": "selected",
                },
                {
                    "sample_id": "sample",
                    "input_type": "rna",
                    "sequence_id": "sample__c000003",
                    "parent_sequence_id": "",
                    "record_type": "input_contig",
                    "coordinates": "",
                    "original_length": "800",
                    "refined_length": "800",
                    "boundary_source": "",
                    "supporting_boundary_tools": "",
                    "boundary_status": "unchanged",
                },
                {
                    "sample_id": "sample",
                    "input_type": "rna",
                    "sequence_id": "sample__c000004",
                    "parent_sequence_id": "",
                    "record_type": "input_contig",
                    "coordinates": "",
                    "original_length": "500",
                    "refined_length": "500",
                    "boundary_source": "",
                    "supporting_boundary_tools": "",
                    "boundary_status": "unchanged",
                },
            ],
        )
        lineages = {
            "graph": "Species A;GenusA;FamilyA;OrderA;ClassA;PhylumA;KingdomA;RealmA",
            "fallback": "-;-;-;OrderA;ClassA;PhylumA;KingdomA;RealmA",
            "viriform": "-;-;Rhodogtaviriformidae;-;-;-;-;-",
            "reference": "Species Ref;GenusRef;FamilyRef;OrderRef;ClassRef;PhylumRef;KingdomRef;RealmRef",
        }
        best_rows = [
            {
                "Genome_ID": "sample__c000001",
                "lineage": lineages["graph"],
                "lineage_score/participation_index": "2.5",
                "Confidence_level": "High-confidence",
            },
            {
                "Genome_ID": "sample__c000002|provirus_101_700",
                "lineage": lineages["fallback"],
                "lineage_score/participation_index": "0.3333333333333333",
                "Confidence_level": "UniRef90-based",
            },
            {
                "Genome_ID": "sample__c000003",
                "lineage": lineages["viriform"],
                "lineage_score/participation_index": "0.7",
                "Confidence_level": "Medium-confidence",
            },
            {
                "Genome_ID": "REFERENCE_ACCESSION",
                "lineage": lineages["reference"],
                "lineage_score/participation_index": "4.0",
                "Confidence_level": "High-confidence",
            },
        ]
        write_table(paths["best.tsv"], list(best_rows[0]), best_rows)
        all_rows = [
            {
                "Genome_ID": "sample__c000001",
                "lineage": lineages["graph"],
                "lineage_score/participation_index": "2.5",
            },
            {
                "Genome_ID": "sample__c000001",
                "lineage": "-;GenusA;FamilyA;OrderA;ClassA;PhylumA;KingdomA;RealmA",
                "lineage_score/participation_index": "2.8",
            },
            {
                "Genome_ID": "sample__c000003",
                "lineage": lineages["viriform"],
                "lineage_score/participation_index": "0.7",
            },
            {
                "Genome_ID": "REFERENCE_ACCESSION",
                "lineage": lineages["reference"],
                "lineage_score/participation_index": "4.0",
            },
        ]
        write_table(paths["all.tsv"], list(all_rows[0]), all_rows)
        fallback_rows = [
            {
                "genome_id": "sample__c000002|provirus_101_700",
                "taxa_name": "OrderA",
                "participation_index": "0.3333333333333333",
                "taxon_level": "order",
            },
            {
                "genome_id": "REFERENCE_ACCESSION",
                "taxa_name": "OrderRef",
                "participation_index": "1.0",
                "taxon_level": "order",
            },
        ]
        write_table(paths["fallback.tsv"], list(fallback_rows[0]), fallback_rows)
        metadata = {
            "sample_id": "sample",
            "input_type": "rna",
            "refined_sequence_count": "4",
            "best_lineage_row_count": "4",
            "all_lineage_row_count": "4",
            "uniref90_fallback_row_count": "2",
            "include_low_confidence": "false",
            "threads": "8",
            "vitap_version": "1.12",
            "database_path": "/db/DB_VMR-MSL41",
            "database_release": "VMR-MSL41",
            "run_status": "completed",
            "input_fasta": "sample.refined_candidates.fasta",
        }
        write_table(paths["metadata.tsv"], list(metadata), [metadata])
        vmr_rows = [
            {
                "Virus GENBANK accession": "A",
                "Realm": "RealmA",
                "Kingdom": "KingdomA",
                "Phylum": "PhylumA",
                "Class": "ClassA",
                "Order": "OrderA",
                "Family": "FamilyA",
                "Genus": "GenusA",
                "Species": "Species A",
            },
            {
                "Virus GENBANK accession": "R",
                "Realm": "RealmRef",
                "Kingdom": "KingdomRef",
                "Phylum": "PhylumRef",
                "Class": "ClassRef",
                "Order": "OrderRef",
                "Family": "FamilyRef",
                "Genus": "GenusRef",
                "Species": "Species Ref",
            },
            {
                "Virus GENBANK accession": "V",
                "Realm": "-",
                "Kingdom": "-",
                "Phylum": "-",
                "Class": "-",
                "Order": "-",
                "Family": "Rhodogtaviriformidae",
                "Genus": "-",
                "Species": "-",
            },
        ]
        write_table(paths["vmr.csv"], list(vmr_rows[0]), vmr_rows, delimiter=",")
        return paths

    def run_standardizer(self, paths: dict[str, Path]) -> None:
        argv = [
            "standardize_vitap.py",
            "--sample-id",
            "sample",
            "--input-type",
            "rna",
            "--region-map",
            str(paths["regions.tsv"]),
            "--best-lineages",
            str(paths["best.tsv"]),
            "--all-lineages",
            str(paths["all.tsv"]),
            "--uniref-fallback",
            str(paths["fallback.tsv"]),
            "--run-metadata",
            str(paths["metadata.tsv"]),
            "--vmr",
            str(paths["vmr.csv"]),
            "--output-evidence",
            str(paths["evidence.tsv"]),
            "--output-audit",
            str(paths["audit.tsv"]),
        ]
        with patch.object(sys, "argv", argv):
            standardize_vitap.main()

    def test_standardizes_best_rows_and_excludes_inserted_references(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            self.run_standardizer(paths)
            evidence = read_tsv(paths["evidence.tsv"])
            audit = read_tsv(paths["audit.tsv"])

            self.assertEqual(len(evidence), 3)
            self.assertNotIn("REFERENCE_ACCESSION", {row["sequence_id"] for row in evidence})
            graph = next(row for row in evidence if row["sequence_id"] == "sample__c000001")
            self.assertEqual(graph["classification"], "virus")
            self.assertEqual(graph["score_type"], "vitap_lineage_score")
            self.assertEqual(graph["classification_rank"], "species")
            self.assertEqual(graph["r__Realm"], "r__RealmA")
            self.assertEqual(graph["s__Species"], "s__Species A")

            fallback = next(row for row in evidence if "provirus" in row["sequence_id"])
            self.assertEqual(fallback["parent_sequence_id"], "sample__c000002")
            self.assertEqual(fallback["record_type"], "provirus")
            self.assertEqual(fallback["coordinates"], "101-700")
            self.assertEqual(fallback["score_type"], "vitap_uniref90_participation_index")
            self.assertEqual(fallback["classification_rank"], "order")
            self.assertEqual(fallback["g__Genus"], "g__unclassified")

            viriform = next(row for row in evidence if row["sequence_id"] == "sample__c000003")
            self.assertEqual(viriform["classification"], "viriform")
            self.assertEqual(viriform["d__Domain"], "d__unclassified")
            self.assertEqual(viriform["f__Family"], "f__Rhodogtaviriformidae")
            self.assertEqual(viriform["classification_rank"], "family")

            self.assertEqual(len(audit), 5)
            no_call = next(row for row in audit if row["sequence_id"] == "sample__c000004")
            self.assertEqual(no_call["decision"], "no_best_assignment")
            reference = next(row for row in audit if row["sequence_id"] == "REFERENCE_ACCESSION")
            self.assertEqual(reference["decision"], "excluded_non_input_reference")
            self.assertEqual(reference["in_refinement_map"], "false")
            self.assertEqual(reference["accepted"], "false")
            graph_audit = next(row for row in audit if row["sequence_id"] == "sample__c000001")
            self.assertEqual(graph_audit["alternative_lineage_count"], "2")

    def test_rejects_taxon_absent_from_claimed_vmr_rank(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            rows = read_tsv(paths["best.tsv"])
            rows[0]["lineage"] = rows[0]["lineage"].replace("FamilyA", "FamilyMissing")
            write_table(paths["best.tsv"], list(rows[0]), rows)
            with self.assertRaisesRegex(SystemExit, "absent from VMR family"):
                self.run_standardizer(paths)

    def test_zero_candidates_writes_header_only_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = self.build_inputs(Path(temporary_directory))
            region_columns = list(read_tsv(paths["regions.tsv"])[0])
            write_table(paths["regions.tsv"], region_columns, [])
            write_table(paths["best.tsv"], sorted(standardize_vitap.BEST_REQUIRED), [])
            write_table(paths["all.tsv"], sorted(standardize_vitap.ALL_REQUIRED), [])
            write_table(paths["fallback.tsv"], sorted(standardize_vitap.FALLBACK_REQUIRED), [])
            metadata_rows = read_tsv(paths["metadata.tsv"])
            metadata_rows[0].update(
                {
                    "refined_sequence_count": "0",
                    "best_lineage_row_count": "0",
                    "all_lineage_row_count": "0",
                    "uniref90_fallback_row_count": "0",
                    "run_status": "skipped_no_refined_candidates",
                }
            )
            write_table(paths["metadata.tsv"], list(metadata_rows[0]), metadata_rows)
            self.run_standardizer(paths)
            self.assertEqual(read_tsv(paths["evidence.tsv"]), [])
            self.assertEqual(read_tsv(paths["audit.tsv"]), [])


if __name__ == "__main__":
    unittest.main()
