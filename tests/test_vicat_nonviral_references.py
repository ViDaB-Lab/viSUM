import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from bin import classify_vicat_nonviral_references as classifier
from bin import download_vicat_refseq_metadata as downloader


ROOT = Path(__file__).resolve().parents[1]


FEATURE_HEADER = [
    "# feature", "class", "assembly", "assembly_unit", "seq_type",
    "chromosome", "genomic_accession", "start", "end", "strand",
    "product_accession", "non-redundant_refseq", "related_accession", "name",
]


def write_feature_table(path: Path, rows: list[list[str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FEATURE_HEADER)
        writer.writerows(rows)


def feature_row(protein: str, sequence_type: str, replicon: str) -> list[str]:
    return [
        "CDS", "with_protein", "GCF_1", "Primary Assembly", sequence_type,
        "", replicon, "1", "300", "+", protein, protein, "", f"name-{protein}",
    ]


class NonviralReferenceTests(unittest.TestCase):
    def test_classifier_preserves_unplaced_but_only_chromosome_is_flank_eligible(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "assembly_accession\tcellular_group\torganism_name\ttaxid\n"
                "GCF_1\tbacteria\tTest organism\t123\n",
                encoding="utf-8",
            )
            package_data = root / "package" / "ncbi_dataset" / "data" / "GCF_1"
            package_data.mkdir(parents=True)
            proteins = {
                "P_CHROM": "MCHROM",
                "P_UNPLACED": "MUNPLACED",
                "P_PLASMID": "MPLASMID",
                "P_PLASTID": "MPLASTID",
                "P_MITO": "MMITO",
                "P_SHARED": "MSHARED",
                "P_UNKNOWN": "MUNKNOWN",
                "P_MISSING": "MMISSING",
            }
            (package_data / "protein.faa").write_text(
                "".join(f">{protein} description\n{sequence}\n" for protein, sequence in proteins.items()),
                encoding="utf-8",
            )
            metadata = root / "metadata" / "bacteria" / "GCF_1"
            metadata.mkdir(parents=True)
            write_feature_table(
                metadata / "GCF_1.feature_table.txt.gz",
                [
                    feature_row("P_CHROM", "chromosome", "NC_1"),
                    feature_row("P_UNPLACED", "unplaced scaffold", "NW_1"),
                    feature_row("P_PLASMID", "plasmid", "NZ_1"),
                    feature_row("P_PLASTID", "chloroplast", "NC_2"),
                    feature_row("P_MITO", "mitochondrion", "NC_3"),
                    feature_row("P_SHARED", "chromosome", "NC_1"),
                    feature_row("P_SHARED", "plasmid", "NZ_1"),
                    feature_row("P_UNKNOWN", "contig", "NZ_2"),
                ],
            )
            output = root / "output"
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "bin" / "classify_vicat_nonviral_references.py"),
                    "--manifest", str(manifest),
                    "--package-root", str(root / "package"),
                    "--metadata-root", str(root / "metadata"),
                    "--output-dir", str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with gzip.open(
                output / "vicat_nonviral_reference_metadata.tsv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                rows = {row["source_protein_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

            self.assertEqual(set(rows), set(proteins) - {"P_UNKNOWN"})
            self.assertEqual(rows["P_CHROM"]["reference_class"], "CELLULAR_CHROMOSOME")
            self.assertEqual(rows["P_CHROM"]["provirus_flank_eligible"], "true")
            self.assertEqual(rows["P_UNPLACED"]["reference_class"], "CELLULAR_UNPLACED")
            self.assertEqual(rows["P_UNPLACED"]["provirus_flank_eligible"], "false")
            self.assertEqual(rows["P_PLASMID"]["reference_class"], "PLASMID")
            self.assertEqual(rows["P_PLASTID"]["reference_class"], "PLASTID")
            self.assertEqual(rows["P_MITO"]["reference_class"], "MITOCHONDRIAL")
            self.assertEqual(rows["P_SHARED"]["reference_class"], "SHARED_NONVIRAL")
            self.assertEqual(
                rows["P_SHARED"]["source_classes"],
                "CELLULAR_CHROMOSOME,PLASMID",
            )
            self.assertEqual(rows["P_SHARED"]["provirus_flank_eligible"], "false")
            self.assertEqual(rows["P_MISSING"]["reference_class"], "CELLULAR_UNPLACED")
            self.assertEqual(rows["P_MISSING"]["provirus_flank_eligible"], "false")
            self.assertEqual(
                rows["P_MISSING"]["classification_note"],
                "catalog_fallback_no_feature_mapping",
            )

            with gzip.open(
                output / "vicat_nonviral_exclusions.tsv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                excluded = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(excluded), 1)
            self.assertEqual(excluded[0]["source_protein_id"], "P_UNKNOWN")
            self.assertEqual(excluded[0]["reason"], "unsupported_sequence_type:contig")

            expected_fasta_counts = {
                "cellular_chromosome": 1,
                "cellular_unplaced": 2,
                "plasmid": 1,
                "plastid": 1,
                "mitochondrial": 1,
                "shared_nonviral": 1,
            }
            for label, expected_count in expected_fasta_counts.items():
                with gzip.open(output / f"{label}.faa.gz", "rt", encoding="utf-8") as handle:
                    self.assertEqual(handle.read().count(">"), expected_count)

    def test_linkage_and_unlocalized_chromosome_are_flank_eligible(self):
        for sequence_type in ("linkage group", "unlocalized scaffold on chromosome"):
            with self.subTest(sequence_type=sequence_type):
                self.assertEqual(
                    classifier.reference_class(sequence_type),
                    "CELLULAR_CHROMOSOME",
                )

    def test_organelle_and_plasmid_sequence_type_variants(self):
        expected = {
            "apicoplast": "PLASTID",
            "unlocalized scaffold on plasmid": "PLASMID",
            "unlocalized scaffold on mitochondrion": "MITOCHONDRIAL",
        }
        for sequence_type, reference_class in expected.items():
            with self.subTest(sequence_type=sequence_type):
                self.assertEqual(
                    classifier.reference_class(sequence_type),
                    reference_class,
                )

    def test_cluster_metadata_preserves_classes_and_cluster_sizes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata.tsv.gz"
            membership = root / "membership.tsv.gz"
            output_tsv = root / "representatives.tsv.gz"
            output_parquet = root / "representatives.parquet"
            output_summary = root / "summary.tsv"

            metadata_columns = [
                "reference_id", "source_protein_id", "reference_class",
                "source_classes", "cellular_group", "source_accession",
                "replicon_accessions", "replicon_types",
                "provirus_flank_eligible", "classification_note",
                "organism_name", "taxid",
            ]
            metadata_rows = [
                [
                    "NONVIRAL|CELLULAR_CHROMOSOME|GCF_1|P1", "P1",
                    "CELLULAR_CHROMOSOME", "CELLULAR_CHROMOSOME", "bacteria",
                    "GCF_1", "NC_1", "chromosome", "true", "", "Organism 1", "1",
                ],
                [
                    "NONVIRAL|CELLULAR_CHROMOSOME|GCF_2|P2", "P2",
                    "CELLULAR_CHROMOSOME", "CELLULAR_CHROMOSOME", "bacteria",
                    "GCF_2", "NC_2", "chromosome", "true", "", "Organism 2", "2",
                ],
                [
                    "NONVIRAL|PLASMID|GCF_1|P3", "P3", "PLASMID", "PLASMID",
                    "bacteria", "GCF_1", "NZ_1", "plasmid", "false", "",
                    "Organism 1", "1",
                ],
            ]
            with gzip.open(metadata, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(metadata_columns)
                writer.writerows(metadata_rows)

            chromosome_representative = metadata_rows[0][0]
            plasmid_representative = metadata_rows[2][0]
            with gzip.open(membership, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["reference_class", "representative_id", "member_id"])
                writer.writerow(
                    ["CELLULAR_CHROMOSOME", chromosome_representative, metadata_rows[0][0]]
                )
                writer.writerow(
                    ["CELLULAR_CHROMOSOME", chromosome_representative, metadata_rows[1][0]]
                )
                writer.writerow(["PLASMID", plasmid_representative, metadata_rows[2][0]])

            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "bin" / "prepare_vicat_nonviral_cluster_metadata.py"),
                    "--metadata", str(metadata),
                    "--membership", str(membership),
                    "--output-tsv", str(output_tsv),
                    "--output-parquet", str(output_parquet),
                    "--output-summary", str(output_summary),
                    "--work-database", str(root / "work.duckdb"),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            with gzip.open(output_tsv, "rt", encoding="utf-8", newline="") as handle:
                rows = {
                    row["reference_id"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(set(rows), {chromosome_representative, plasmid_representative})
            self.assertEqual(rows[chromosome_representative]["cluster_member_count"], "2")
            self.assertEqual(
                rows[chromosome_representative]["reference_class"],
                "CELLULAR_CHROMOSOME",
            )
            self.assertEqual(rows[chromosome_representative]["provirus_flank_eligible"], "true")
            self.assertEqual(rows[plasmid_representative]["provirus_flank_eligible"], "false")

            with output_summary.open(encoding="utf-8", newline="") as handle:
                summary = {
                    row["reference_class"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(summary["CELLULAR_CHROMOSOME"]["input_proteins"], "2")
            self.assertEqual(summary["CELLULAR_CHROMOSOME"]["representatives"], "1")
            self.assertEqual(summary["PLASMID"]["representatives"], "1")

    def test_downloader_derives_refseq_urls_and_reuses_valid_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "GCF_1_ASM1"
            source.mkdir()
            write_feature_table(
                source / "GCF_1_ASM1_feature_table.txt.gz",
                [feature_row("P1", "chromosome", "NC_1")],
            )
            (source / "GCF_1_ASM1_assembly_report.txt").write_text("# assembly report\n", encoding="utf-8")
            assembly = downloader.Assembly("GCF_1", "bacteria", source.as_uri())
            output = root / "output"

            first = downloader.download_assembly(assembly, output, retries=1, timeout=10)
            second = downloader.download_assembly(assembly, output, retries=1, timeout=10)

            self.assertEqual(first["feature_status"], "downloaded")
            self.assertEqual(first["report_status"], "downloaded")
            self.assertEqual(second["feature_status"], "cached")
            self.assertEqual(second["report_status"], "cached")


if __name__ == "__main__":
    unittest.main()
