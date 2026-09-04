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
