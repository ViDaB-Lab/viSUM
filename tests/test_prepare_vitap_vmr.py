import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "prepare_vitap_vmr.py"


class PrepareVitapVmrTests(unittest.TestCase):
    def run_helper(self, source: Path, destination: Path, metadata: Path, *extra: str):
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--input",
                str(source),
                "--output-csv",
                str(destination),
                "--metadata",
                str(metadata),
                *extra,
            ],
            text=True,
            capture_output=True,
        )

    def test_prepares_nine_columns_and_audits_vitap_accession_syntax(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            source = directory_path / "VMR_MSL41.csv"
            destination = directory_path / "prepared.csv"
            metadata = directory_path / "metadata.tsv"
            with source.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    [
                        "Virus GENBANK accession",
                        "Realm",
                        "Kingdom",
                        "Phylum",
                        "Class",
                        "Order",
                        "Family",
                        "Genus",
                        "Species",
                        "Ignored column",
                    ]
                )
                writer.writerow(
                    [
                        "DNA-S: MF926439; segment: AE006468 (844298.877981)",
                        "R",
                        "K",
                        "P",
                        "C",
                        "O",
                        "F",
                        "G",
                        "S",
                        "ignored",
                    ]
                )
                writer.writerow(["", "", "", "", "", "", "", "", "", ""])

            result = self.run_helper(source, destination, metadata)
            self.assertEqual(result.returncode, 0, result.stderr)
            with destination.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(list(rows[0]), [
                "Virus GENBANK accession", "Realm", "Kingdom", "Phylum",
                "Class", "Order", "Family", "Genus", "Species",
            ])
            self.assertEqual(len(rows), 1)

            with metadata.open(encoding="utf-8") as handle:
                audit = list(csv.DictReader(handle, delimiter="\t"))[0]
            self.assertEqual(audit["release_label"], "VMR-MSL41")
            self.assertEqual(audit["accession_entries"], "2")
            self.assertEqual(audit["labeled_entries"], "2")
            self.assertEqual(audit["coordinate_entries"], "1")
            self.assertEqual(audit["skipped_rows_without_accession"], "1")

    def test_rejects_genuinely_malformed_accession(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            source = directory_path / "custom.csv"
            destination = directory_path / "prepared.csv"
            metadata = directory_path / "metadata.tsv"
            with source.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "Virus GENBANK accession", "Realm", "Kingdom", "Phylum",
                    "Class", "Order", "Family", "Genus", "Species",
                ])
                writer.writeheader()
                writer.writerow({"Virus GENBANK accession": "not a valid accession"})

            result = self.run_helper(
                source, destination, metadata, "--label", "VMR-custom"
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("malformed accession", result.stderr)


if __name__ == "__main__":
    unittest.main()
