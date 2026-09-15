import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "prepare_vicat_metadata.py"


class TestPrepareVicatMetadata(unittest.TestCase):
    def test_compacts_columns_and_validates_row_count(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "metadata.tsv.gz"
            output = root / "compact.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    [
                        "uvig",
                        "unused",
                        "votu",
                        "ictv_taxonomy",
                        "ictv_taxonomy_method",
                        "genome_type",
                    ]
                )
                writer.writerow(["uvig1", "drop", "votu1", "r__A", "geNomad", "dsDNA"])
                writer.writerow(["uvig2", "drop", "votu2", "r__B", "MMseqs", "ssRNA"])

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input",
                    str(source),
                    "--output",
                    str(output),
                    "--expected-rows",
                    "2",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            with gzip.open(output, "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(
                list(rows[0]),
                [
                    "uvig",
                    "votu",
                    "ictv_taxonomy",
                    "ictv_taxonomy_method",
                    "genome_type",
                ],
            )
            self.assertEqual(rows[1]["uvig"], "uvig2")

    def test_rejects_missing_required_column(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "metadata.tsv.gz"
            output = root / "compact.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                handle.write("uvig\tvotu\tictv_taxonomy\nuvig1\tvotu1\tr__A\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input",
                    str(source),
                    "--output",
                    str(output),
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("ictv_taxonomy_method", result.stderr)


if __name__ == "__main__":
    unittest.main()
