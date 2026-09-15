import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class VicatCellularTierTests(unittest.TestCase):
    def test_tiers_are_balanced_deterministic_and_nested(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "proteins.faa"
            metadata = root / "metadata.tsv"
            output = root / "tiers"
            with fasta.open("w", encoding="utf-8") as fasta_handle, metadata.open("w", encoding="utf-8", newline="") as metadata_handle:
                writer = csv.writer(metadata_handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["protein_id", "cellular_group", "source_accession"])
                for group in ("bacteria", "eukaryota", "archaea"):
                    for index in range(5):
                        protein_id = f"{group}_{index}"
                        writer.writerow([protein_id, group, f"source_{index}"])
                        fasta_handle.write(f">{protein_id} description\nMPEPTIDE{index}\n")

            subprocess.run(
                [sys.executable, str(ROOT / "bin" / "prepare_vicat_cellular_tiers.py"),
                 "--cellular-proteins", str(fasta), "--cellular-metadata", str(metadata),
                 "--output-dir", str(output), "--tiers", "small:3,medium:6,large:9"],
                check=True, capture_output=True, text=True,
            )

            memberships = {}
            for name, size in (("small", 3), ("medium", 6), ("large", 9)):
                with gzip.open(output / f"cellular_{name}.metadata.tsv.gz", "rt", encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(size, len(rows))
                self.assertEqual({"bacteria", "eukaryota", "archaea"}, {row["cellular_group"] for row in rows})
                memberships[name] = {row["protein_id"] for row in rows}
            self.assertLessEqual(memberships["small"], memberships["medium"])
            self.assertLessEqual(memberships["medium"], memberships["large"])


if __name__ == "__main__":
    unittest.main()
