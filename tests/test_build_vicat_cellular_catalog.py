import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class CellularCatalogTests(unittest.TestCase):
    def test_builds_matching_catalog(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.tsv"
            manifest.write_text("assembly_accession\tcellular_group\nGCF_1\tbacteria\n", encoding="utf-8")
            data = root / "package" / "ncbi_dataset" / "data" / "GCF_1"
            data.mkdir(parents=True)
            (data / "protein.faa").write_text(">WP_1 description\nMPEPTIDE\n>WP_2\nMTEST\n", encoding="utf-8")
            fasta = root / "cellular.faa.gz"
            metadata = root / "cellular.tsv.gz"
            subprocess.run(
                [sys.executable, str(ROOT / "bin" / "build_vicat_cellular_catalog.py"), "--manifest", str(manifest),
                 "--package-root", str(root / "package"), "--output-fasta", str(fasta),
                 "--output-metadata", str(metadata)], check=True, capture_output=True, text=True
            )
            with gzip.open(fasta, "rt", encoding="utf-8") as handle:
                self.assertIn(">CELLULAR|GCF_1|WP_1", handle.read())
            with gzip.open(metadata, "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["cellular_group"], "bacteria")


if __name__ == "__main__":
    unittest.main()
