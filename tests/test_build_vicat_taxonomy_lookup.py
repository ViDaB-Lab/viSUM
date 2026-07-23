import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import duckdb
except ImportError:  # pragma: no cover
    duckdb = None


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "bin" / "build_vicat_taxonomy_lookup.py"


@unittest.skipUnless(duckdb is not None, "duckdb is not installed")
class TestViCATTaxonomyLookupBuilder(unittest.TestCase):
    def test_votu_balancing_and_rank_conflicts(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            members = root / "members.tsv"
            metadata = root / "metadata.tsv"
            output = root / "output"
            work = root / "work"

            members.write_text(
                "repA|gene1\tuvigA1|gene1\n"
                "repA|gene1\tuvigA2|gene9\n"
                "repA|gene1\tuvigA3|gene3\n"
                "repB|gene1\tuvigB1|gene1\n"
                "repC|gene1\tuvigC1|gene1\n"
                "repD|gene1\tuvigD1|gene1\n"
                "repD|gene1\tuvigD2|gene2\n",
                encoding="utf-8",
            )
            with metadata.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    [
                        "uvig",
                        "votu",
                        "ictv_taxonomy",
                        "ictv_taxonomy_method",
                    ]
                )
                writer.writerows(
                    [
                        ["uvigA1", "votu1", self.lineage("FamA", "GenA"), "geNomad"],
                        ["uvigA2", "votu1", self.lineage("FamA", "GenA"), "geNomad"],
                        ["uvigA3", "votu2", self.lineage("FamA", "GenB"), "MMseqs"],
                        ["uvigB1", "votu3", self.lineage("FamB", "GenC", "SpecC"), "MMseqs"],
                        ["uvigC1", "votu4", "\\N", "\\N"],
                        ["uvigD1", "votu5", self.lineage("FamD", "GenD"), "geNomad"],
                        ["uvigD2", "votu5", self.lineage("FamE", "GenE"), "MMseqs"],
                    ]
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(BUILDER),
                    "--cluster-members",
                    str(members),
                    "--metadata",
                    str(metadata),
                    "--output-dir",
                    str(output),
                    "--work-dir",
                    str(work),
                    "--threads",
                    "2",
                    "--memory-limit",
                    "1GB",
                    "--prefix",
                    "test",
                    "--expected-member-count",
                    "7",
                    "--expected-representative-count",
                    "4",
                ],
                text=True,
                capture_output=True,
                env=os.environ.copy(),
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)

            relation = duckdb.sql(
                f"SELECT * FROM read_parquet('{(output / 'test.vicat_taxonomy_lookup.parquet').as_posix()}')"
            )
            columns = [column[0] for column in relation.description]
            rows = {row[0]: dict(zip(columns, row)) for row in relation.fetchall()}

            rep_a = rows["repA|gene1"]
            self.assertEqual(rep_a["member_votu_count"], 2)
            self.assertEqual(rep_a["f__Family"], "f__FamA")
            self.assertEqual(rep_a["g__Genus"], "g__unclassified")
            self.assertEqual(rep_a["taxonomy_conflict_rank"], "genus")

            rep_b = rows["repB|gene1"]
            self.assertEqual(rep_b["s__Species"], "s__SpecC")
            self.assertEqual(rep_b["classification_rank"], "species")

            rep_c = rows["repC|gene1"]
            self.assertEqual(rep_c["classified_votu_count"], 0)
            self.assertEqual(rep_c["classification_rank"], "domain")

            rep_d = rows["repD|gene1"]
            self.assertEqual(rep_d["o__Order"], "o__OrderA")
            self.assertEqual(rep_d["f__Family"], "f__unclassified")
            self.assertEqual(rep_d["taxonomy_conflict_rank"], "family")

    @staticmethod
    def lineage(family: str, genus: str, species: str = "") -> str:
        return (
            "r__RealmA;k__KingdomA;p__PhylumA;c__ClassA;o__OrderA;"
            f"f__{family};g__{genus};s__{species}"
        )


if __name__ == "__main__":
    unittest.main()
