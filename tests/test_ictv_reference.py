import csv
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
ICTV_REFERENCE = ROOT / "assets" / "ICTV_VMR_MSL41.csv"
RANKS = ("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


class IctvReferenceTests(unittest.TestCase):
    def test_bundled_msl41_is_a_nonredundant_rank_reference(self) -> None:
        with ICTV_REFERENCE.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            self.assertEqual(tuple(reader.fieldnames or ()), RANKS)
            lineages = [tuple(row[rank].strip() for rank in RANKS) for row in reader]

        self.assertGreater(len(lineages), 17_000)
        self.assertEqual(len(lineages), len(set(lineages)))
        self.assertTrue(all(lineage[-1] for lineage in lineages))

        realms = {lineage[0] for lineage in lineages}
        for current_realm in ("Efunaviria", "Floreoviria", "Pleomoviria", "Volvereviria"):
            self.assertIn(current_realm, realms)


if __name__ == "__main__":
    unittest.main()
