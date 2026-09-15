from pathlib import Path
import csv
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "prepare_tesorter_input.py"


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if identifier:
                records[identifier] = "".join(sequence)
            identifier = line[1:].split()[0]
            sequence = []
        else:
            sequence.append(line.strip())
    if identifier:
        records[identifier] = "".join(sequence)
    return records


class PrepareTEsorterInputTests(unittest.TestCase):
    def test_only_oversized_sequences_are_windowed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fasta"
            prepared = root / "prepared.fasta"
            mapping = root / "mapping.tsv"
            source.write_text(
                ">short description\nACGTACGT\n"
                ">long\n" + "A" * 21 + "\n",
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--input", str(source),
                    "--output-fasta", str(prepared),
                    "--output-map", str(mapping),
                    "--max-length", "10",
                    "--overlap", "2",
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
            )

            records = read_fasta(prepared)
            self.assertEqual(records["short"], "ACGTACGT")
            long_ids = [identifier for identifier in records if identifier != "short"]
            self.assertEqual(len(long_ids), 3)
            self.assertTrue(all(len(records[identifier]) <= 10 for identifier in long_ids))

            with mapping.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            short_row = next(row for row in rows if row["tesorter_sequence_id"] == "short")
            self.assertEqual(short_row["sequence_id"], "short")
            self.assertEqual(short_row["was_split"], "false")
            long_rows = [row for row in rows if row["sequence_id"] == "long"]
            self.assertEqual(len(long_rows), 3)
            self.assertTrue(all(row["was_split"] == "true" for row in long_rows))
            self.assertEqual(
                [(row["window_start"], row["window_end"]) for row in long_rows],
                [("1", "10"), ("9", "18"), ("17", "21")],
            )


if __name__ == "__main__":
    unittest.main()
