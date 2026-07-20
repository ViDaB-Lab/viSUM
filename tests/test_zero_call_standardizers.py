import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "bin"))

import standardize_genomad
import standardize_virsorter2


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ZeroCallStandardizerTests(unittest.TestCase):
    def test_genomad_zero_viruses_preserves_plasmid_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            header_map = directory / "header_map.tsv"
            virus_summary = directory / "virus_summary.tsv"
            plasmid_summary = directory / "plasmid_summary.tsv"
            metadata = directory / "metadata.tsv"
            output = directory / "evidence.tsv"

            write_tsv(
                header_map,
                ["sample_id", "sequence_id", "record_type"],
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001",
                        "record_type": "input_contig",
                    }
                ],
            )
            write_tsv(
                virus_summary,
                sorted(standardize_genomad.VIRUS_REQUIRED),
                [],
            )
            write_tsv(
                plasmid_summary,
                sorted(standardize_genomad.PLASMID_REQUIRED),
                [
                    {
                        "seq_name": "sample__c000001",
                        "length": "10",
                        "topology": "No terminal repeats",
                        "n_genes": "1",
                        "genetic_code": "11",
                        "plasmid_score": "0.9000",
                        "fdr": "NA",
                        "n_hallmarks": "1",
                        "marker_enrichment": "2.0000",
                    }
                ],
            )
            write_tsv(
                metadata,
                [
                    "sample_id",
                    "input_type",
                    "genomad_version",
                    "run_status",
                    "virus_call_count",
                    "plasmid_call_count",
                ],
                [
                    {
                        "sample_id": "sample",
                        "input_type": "dna",
                        "genomad_version": "1.12.0",
                        "run_status": "completed_no_viruses_detected",
                        "virus_call_count": "0",
                        "plasmid_call_count": "1",
                    }
                ],
            )

            argv = [
                "standardize_genomad.py",
                "--sample-id",
                "sample",
                "--input-type",
                "dna",
                "--header-map",
                str(header_map),
                "--virus-summary",
                str(virus_summary),
                "--plasmid-summary",
                str(plasmid_summary),
                "--run-metadata",
                str(metadata),
                "--output",
                str(output),
            ]
            with patch.object(sys, "argv", argv):
                standardize_genomad.main()

            rows = read_tsv(output)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["classification"], "plasmid")
            self.assertEqual(rows[0]["sequence_id"], "sample__c000001")

    def test_virsorter2_zero_viruses_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            header_map = directory / "header_map.tsv"
            score_table = directory / "score.tsv"
            boundary_table = directory / "boundary.tsv"
            metadata = directory / "metadata.tsv"
            output = directory / "evidence.tsv"

            write_tsv(
                header_map,
                ["sample_id", "sequence_id", "length"],
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001",
                        "length": "2000",
                    }
                ],
            )
            write_tsv(
                score_table,
                sorted(
                    standardize_virsorter2.SCORE_REQUIRED
                    | {"dsDNAphage", "ssDNA"}
                ),
                [],
            )
            write_tsv(
                boundary_table,
                sorted(standardize_virsorter2.BOUNDARY_REQUIRED),
                [],
            )
            write_tsv(
                metadata,
                [
                    "sample_id",
                    "input_type",
                    "virsorter2_version",
                    "classifier_groups",
                    "min_length",
                    "min_score",
                    "run_status",
                    "virus_call_count",
                ],
                [
                    {
                        "sample_id": "sample",
                        "input_type": "dna",
                        "virsorter2_version": "2.2.4",
                        "classifier_groups": "dsDNAphage,ssDNA",
                        "min_length": "1500",
                        "min_score": "0.5",
                        "run_status": "completed_no_viruses_detected",
                        "virus_call_count": "0",
                    }
                ],
            )

            argv = [
                "standardize_virsorter2.py",
                "--sample-id",
                "sample",
                "--input-type",
                "dna",
                "--header-map",
                str(header_map),
                "--score-table",
                str(score_table),
                "--boundary-table",
                str(boundary_table),
                "--run-metadata",
                str(metadata),
                "--output",
                str(output),
            ]
            with patch.object(sys, "argv", argv):
                standardize_virsorter2.main()

            self.assertEqual(read_tsv(output), [])


if __name__ == "__main__":
    unittest.main()
