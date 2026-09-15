import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "bin"))

import standardize_deep6 as standardizer


HEADER_COLUMNS = ["sample_id", "sequence_id", "record_type", "length"]
SCORE_COLUMNS = ["name", "length", *standardizer.DEEP6_CLASSES]
METADATA_COLUMNS = [
    "sample_id",
    "input_type",
    "minimum_length",
    "prediction_count",
    "deep6_version",
    "deep6_revision",
    "raw_score_file",
]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def class_scores(winner: str) -> dict[str, str]:
    scores = {deep6_class: "0.02" for deep6_class in standardizer.DEEP6_CLASSES}
    scores[winner] = "0.90"
    return scores


class StandardizeDeep6Tests(unittest.TestCase):
    def run_standardizer(
        self,
        directory: Path,
        header_rows: list[dict[str, str]],
        score_rows: list[dict[str, str]],
        minimum_score: str = "0.7",
        median_multiplier: str = "1.25",
        prediction_count: int | None = None,
    ) -> tuple[list[dict[str, str]], list[str]]:
        header_map = directory / "header_map.tsv"
        scores = directory / "scores.tsv"
        metadata = directory / "metadata.tsv"
        output = directory / "evidence.tsv"

        write_tsv(header_map, HEADER_COLUMNS, header_rows)
        write_tsv(scores, SCORE_COLUMNS, score_rows)
        write_tsv(
            metadata,
            METADATA_COLUMNS,
            [
                {
                    "sample_id": "sample",
                    "input_type": "rna",
                    "minimum_length": "250",
                    "prediction_count": str(
                        len(score_rows)
                        if prediction_count is None
                        else prediction_count
                    ),
                    "deep6_version": "1",
                    "deep6_revision": "tested-revision",
                    "raw_score_file": "sample.deep6_scores.tsv",
                }
            ],
        )

        argv = [
            "standardize_deep6.py",
            "--sample-id",
            "sample",
            "--input-type",
            "rna",
            "--header-map",
            str(header_map),
            "--score-table",
            str(scores),
            "--run-metadata",
            str(metadata),
            "--minimum-score",
            minimum_score,
            "--median-multiplier",
            median_multiplier,
            "--output",
            str(output),
        ]
        with patch.object(sys, "argv", argv):
            standardizer.main()

        with output.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return list(reader), list(reader.fieldnames or [])

    def test_all_deep6_classes_are_mapped(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            header_rows = []
            score_rows = []
            for index, deep6_class in enumerate(standardizer.DEEP6_CLASSES, start=1):
                sequence_id = f"sample__c{index:06d}"
                header_rows.append(
                    {
                        "sample_id": "sample",
                        "sequence_id": sequence_id,
                        "record_type": "input_contig",
                        "length": "500",
                    }
                )
                score_rows.append(
                    {
                        "name": sequence_id,
                        "length": "500",
                        **class_scores(deep6_class),
                    }
                )

            rows, columns = self.run_standardizer(
                directory,
                header_rows,
                score_rows,
            )

            self.assertEqual(columns, standardizer.OUTPUT_COLUMNS)
            self.assertEqual(len(rows), 6)
            by_class = {row["deep6_class"]: row for row in rows}
            for deep6_class in ("duplo", "mono", "ribo", "vari"):
                self.assertEqual(by_class[deep6_class]["classification"], "virus")
                self.assertEqual(by_class[deep6_class]["d__Domain"], "d__Viruses")
                self.assertEqual(
                    by_class[deep6_class]["r__Realm"],
                    f"r__{standardizer.REALM_MAP[deep6_class]}",
                )

            self.assertEqual(by_class["euk"]["classification"], "cellular")
            self.assertEqual(by_class["euk"]["d__Domain"], "d__Eukaryota")
            self.assertEqual(by_class["pro"]["classification"], "cellular")
            self.assertEqual(by_class["pro"]["d__Domain"], "d__Prokaryota")
            for row in rows:
                self.assertEqual(row["score"], "0.90")
                self.assertEqual(row["score_type"], "deep6_top_class_score")
                self.assertEqual(row["record_type"], "input_contig")
                self.assertEqual(row["evidence_strength"], "qualified")
                self.assertEqual(
                    row["strength_basis"],
                    "deep6_top_score_and_median_thresholds",
                )

            for deep6_class in ("duplo", "mono", "ribo", "vari"):
                self.assertEqual(by_class[deep6_class]["evidence_strength"], "qualified")

    def test_zero_predictions_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows, columns = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001",
                        "record_type": "input_contig",
                        "length": "200",
                    }
                ],
                [],
            )

            self.assertEqual(rows, [])
            self.assertEqual(columns, standardizer.OUTPUT_COLUMNS)

    def test_zero_confident_predictions_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows, columns = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001",
                        "record_type": "input_contig",
                        "length": "500",
                    }
                ],
                [
                    {
                        "name": "sample__c000001",
                        "length": "500",
                        "duplo": "0.40",
                        "euk": "0.12",
                        "mono": "0.12",
                        "pro": "0.12",
                        "ribo": "0.12",
                        "vari": "0.12",
                    }
                ],
            )

            self.assertEqual(rows, [])
            self.assertEqual(columns, standardizer.OUTPUT_COLUMNS)

    def test_median_multiplier_filters_ambiguous_prediction(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows, _ = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001",
                        "record_type": "input_contig",
                        "length": "500",
                    }
                ],
                [
                    {
                        "name": "sample__c000001",
                        "length": "500",
                        "duplo": "0.30",
                        "euk": "0.20",
                        "mono": "0.15",
                        "pro": "0.13",
                        "ribo": "0.12",
                        "vari": "0.10",
                    }
                ],
                minimum_score="0.2",
                median_multiplier="2.5",
            )

            self.assertEqual(rows, [])

    def test_empty_scores_are_rejected_when_eligible_sequences_exist(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            with self.assertRaisesRegex(ValueError, "coverage mismatch"):
                self.run_standardizer(
                    Path(temp_directory),
                    [
                        {
                            "sample_id": "sample",
                            "sequence_id": "sample__c000001",
                            "record_type": "input_contig",
                            "length": "500",
                        }
                    ],
                    [],
                )


if __name__ == "__main__":
    unittest.main()
