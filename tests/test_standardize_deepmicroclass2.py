import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "bin"))

import standardize_deepmicroclass2 as standardizer


HEADER_COLUMNS = ["sample_id", "sequence_id", "record_type", "length"]
SCORE_COLUMNS = [
    "contig",
    "label",
    "confidence",
    *standardizer.DEEPMICROCLASS2_CLASSES,
]
METADATA_COLUMNS = [
    "sample_id",
    "input_type",
    "model_mode",
    "minimum_length",
    "input_sequence_count",
    "eligible_sequence_count",
    "prediction_count",
    "class_thresholds",
    "deepmicroclass2_revision",
    "device",
    "run_status",
    "score_file",
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


def confident_scores(winner: str) -> dict[str, str]:
    winner_scores = {
        "arc": 0.70,
        "bac": 0.50,
        "chlor": 0.50,
        "euk": 0.50,
        "eukvir": 0.96,
        "mit": 0.40,
        "pls": 0.85,
        "prokvir": 0.999,
    }
    winning_score = winner_scores[winner]
    other_score = (1.0 - winning_score) / 7.0
    values = {
        class_name: winning_score if class_name == winner else other_score
        for class_name in standardizer.DEEPMICROCLASS2_CLASSES
    }
    return {class_name: f"{value:.6f}" for class_name, value in values.items()}


class StandardizeDeepMicroClass2Tests(unittest.TestCase):
    def run_standardizer(
        self,
        directory: Path,
        header_rows: list[dict[str, str]],
        score_rows: list[dict[str, str]],
        model_mode: str = "8class",
        minimum_length: int = 500,
        eligible_count: int | None = None,
        prediction_count: int | None = None,
    ) -> tuple[list[dict[str, str]], list[str]]:
        header_map = directory / "header_map.tsv"
        scores = directory / "scores.tsv"
        metadata = directory / "metadata.tsv"
        output = directory / "evidence.tsv"

        write_tsv(header_map, HEADER_COLUMNS, header_rows)
        write_tsv(scores, SCORE_COLUMNS, score_rows)
        expected_eligible = sum(
            int(row["length"]) >= minimum_length for row in header_rows
        )
        write_tsv(
            metadata,
            METADATA_COLUMNS,
            [
                {
                    "sample_id": "sample",
                    "input_type": "dna",
                    "model_mode": model_mode,
                    "minimum_length": str(minimum_length),
                    "input_sequence_count": str(len(header_rows)),
                    "eligible_sequence_count": str(
                        expected_eligible
                        if eligible_count is None
                        else eligible_count
                    ),
                    "prediction_count": str(
                        len(score_rows)
                        if prediction_count is None
                        else prediction_count
                    ),
                    "class_thresholds": standardizer.THRESHOLD_METADATA,
                    "deepmicroclass2_revision": "tested-revision",
                    "device": "cpu",
                    "run_status": (
                        "completed"
                        if score_rows
                        else "completed_no_eligible_sequences"
                    ),
                    "score_file": "sample.deepmicroclass2_scores.tsv",
                }
            ],
        )

        argv = [
            "standardize_deepmicroclass2.py",
            "--sample-id",
            "sample",
            "--input-type",
            "dna",
            "--header-map",
            str(header_map),
            "--score-table",
            str(scores),
            "--run-metadata",
            str(metadata),
            "--output",
            str(output),
        ]
        with patch.object(sys, "argv", argv):
            standardizer.main()

        with output.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return list(reader), list(reader.fieldnames or [])

    def test_all_classes_are_mapped_from_unique_top_probabilities(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            header_rows = []
            score_rows = []
            for index, class_name in enumerate(
                standardizer.DEEPMICROCLASS2_CLASSES,
                start=1,
            ):
                sequence_id = f"sample__c{index:06d}"
                probabilities = confident_scores(class_name)
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
                        "contig": sequence_id,
                        "label": class_name,
                        "confidence": probabilities[class_name],
                        **probabilities,
                    }
                )

            rows, columns = self.run_standardizer(
                directory,
                header_rows,
                score_rows,
            )

            self.assertEqual(columns, standardizer.OUTPUT_COLUMNS)
            self.assertEqual(len(rows), 8)
            by_class = {row["deepmicroclass2_class"]: row for row in rows}
            for class_name in ("eukvir", "prokvir"):
                self.assertEqual(by_class[class_name]["classification"], "virus")
                self.assertEqual(by_class[class_name]["d__Domain"], "d__Viruses")
            for class_name in ("arc", "bac", "chlor", "euk", "mit"):
                self.assertEqual(by_class[class_name]["classification"], "cellular")
            self.assertEqual(by_class["arc"]["d__Domain"], "d__Archaea")
            self.assertEqual(by_class["bac"]["d__Domain"], "d__Bacteria")
            self.assertEqual(by_class["euk"]["d__Domain"], "d__Eukaryota")
            self.assertEqual(by_class["chlor"]["d__Domain"], "d__unclassified")
            self.assertEqual(by_class["mit"]["d__Domain"], "d__unclassified")
            self.assertEqual(by_class["pls"]["classification"], "plasmid")
            for row in rows:
                self.assertEqual(
                    row["score_type"],
                    "deepmicroclass2_top_class_probability",
                )

    def test_bacterial_fallback_is_not_cellular_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            probabilities = {
                "arc": "0.0823",
                "bac": "0.0028",
                "chlor": "0.0004",
                "euk": "0.2397",
                "eukvir": "0.4601",
                "mit": "0.0036",
                "pls": "0.0813",
                "prokvir": "0.1298",
            }
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
                        "contig": "sample__c000001",
                        "label": "bac",
                        "confidence": "0.0028",
                        **probabilities,
                    }
                ],
            )

            self.assertEqual(rows, [])

    def test_lower_probability_class_is_never_selected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            probabilities = {
                "arc": "0.0300",
                "bac": "0.0300",
                "chlor": "0.0300",
                "euk": "0.0700",
                "eukvir": "0.5000",
                "mit": "0.3000",
                "pls": "0.0300",
                "prokvir": "0.0100",
            }
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
                        "contig": "sample__c000001",
                        # The upstream priority rule chooses mit because its
                        # threshold is lower, despite eukvir being the top score.
                        "label": "mit",
                        "confidence": "0.3000",
                        **probabilities,
                    }
                ],
            )

            self.assertEqual(rows, [])

    def test_top_passing_class_replaces_lower_upstream_priority_call(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            probabilities = {
                "arc": "0.0020",
                "bac": "0.0020",
                "chlor": "0.0020",
                "euk": "0.6900",
                "eukvir": "0.0020",
                "mit": "0.3000",
                "pls": "0.0010",
                "prokvir": "0.0010",
            }
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
                        "contig": "sample__c000001",
                        # DeepMicroClass2 checks mit before euk and therefore
                        # reports mit even though euk has the higher score.
                        "label": "mit",
                        "confidence": "0.3000",
                        **probabilities,
                    }
                ],
            )

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["deepmicroclass2_class"], "euk")
            self.assertEqual(rows[0]["classification"], "cellular")
            self.assertEqual(rows[0]["score"], "0.6900")

    def test_high_precision_accepts_300_to_499_nt_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            probabilities = confident_scores("eukvir")
            rows, _ = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001",
                        "record_type": "input_contig",
                        "length": "350",
                    }
                ],
                [
                    {
                        "contig": "sample__c000001",
                        "label": "eukvir",
                        "confidence": probabilities["eukvir"],
                        **probabilities,
                    }
                ],
                model_mode="high_precision",
                minimum_length=300,
            )

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["sequence_id"], "sample__c000001")

    def test_no_eligible_sequences_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows, columns = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001",
                        "record_type": "input_contig",
                        "length": "299",
                    }
                ],
                [],
                model_mode="high_precision",
                minimum_length=300,
            )

            self.assertEqual(rows, [])
            self.assertEqual(columns, standardizer.OUTPUT_COLUMNS)

    def test_missing_eligible_prediction_is_rejected(self) -> None:
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
