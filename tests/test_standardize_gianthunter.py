import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "bin"))

import standardize_gianthunter as standardizer


HEADER_COLUMNS = ["sample_id", "sequence_id", "record_type", "length"]
PREDICTION_COLUMNS = [
    "Accession",
    "Length",
    "GiantVirus",
    "PotentialLineage",
    "Score",
]
METADATA_COLUMNS = list(standardizer.METADATA_COLUMNS)
ICTV_COLUMNS = list(standardizer.ICTV_RANKS)
ANNOTATION_COLUMNS = ["Genome", "ORF"]


def write_table(
    path: Path,
    columns: list[str],
    rows: list[dict[str, str]],
    delimiter: str,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=delimiter,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


class StandardizeGiantHunterTests(unittest.TestCase):
    def run_standardizer(
        self,
        directory: Path,
        header_rows: list[dict[str, str]],
        prediction_rows: list[dict[str, str]],
        fasta_sequences: dict[str, str],
        annotation_rows: list[dict[str, str]],
        *,
        eligible_count: int,
        call_count: int,
        run_status: str,
        empty_annotations: bool = False,
        prediction_columns: list[str] = PREDICTION_COLUMNS,
    ) -> list[dict[str, str]]:
        header_map = directory / "header_map.tsv"
        prediction = directory / "prediction.tsv"
        virus_fasta = directory / "viruses.fasta"
        annotations = directory / "annotations.tsv"
        metadata = directory / "metadata.tsv"
        ictv = directory / "ictv.csv"
        output = directory / "evidence.tsv"

        write_table(header_map, HEADER_COLUMNS, header_rows, "\t")
        write_table(prediction, prediction_columns, prediction_rows, "\t")
        with virus_fasta.open("w", encoding="utf-8") as handle:
            for sequence_id, sequence in fasta_sequences.items():
                handle.write(f">{sequence_id}\n{sequence}\n")
        if empty_annotations:
            annotations.write_text("", encoding="utf-8")
        else:
            write_table(
                annotations,
                ANNOTATION_COLUMNS,
                annotation_rows,
                "\t",
            )
        write_table(
            metadata,
            METADATA_COLUMNS,
            [
                {
                    "sample_id": "sample",
                    "input_type": "dna",
                    "input_sequence_count": str(len(header_rows)),
                    "eligible_sequence_count": str(eligible_count),
                    "giant_virus_call_count": str(call_count),
                    "gianthunter_version": "1.1.0",
                    "minimum_length": "3000",
                    "run_status": run_status,
                }
            ],
            "\t",
        )
        write_table(
            ictv,
            ICTV_COLUMNS,
            [
                {
                    "Realm": "Varidnaviria",
                    "Kingdom": "Bamfordvirae",
                    "Phylum": "Nucleocytoviricota",
                    "Class": "Pokkesviricetes",
                    "Order": "Chitovirales",
                    "Family": "Poxviridae",
                    "Genus": "Cervidpoxvirus",
                    "Species": "Cervidpoxvirus muledeerpox",
                },
                {
                    "Realm": "Varidnaviria",
                    "Kingdom": "Bamfordvirae",
                    "Phylum": "Nucleocytoviricota",
                    "Class": "Megaviricetes",
                    "Order": "Imitervirales",
                    "Family": "Schizomimiviridae",
                    "Genus": "Biavirus",
                    "Species": "Biavirus raunefjordenense",
                },
            ],
            ",",
        )

        argv = [
            "standardize_gianthunter.py",
            "--sample-id",
            "sample",
            "--input-type",
            "dna",
            "--header-map",
            str(header_map),
            "--prediction-table",
            str(prediction),
            "--virus-fasta",
            str(virus_fasta),
            "--gene-annotations",
            str(annotations),
            "--run-metadata",
            str(metadata),
            "--ictv-csv",
            str(ictv),
            "--output",
            str(output),
        ]
        with patch.object(sys, "argv", argv):
            standardizer.main()

        with output.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_taxonomy_and_model_calls_use_distinct_score_types(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig_taxonomy",
                        "record_type": "input_contig",
                        "length": "5787",
                    },
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig_model",
                        "record_type": "input_contig",
                        "length": "3543",
                    },
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig_negative",
                        "record_type": "input_contig",
                        "length": "4000",
                    },
                ],
                [
                    {
                        "Accession": "contig_taxonomy",
                        "Length": "5787",
                        "GiantVirus": "GiantVirus",
                        "PotentialLineage": (
                            "superkingdom:Viruses;clade:Varidnaviria;"
                            "kingdom:Bamfordvirae;phylum:Nucleocytoviricota;"
                            "class:Pokkesviricetes;order:Chitovirales;"
                            "family:Poxviridae;genus:Cervidpoxvirus;"
                            "species:Mule deerpox virus"
                        ),
                        "Score": "1.00",
                    },
                    {
                        "Accession": "contig_model",
                        "Length": "3543",
                        "GiantVirus": "GiantVirus",
                        "PotentialLineage": "Unclassified",
                        "Score": "0.969",
                    },
                    {
                        "Accession": "contig_negative",
                        "Length": "4000",
                        "GiantVirus": "Non-GiantVirus",
                        "PotentialLineage": "-",
                        "Score": "-",
                    },
                ],
                {
                    "contig_taxonomy": "A" * 5787,
                    "contig_model": "C" * 3543,
                },
                [
                    {"Genome": "contig_taxonomy", "ORF": "orf1"},
                    {"Genome": "contig_taxonomy", "ORF": "orf2"},
                    {"Genome": "contig_model", "ORF": "orf1"},
                ],
                eligible_count=3,
                call_count=2,
                run_status="completed_with_giant_virus_calls",
            )

        self.assertEqual(len(rows), 2)
        by_id = {row["sequence_id"]: row for row in rows}

        taxonomy = by_id["contig_taxonomy"]
        self.assertEqual(
            taxonomy["score_type"], "gianthunter_weighted_lca_support"
        )
        self.assertEqual(taxonomy["g__Genus"], "g__Cervidpoxvirus")
        self.assertEqual(taxonomy["s__Species"], "s__unclassified")
        self.assertEqual(taxonomy["n_genes"], "2")

        model = by_id["contig_model"]
        self.assertEqual(model["score_type"], "gianthunter_model_score")
        self.assertEqual(model["r__Realm"], "r__Varidnaviria")
        self.assertEqual(model["k__Kingdom"], "k__Bamfordvirae")
        self.assertEqual(model["p__Phylum"], "p__Nucleocytoviricota")
        self.assertEqual(model["c__Class"], "c__unclassified")
        self.assertEqual(model["n_genes"], "1")

    def test_completed_zero_calls_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig1",
                        "record_type": "input_contig",
                        "length": "4000",
                    }
                ],
                [
                    {
                        "Accession": "contig1",
                        "Length": "4000",
                        "GiantVirus": "Non-GiantVirus",
                        "PotentialLineage": "-",
                        "Score": "-",
                    }
                ],
                {},
                [],
                eligible_count=1,
                call_count=0,
                run_status="completed_no_giant_virus_calls",
            )

        self.assertEqual(rows, [])

    def test_skipped_no_eligible_sequences_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig1",
                        "record_type": "input_contig",
                        "length": "1000",
                    }
                ],
                [],
                {},
                [],
                eligible_count=0,
                call_count=0,
                run_status="skipped_no_sequences_meeting_minimum_length",
                empty_annotations=True,
            )

        self.assertEqual(rows, [])

    def test_no_reference_hit_table_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig1",
                        "record_type": "input_contig",
                        "length": "4000",
                    }
                ],
                [
                    {
                        "Accession": "contig1",
                        "Length": "4000",
                        "PotentialLineage": "hits not found in taxonomy files",
                        "Score": "-1",
                        "Genus": "-",
                        "GenusCluster": "-",
                    }
                ],
                {},
                [],
                eligible_count=1,
                call_count=0,
                run_status="completed_no_reference_protein_hits",
                prediction_columns=[
                    "Accession",
                    "Length",
                    "PotentialLineage",
                    "Score",
                    "Genus",
                    "GenusCluster",
                ],
            )

        self.assertEqual(rows, [])

    def test_positive_without_matching_fasta_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            with self.assertRaisesRegex(ValueError, "prediction/FASTA identifiers"):
                self.run_standardizer(
                    Path(temp_directory),
                    [
                        {
                            "sample_id": "sample",
                            "sequence_id": "contig1",
                            "record_type": "input_contig",
                            "length": "4000",
                        }
                    ],
                    [
                        {
                            "Accession": "contig1",
                            "Length": "4000",
                            "GiantVirus": "GiantVirus",
                            "PotentialLineage": "Unclassified",
                            "Score": "0.9",
                        }
                    ],
                    {},
                    [{"Genome": "contig1", "ORF": "orf1"}],
                    eligible_count=1,
                    call_count=1,
                    run_status="completed_with_giant_virus_calls",
                )


if __name__ == "__main__":
    unittest.main()
