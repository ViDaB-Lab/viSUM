import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "bin"))

import standardize_virbot as standardizer


HEADER_COLUMNS = ["sample_id", "sequence_id", "record_type", "length"]
SCORE_COLUMNS = [
    "Contig_acc",
    "RNA-viral_gene_content",
    "Encoded_proteins_num",
    "Likely_taxa",
]
METADATA_COLUMNS = list(standardizer.METADATA_REQUIRED)
ICTV_COLUMNS = list(standardizer.ICTV_RANKS)


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


class StandardizeVirBotTests(unittest.TestCase):
    def run_standardizer(
        self,
        directory: Path,
        header_rows: list[dict[str, str]],
        score_rows: list[dict[str, str]],
        sequences: dict[str, str],
    ) -> list[dict[str, str]]:
        header_map = directory / "header_map.tsv"
        scores = directory / "scores.csv"
        virus_fasta = directory / "viruses.fasta"
        metadata = directory / "metadata.tsv"
        ictv = directory / "ictv.csv"
        output = directory / "evidence.tsv"

        write_table(header_map, HEADER_COLUMNS, header_rows, "\t")
        write_table(scores, SCORE_COLUMNS, score_rows, ",")
        with virus_fasta.open("w", encoding="utf-8") as handle:
            for sequence_id, sequence in sequences.items():
                handle.write(f">{sequence_id}\n{sequence}\n")

        write_table(
            metadata,
            METADATA_COLUMNS,
            [
                {
                    "sample_id": "sample",
                    "input_type": "rna",
                    "input_sequence_count": str(len(header_rows)),
                    "positive_sequence_count": str(len(score_rows)),
                    "sensitive_mode": "false",
                    "taxa_mode": "TOP",
                    "virbot_version": "1.0",
                    "virbot_revision": "tested-revision",
                    "run_status": (
                        "completed" if score_rows else "completed_no_virus_calls"
                    ),
                    "score_file": "sample.virbot_scores.csv",
                    "virus_fasta": "sample.virbot_virus_sequences.fasta",
                }
            ],
            "\t",
        )
        write_table(
            ictv,
            ICTV_COLUMNS,
            [
                {
                    "Realm": "Riboviria",
                    "Kingdom": "Orthornavirae",
                    "Phylum": "Negarnaviricota",
                    "Class": "Bunyaviricetes",
                    "Order": "Hareavirales",
                    "Family": "Nairoviridae",
                    "Genus": "Orthonairovirus",
                    "Species": "Orthonairovirus testense",
                },
                {
                    "Realm": "Riboviria",
                    "Kingdom": "Orthornavirae",
                    "Phylum": "Negarnaviricota",
                    "Class": "Bunyaviricetes",
                    "Order": "Hareavirales",
                    "Family": "Nairoviridae",
                    "Genus": "Norwavirus",
                    "Species": "Norwavirus testense",
                },
                {
                    "Realm": "Riboviria",
                    "Kingdom": "Orthornavirae",
                    "Phylum": "Pisuviricota",
                    "Class": "Pisoniviricetes",
                    "Order": "Picornavirales",
                    "Family": "Picornaviridae",
                    "Genus": "Enterovirus",
                    "Species": "Enterovirus testense",
                },
            ],
            ",",
        )

        argv = [
            "standardize_virbot.py",
            "--sample-id",
            "sample",
            "--input-type",
            "rna",
            "--header-map",
            str(header_map),
            "--score-table",
            str(scores),
            "--virus-fasta",
            str(virus_fasta),
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

    def test_old_virbot_lineage_is_remapped_from_current_family(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig1",
                        "record_type": "input_contig",
                        "length": "6",
                    }
                ],
                [
                    {
                        "Contig_acc": "contig1",
                        "RNA-viral_gene_content": "1.0",
                        "Encoded_proteins_num": "1",
                        "Likely_taxa": (
                            "Viruses; Riboviria; Orthornavirae; "
                            "Negarnaviricota; Ellioviricetes; Bunyavirales; "
                            "Nairoviridae."
                        ),
                    }
                ],
                {"contig1": "ACGTAC"},
            )

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["classification"], "virus")
        self.assertEqual(row["score_type"], "virbot_rna_viral_gene_fraction")
        self.assertEqual(row["c__Class"], "c__Bunyaviricetes")
        self.assertEqual(row["o__Order"], "o__Hareavirales")
        self.assertEqual(row["f__Family"], "f__Nairoviridae")
        self.assertEqual(row["g__Genus"], "g__unclassified")

    def test_realm_only_lineage_leaves_lower_ranks_unclassified(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig1",
                        "record_type": "input_contig",
                        "length": "4",
                    }
                ],
                [
                    {
                        "Contig_acc": "contig1",
                        "RNA-viral_gene_content": "1.0",
                        "Encoded_proteins_num": "1",
                        "Likely_taxa": "Viruses; Riboviria.",
                    }
                ],
                {"contig1": "ACGT"},
            )

        self.assertEqual(rows[0]["d__Domain"], "d__Viruses")
        self.assertEqual(rows[0]["r__Realm"], "r__Riboviria")
        self.assertEqual(rows[0]["k__Kingdom"], "k__unclassified")
        self.assertEqual(rows[0]["s__Species"], "s__unclassified")

    def test_unrecognized_taxon_is_not_passed_through(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig1",
                        "record_type": "input_contig",
                        "length": "4",
                    }
                ],
                [
                    {
                        "Contig_acc": "contig1",
                        "RNA-viral_gene_content": "0.5",
                        "Encoded_proteins_num": "2",
                        "Likely_taxa": "Viruses; Fabricatedvirus.",
                    }
                ],
                {"contig1": "ACGT"},
            )

        self.assertEqual(rows[0]["d__Domain"], "d__Viruses")
        for column in standardizer.TAXONOMY_COLUMNS[1:]:
            self.assertTrue(rows[0][column].endswith("__unclassified"))
        self.assertNotIn("Fabricatedvirus", rows[0].values())

    def test_zero_calls_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "contig1",
                        "record_type": "input_contig",
                        "length": "4",
                    }
                ],
                [],
                {},
            )

        self.assertEqual(rows, [])

    def test_nonviral_lineage_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            with self.assertRaisesRegex(ValueError, "does not begin with Viruses"):
                self.run_standardizer(
                    Path(temp_directory),
                    [
                        {
                            "sample_id": "sample",
                            "sequence_id": "contig1",
                            "record_type": "input_contig",
                            "length": "4",
                        }
                    ],
                    [
                        {
                            "Contig_acc": "contig1",
                            "RNA-viral_gene_content": "1.0",
                            "Encoded_proteins_num": "1",
                            "Likely_taxa": "Eukaryota; Fabricatedvirus.",
                        }
                    ],
                    {"contig1": "ACGT"},
                )


if __name__ == "__main__":
    unittest.main()
