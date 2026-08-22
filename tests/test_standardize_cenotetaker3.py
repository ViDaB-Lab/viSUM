import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "bin"))

import standardize_cenotetaker3 as standardizer


SUMMARY_COLUMNS = [
    "contig",
    "input_name",
    "organism",
    "virus_seq_length",
    "end_feature",
    "gene_count",
    "virion_hallmark_count",
    "rep_hallmark_count",
    "RDRP_hallmark_count",
    "virion_hallmark_genes",
    "rep_hallmark_genes",
    "RDRP_hallmark_genes",
    "taxonomy_hierarchy",
    "ORF_caller",
    "gcode",
    "avg_read_depth",
]

PRUNE_COLUMNS = [
    "contig",
    "contig_length",
    "chunk_length",
    "chunk_name",
    "chunk_start",
    "chunk_stop",
]

GENE_COLUMNS = [
    "contig",
    "gene_start",
    "gene_stop",
    "gene_name",
    "gene_orient",
    "contig_length",
    "dtr_seq",
    "evidence_acession",
    "evidence_description",
    "Evidence_source",
    "vscore_category",
    "chunk_name",
    "chunk_length",
    "chunk_start",
    "chunk_stop",
]

METADATA_COLUMNS = [
    "sample_id",
    "input_type",
    "cenotetaker3_version",
    "run_status",
    "virus_call_count",
]


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class StandardizeCenoteTaker3Tests(unittest.TestCase):
    def run_standardizer(
        self,
        directory: Path,
        summary_rows: list[dict[str, str]],
        prune_rows: list[dict[str, str]],
        gene_rows: list[dict[str, str]],
        fasta_text: str,
        status: str,
    ) -> list[dict[str, str]]:
        header_map = directory / "header_map.tsv"
        summary = directory / "summary.tsv"
        fasta = directory / "viruses.fna"
        prune = directory / "prune.tsv"
        genes = directory / "genes.tsv"
        metadata = directory / "metadata.tsv"
        output = directory / "evidence.tsv"

        write_tsv(
            header_map,
            ["sample_id", "sequence_id", "record_type", "length"],
            [
                {
                    "sample_id": "sample",
                    "sequence_id": "sample__c000001",
                    "record_type": "input_contig",
                    "length": "10",
                },
                {
                    "sample_id": "sample",
                    "sequence_id": "sample__c000002",
                    "record_type": "input_contig",
                    "length": "12",
                },
            ],
        )
        write_tsv(summary, SUMMARY_COLUMNS, summary_rows)
        write_tsv(prune, PRUNE_COLUMNS, prune_rows)
        write_tsv(genes, GENE_COLUMNS, gene_rows)
        fasta.write_text(fasta_text, encoding="utf-8")
        write_tsv(
            metadata,
            METADATA_COLUMNS,
            [
                {
                    "sample_id": "sample",
                    "input_type": "dna",
                    "cenotetaker3_version": "3.4.4",
                    "run_status": status,
                    "virus_call_count": str(len(summary_rows)),
                }
            ],
        )

        argv = [
            "standardize_cenotetaker3.py",
            "--sample-id",
            "sample",
            "--input-type",
            "dna",
            "--header-map",
            str(header_map),
            "--virus-summary",
            str(summary),
            "--virus-fasta",
            str(fasta),
            "--prune-summary",
            str(prune),
            "--gene-annotations",
            str(genes),
            "--run-metadata",
            str(metadata),
            "--output",
            str(output),
        ]
        with patch.object(sys, "argv", argv):
            standardizer.main()

        with output.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_full_and_pruned_calls_are_standardized(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            base_summary = {
                "organism": "Generated pseudo-species that must be ignored",
                "end_feature": "None",
                "virion_hallmark_genes": "Capsid protein|Terminase",
                "rep_hallmark_genes": "",
                "RDRP_hallmark_genes": "",
                "ORF_caller": "prodigal-gv",
                "gcode": "11",
                "avg_read_depth": "NaN",
            }
            summary_rows = [
                {
                    **base_summary,
                    "contig": "ct3_run_1",
                    "input_name": "sample__c000001",
                    "virus_seq_length": "10",
                    "gene_count": "2",
                    "virion_hallmark_count": "2",
                    "rep_hallmark_count": "0",
                    "RDRP_hallmark_count": "0",
                    "taxonomy_hierarchy": "-_Viruses;-_Varidnaviria;k_Bamfordvirae;f_Eupolintoviridae",
                },
                {
                    **base_summary,
                    "contig": "ct3_run_2@C8",
                    "input_name": "sample__c000002",
                    "virus_seq_length": "6",
                    "gene_count": "1",
                    "virion_hallmark_count": "0",
                    "rep_hallmark_count": "1",
                    "RDRP_hallmark_count": "0",
                    "taxonomy_hierarchy": "unclassified virus",
                },
            ]
            prune_rows = [
                {
                    "contig": "ct3_run_2",
                    "contig_length": "12",
                    "chunk_length": "6",
                    "chunk_name": "C8",
                    "chunk_start": "7",
                    "chunk_stop": "13",
                }
            ]
            gene_rows = [
                {
                    "contig": "ct3_run_1",
                    "gene_start": "1",
                    "gene_stop": "5",
                    "gene_name": "ct3_run_1_1",
                    "chunk_name": "NaN",
                },
                {
                    "contig": "ct3_run_1",
                    "gene_start": "6",
                    "gene_stop": "10",
                    "gene_name": "ct3_run_1_2",
                    "chunk_name": "NaN",
                },
                {
                    "contig": "ct3_run_2",
                    "gene_start": "1",
                    "gene_stop": "5",
                    "gene_name": "ct3_run_2_1",
                    "chunk_name": "C8",
                },
            ]

            rows = self.run_standardizer(
                directory,
                summary_rows,
                prune_rows,
                gene_rows,
                ">ct3_run_1\nAAAAAAAAAA\n>ct3_run_2@C8\nCCCCC\n",
                "completed_with_virus_calls",
            )

            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["sequence_id"], "sample__c000001")
            self.assertEqual(rows[0]["score"], "")
            self.assertEqual(rows[0]["r__Realm"], "r__Varidnaviria")
            self.assertEqual(rows[0]["f__Family"], "f__Eupolintoviridae")
            self.assertEqual(rows[0]["g__Genus"], "g__unclassified")
            self.assertEqual(rows[0]["s__Species"], "s__unclassified")
            self.assertEqual(rows[0]["evidence_strength"], "strong")
            self.assertEqual(
                rows[0]["strength_basis"],
                "cenotetaker3_two_hallmark_linear_recommendation",
            )

            self.assertEqual(
                rows[1]["sequence_id"], "sample__c000002|provirus_8_12"
            )
            self.assertEqual(rows[1]["parent_sequence_id"], "sample__c000002")
            self.assertEqual(rows[1]["record_type"], "provirus")
            self.assertEqual(rows[1]["coordinates"], "8-12")
            self.assertEqual(rows[1]["length"], "5")
            self.assertEqual(rows[1]["d__Domain"], "d__Viruses")
            self.assertEqual(rows[1]["r__Realm"], "r__unclassified")
            self.assertEqual(rows[1]["evidence_strength"], "qualified")
            self.assertEqual(
                rows[1]["strength_basis"],
                "cenotetaker3_default_hallmark_requirement",
            )

    def test_zero_call_run_writes_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [],
                [],
                [],
                "",
                "completed_no_viruses_detected",
            )
            self.assertEqual(rows, [])

    def test_repeated_functional_gene_names_at_distinct_loci_are_counted(self) -> None:
        gene_rows = [
            {
                "contig": "ct3_run_1",
                "gene_start": "1",
                "gene_stop": "3",
                "gene_name": "tRNA-Met",
                "chunk_name": "NaN",
            },
            {
                "contig": "ct3_run_1",
                "gene_start": "6",
                "gene_stop": "8",
                "gene_name": "tRNA-Met",
                "chunk_name": "NaN",
            },
            # Exact duplicate annotation rows must not inflate gene_count.
            {
                "contig": "ct3_run_1",
                "gene_start": "6",
                "gene_stop": "8",
                "gene_name": "tRNA-Met",
                "chunk_name": "NaN",
            },
        ]

        self.assertEqual(
            standardizer.load_gene_counts(gene_rows, {"ct3_run_1": 10}),
            {"ct3_run_1": 2},
        )

    def test_zero_hallmark_exploratory_call_is_not_formal_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            rows = self.run_standardizer(
                Path(temp_directory),
                [
                    {
                        "contig": "ct3_run_1",
                        "input_name": "sample__c000001",
                        "organism": "unclassified virus",
                        "virus_seq_length": "10",
                        "end_feature": "None",
                        "gene_count": "1",
                        "virion_hallmark_count": "0",
                        "rep_hallmark_count": "0",
                        "RDRP_hallmark_count": "0",
                        "virion_hallmark_genes": "",
                        "rep_hallmark_genes": "",
                        "RDRP_hallmark_genes": "",
                        "taxonomy_hierarchy": "unclassified virus",
                        "ORF_caller": "prodigal-gv",
                        "gcode": "11",
                        "avg_read_depth": "NaN",
                    }
                ],
                [],
                [
                    {
                        "contig": "ct3_run_1",
                        "gene_start": "1",
                        "gene_stop": "10",
                        "gene_name": "ct3_run_1_1",
                        "chunk_name": "NaN",
                    }
                ],
                ">ct3_run_1\nAAAAAAAAAA\n",
                "completed_with_virus_calls",
            )
            self.assertEqual(rows, [])


if __name__ == "__main__":
    unittest.main()
