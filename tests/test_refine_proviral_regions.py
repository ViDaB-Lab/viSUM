import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import refine_proviral_regions
from evidence_schema import CORE_EVIDENCE_COLUMNS


def write_evidence(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=CORE_EVIDENCE_COLUMNS, delimiter="\t"
        )
        writer.writeheader()
        for row in rows:
            complete = {column: "" for column in CORE_EVIDENCE_COLUMNS}
            complete.update(row)
            writer.writerow(complete)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8") as handle:
        return dict(refine_proviral_regions.read_fasta(handle))


class RefineProviralRegionsTests(unittest.TestCase):
    def run_refiner(
        self,
        directory: Path,
        evidence_paths: list[Path],
        input_type: str = "dna",
        allow_ct3_only_refinement: bool = False,
        vicat_support_min_overlap_fraction: float = 0.5,
    ) -> dict[str, Path]:
        paths = {
            "fasta": directory / "candidates.fasta",
            "output": directory / "refined.fasta",
            "map": directory / "map.tsv",
            "audit": directory / "audit.tsv",
            "summary": directory / "summary.tsv",
        }
        argv = [
            "refine_proviral_regions.py",
            "--sample-id",
            "sample",
            "--input-type",
            input_type,
            "--candidate-fasta",
            str(paths["fasta"]),
            "--evidence",
            *[str(path) for path in evidence_paths],
            "--output-fasta",
            str(paths["output"]),
            "--output-map",
            str(paths["map"]),
            "--output-audit",
            str(paths["audit"]),
            "--output-summary",
            str(paths["summary"]),
            "--vicat-support-min-overlap-fraction",
            str(vicat_support_min_overlap_fraction),
        ]
        if allow_ct3_only_refinement:
            argv.append("--allow-ct3-only-refinement")
        with patch.object(sys, "argv", argv):
            refine_proviral_regions.main()
        return paths

    def test_priority_conflicts_and_unchanged_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            (directory / "candidates.fasta").write_text(
                ">sample__c000001\n" + "ACGT" * 30 + "\n"
                ">sample__c000002\n" + "T" * 50 + "\n",
                encoding="utf-8",
            )
            genomad = directory / "genomad.tsv"
            checkv = directory / "checkv.tsv"
            ct3 = directory / "ct3.tsv"
            write_evidence(
                genomad,
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001|provirus_21_80",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "provirus",
                        "coordinates": "21-80",
                        "tool": "genomad",
                        "classification": "virus",
                    }
                ],
            )
            write_evidence(
                checkv,
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001|provirus_25_75",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "provirus",
                        "coordinates": "25-75",
                        "tool": "checkv",
                        "classification": "virus",
                    }
                ],
            )
            write_evidence(
                ct3,
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001|provirus_30_70",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "provirus",
                        "coordinates": "30-70",
                        "tool": "cenotetaker3",
                        "classification": "virus",
                    }
                ],
            )
            paths = self.run_refiner(directory, [genomad, checkv, ct3])

            fasta = read_fasta(paths["output"])
            selected_id = "sample__c000001|viral_region_21_80"
            self.assertEqual(len(fasta[selected_id]), 60)
            self.assertEqual(len(fasta["sample__c000002"]), 50)
            self.assertNotIn("sample__c000001", fasta)

            mapping = {row["sequence_id"]: row for row in read_tsv(paths["map"])}
            self.assertEqual(mapping[selected_id]["boundary_source"], "genomad")
            self.assertEqual(
                mapping[selected_id]["supporting_boundary_tools"],
                "cenotetaker3,checkv,genomad",
            )
            self.assertEqual(
                mapping[selected_id]["boundary_status"],
                "selected_boundary_conflict",
            )
            audit = read_tsv(paths["audit"])
            self.assertEqual(len(audit), 3)
            self.assertEqual(sum(row["selected"] == "true" for row in audit), 1)

    def test_nonoverlapping_calls_create_multiple_regions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            (directory / "candidates.fasta").write_text(
                ">sample__c000001\n" + "A" * 100 + "\n",
                encoding="utf-8",
            )
            evidence = directory / "genomad.tsv"
            write_evidence(
                evidence,
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "a",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "provirus",
                        "coordinates": "1-20",
                        "tool": "genomad",
                        "classification": "virus",
                    },
                    {
                        "sample_id": "sample",
                        "sequence_id": "b",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "provirus",
                        "coordinates": "81-100",
                        "tool": "genomad",
                        "classification": "virus",
                    },
                ],
            )
            paths = self.run_refiner(directory, [evidence])
            fasta = read_fasta(paths["output"])
            self.assertEqual(len(fasta), 2)
            self.assertTrue(all(len(sequence) == 20 for sequence in fasta.values()))

    def test_no_boundaries_preserves_candidate_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            original = ">sample__c000001\nACGTACGT\n"
            (directory / "candidates.fasta").write_text(original, encoding="utf-8")
            paths = self.run_refiner(directory, [])
            self.assertEqual(read_fasta(paths["output"]), {"sample__c000001": "ACGTACGT"})
            self.assertEqual(read_tsv(paths["audit"]), [])
            self.assertEqual(
                read_tsv(paths["map"])[0]["boundary_status"],
                "unchanged_no_provirus_call",
            )

    def test_ct3_only_boundary_is_audited_but_does_not_refine_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            original_sequence = "ACGT" * 30
            (directory / "candidates.fasta").write_text(
                f">sample__c000001\n{original_sequence}\n", encoding="utf-8"
            )
            ct3 = directory / "ct3.tsv"
            write_evidence(
                ct3,
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001|provirus_21_80",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "provirus",
                        "coordinates": "21-80",
                        "tool": "cenotetaker3",
                        "classification": "virus",
                    },
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001|provirus_91_110",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "provirus",
                        "coordinates": "91-110",
                        "tool": "cenotetaker3",
                        "classification": "virus",
                    },
                ],
            )
            paths = self.run_refiner(directory, [ct3])

            self.assertEqual(
                read_fasta(paths["output"]), {"sample__c000001": original_sequence}
            )
            mapping = read_tsv(paths["map"])[0]
            self.assertEqual(
                mapping["boundary_status"],
                "unchanged_ct3_only_boundary_not_allowed",
            )
            audit = read_tsv(paths["audit"])
            self.assertEqual(len(audit), 2)
            self.assertTrue(all(row["selected"] == "false" for row in audit))
            self.assertTrue(all(row["selected_tool"] == "" for row in audit))
            self.assertTrue(
                all(
                    row["boundary_status"] == "not_selected_ct3_only_default"
                    for row in audit
                )
            )
            summary = read_tsv(paths["summary"])[0]
            self.assertEqual(summary["refined_parent_count"], "0")
            self.assertEqual(summary["ct3_only_boundary_call_count"], "2")
            self.assertEqual(summary["ct3_only_locus_skipped_count"], "2")
            self.assertEqual(summary["allow_ct3_only_refinement"], "false")

    def test_ct3_only_boundary_can_refine_when_explicitly_enabled(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            (directory / "candidates.fasta").write_text(
                ">sample__c000001\n" + "A" * 100 + "\n", encoding="utf-8"
            )
            ct3 = directory / "ct3.tsv"
            write_evidence(
                ct3,
                [
                    {
                        "sample_id": "sample",
                        "sequence_id": "sample__c000001|provirus_21_80",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "provirus",
                        "coordinates": "21-80",
                        "tool": "cenotetaker3",
                        "classification": "virus",
                    }
                ],
            )
            paths = self.run_refiner(
                directory, [ct3], allow_ct3_only_refinement=True
            )

            fasta = read_fasta(paths["output"])
            self.assertEqual(
                fasta, {"sample__c000001|viral_region_21_80": "A" * 60}
            )
            self.assertEqual(read_tsv(paths["audit"])[0]["selected"], "true")
            summary = read_tsv(paths["summary"])[0]
            self.assertEqual(summary["refined_parent_count"], "1")
            self.assertEqual(summary["ct3_only_boundary_call_count"], "1")
            self.assertEqual(summary["allow_ct3_only_refinement"], "true")

    def test_vicat_advisory_region_does_not_merge_independent_ct3_loci(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            (directory / "candidates.fasta").write_text(
                ">sample__c000001\n" + "A" * 70000 + "\n", encoding="utf-8"
            )
            ct3 = directory / "ct3.tsv"
            vicat = directory / "vicat.tsv"
            common = {
                "sample_id": "sample",
                "parent_sequence_id": "sample__c000001",
                "record_type": "provirus",
                "classification": "virus",
            }
            write_evidence(
                ct3,
                [
                    {
                        **common,
                        "sequence_id": "sample__c000001|provirus_7367_47239",
                        "coordinates": "7367-47239",
                        "tool": "cenotetaker3",
                    },
                    {
                        **common,
                        "sequence_id": "sample__c000001|provirus_51272_63319",
                        "coordinates": "51272-63319",
                        "tool": "cenotetaker3",
                    },
                ],
            )
            write_evidence(
                vicat,
                [
                    {
                        **common,
                        "sequence_id": "sample__c000001|vicat_provirus_45000_53000",
                        "coordinates": "45000-53000",
                        "tool": "vicat",
                    }
                ],
            )

            paths = self.run_refiner(
                directory,
                [ct3, vicat],
                allow_ct3_only_refinement=True,
            )

            fasta = read_fasta(paths["output"])
            self.assertEqual(
                set(fasta),
                {
                    "sample__c000001|viral_region_7367_47239",
                    "sample__c000001|viral_region_51272_63319",
                },
            )
            self.assertEqual(
                len(fasta["sample__c000001|viral_region_7367_47239"]), 39873
            )
            self.assertEqual(
                len(fasta["sample__c000001|viral_region_51272_63319"]), 12048
            )
            mapping = read_tsv(paths["map"])
            self.assertTrue(
                all(row["boundary_source"] == "cenotetaker3" for row in mapping)
            )
            self.assertEqual(
                read_tsv(paths["summary"])[0]["refined_region_count"], "2"
            )

    def test_vicat_boundary_is_advisory_and_supports_genomad(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            original = "A" * 120
            (directory / "candidates.fasta").write_text(
                f">sample__c000001\n{original}\n", encoding="utf-8"
            )
            vicat = directory / "vicat.tsv"
            write_evidence(
                vicat,
                [{
                    "sample_id": "sample",
                    "sequence_id": "sample__c000001|vicat_provirus_21_80",
                    "parent_sequence_id": "sample__c000001",
                    "record_type": "provirus",
                    "coordinates": "21-80",
                    "tool": "vicat",
                    "classification": "virus",
                }],
            )
            advisory_only = self.run_refiner(directory, [vicat])
            self.assertEqual(
                read_fasta(advisory_only["output"]),
                {"sample__c000001": original},
            )
            self.assertEqual(
                read_tsv(advisory_only["map"])[0]["boundary_status"],
                "unchanged_vicat_advisory_only",
            )
            summary = read_tsv(advisory_only["summary"])[0]
            self.assertEqual(summary["vicat_advisory_boundary_call_count"], "1")
            self.assertEqual(summary["vicat_only_locus_skipped_count"], "1")

            genomad = directory / "genomad.tsv"
            write_evidence(
                genomad,
                [{
                    "sample_id": "sample",
                    "sequence_id": "sample__c000001|provirus_25_75",
                    "parent_sequence_id": "sample__c000001",
                    "record_type": "provirus",
                    "coordinates": "25-75",
                    "tool": "genomad",
                    "classification": "virus",
                }],
            )
            supported_directory = directory / "supported"
            supported_directory.mkdir()
            (supported_directory / "candidates.fasta").write_text(
                f">sample__c000001\n{original}\n", encoding="utf-8"
            )
            supported = self.run_refiner(supported_directory, [vicat, genomad])
            mapping = read_tsv(supported["map"])[0]
            self.assertEqual(mapping["boundary_source"], "genomad")
            self.assertEqual(mapping["supporting_boundary_tools"], "genomad,vicat")
            self.assertEqual(mapping["coordinates"], "25-75")

    def test_checkv_requires_ct3_or_thresholded_vicat_support(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            sequence = "A" * 200
            (directory / "candidates.fasta").write_text(
                f">sample__c000001\n{sequence}\n", encoding="utf-8"
            )
            checkv = directory / "checkv.tsv"
            write_evidence(checkv, [{
                "sample_id": "sample",
                "sequence_id": "sample__c000001|provirus_51_150",
                "parent_sequence_id": "sample__c000001",
                "record_type": "provirus",
                "coordinates": "51-150",
                "tool": "checkv",
                "classification": "virus",
            }])

            checkv_only = self.run_refiner(directory, [checkv])
            self.assertEqual(
                read_fasta(checkv_only["output"]), {"sample__c000001": sequence}
            )
            self.assertEqual(
                read_tsv(checkv_only["audit"])[0]["boundary_status"],
                "not_selected_checkv_only_candidate",
            )

            supported_directory = directory / "supported"
            supported_directory.mkdir()
            (supported_directory / "candidates.fasta").write_text(
                f">sample__c000001\n{sequence}\n", encoding="utf-8"
            )
            vicat = directory / "vicat.tsv"
            write_evidence(vicat, [{
                "sample_id": "sample",
                "sequence_id": "sample__c000001|vicat_provirus_101_180",
                "parent_sequence_id": "sample__c000001",
                "record_type": "provirus",
                "coordinates": "101-180",
                "tool": "vicat",
                "classification": "virus",
            }])
            supported = self.run_refiner(supported_directory, [checkv, vicat])
            self.assertEqual(
                read_fasta(supported["output"]),
                {"sample__c000001|viral_region_51_150": "A" * 100},
            )
            self.assertEqual(
                read_tsv(supported["map"])[0]["boundary_status"],
                "selected_checkv_with_vicat_support",
            )

    def test_checkv_vicat_support_below_overlap_threshold_does_not_trim(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            sequence = "A" * 200
            (directory / "candidates.fasta").write_text(
                f">sample__c000001\n{sequence}\n", encoding="utf-8"
            )
            checkv = directory / "checkv.tsv"
            vicat = directory / "vicat.tsv"
            common = {
                "sample_id": "sample",
                "parent_sequence_id": "sample__c000001",
                "record_type": "provirus",
                "classification": "virus",
            }
            write_evidence(checkv, [{
                **common,
                "sequence_id": "sample__c000001|provirus_51_100",
                "coordinates": "51-100",
                "tool": "checkv",
            }])
            write_evidence(vicat, [{
                **common,
                "sequence_id": "sample__c000001|vicat_provirus_100_200",
                "coordinates": "100-200",
                "tool": "vicat",
            }])
            paths = self.run_refiner(directory, [checkv, vicat])
            self.assertEqual(read_fasta(paths["output"]), {"sample__c000001": sequence})
            self.assertTrue(all(
                row["boundary_status"] == "not_selected_vicat_overlap_below_threshold"
                for row in read_tsv(paths["audit"])
            ))


if __name__ == "__main__":
    unittest.main()
