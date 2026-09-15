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
    extra_columns = sorted(
        {column for row in rows for column in row}.difference(CORE_EVIDENCE_COLUMNS)
    )
    columns = [*CORE_EVIDENCE_COLUMNS, *extra_columns]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t"
        )
        writer.writeheader()
        for row in rows:
            complete = {column: "" for column in columns}
            complete.update(row)
            writer.writerow(complete)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8") as handle:
        return dict(refine_proviral_regions.read_fasta(handle))


class RefineProviralRegionsTests(unittest.TestCase):
    def test_advisory_vicat_cannot_restore_ct3_trimming_authority(self) -> None:
        ct3 = refine_proviral_regions.BoundaryCall("parent", "cenotetaker3", "ct3", 10, 90)
        vicat = refine_proviral_regions.BoundaryCall("parent", "vicat", "vicat", 20, 80)
        self.assertIsNone(refine_proviral_regions.select_boundary([ct3, vicat], False, 0.5))
        self.assertEqual(refine_proviral_regions.select_boundary([ct3, vicat], True, 0.5), ct3)

    def test_disabling_ct3_authority_preserves_authoritative_boundaries(self) -> None:
        ct3 = refine_proviral_regions.BoundaryCall("parent", "cenotetaker3", "ct3", 10, 90)
        checkv = refine_proviral_regions.BoundaryCall("parent", "checkv", "checkv", 20, 80)
        genomad = refine_proviral_regions.BoundaryCall("parent", "genomad", "genomad", 15, 85)
        vicat = refine_proviral_regions.BoundaryCall("parent", "vicat", "vicat", 30, 70)
        for enabled in (False, True):
            self.assertEqual(refine_proviral_regions.select_boundary([ct3, checkv], enabled, 0.5), checkv)
            self.assertEqual(refine_proviral_regions.select_boundary([ct3, checkv, genomad], enabled, 0.5), genomad)
            self.assertEqual(refine_proviral_regions.select_boundary([checkv, vicat], enabled, 0.5), checkv)

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

    def test_exact_full_span_call_preserves_input_contig_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            sequence = "A" * 100
            (directory / "candidates.fasta").write_text(
                f">sample__c000001\n{sequence}\n", encoding="utf-8"
            )
            genomad = directory / "genomad.tsv"
            write_evidence(
                genomad,
                [{
                    "sample_id": "sample",
                    "sequence_id": "sample__c000001|provirus_1_100",
                    "parent_sequence_id": "sample__c000001",
                    "record_type": "provirus",
                    "coordinates": "1-100",
                    "tool": "genomad",
                    "classification": "virus",
                }],
            )

            paths = self.run_refiner(directory, [genomad])

            self.assertEqual(
                read_fasta(paths["output"]), {"sample__c000001": sequence}
            )
            mapping = read_tsv(paths["map"])[0]
            self.assertEqual(mapping["record_type"], "input_contig")
            self.assertEqual(mapping["coordinates"], "")
            self.assertEqual(
                mapping["boundary_status"],
                "selected_full_span_preserved_contig",
            )
            self.assertEqual(
                read_tsv(paths["audit"])[0]["boundary_status"],
                "selected_full_span_preserved_contig",
            )
            summary = read_tsv(paths["summary"])[0]
            self.assertEqual(summary["refined_parent_count"], "0")
            self.assertEqual(summary["full_span_call_preserved_count"], "1")

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

    def test_whole_contig_consensus_vetoes_six_observed_ct3_fragment_patterns(self) -> None:
        observed_patterns = {
            "sample__virus132": (159378, [(1, 133207), (136248, 159378)]),
            "sample__virus212": (609674, [(1, 85464), (92550, 215229), (234967, 445756), (553739, 609674)]),
            "sample__virus213": (127011, [(11025, 50008), (53941, 110648), (114386, 127011)]),
            "sample__virus214": (437255, [(1, 32959), (39392, 152333), (156489, 247597), (252980, 279247), (301071, 406169)]),
            "sample__virus215": (487887, [(1, 11878), (22673, 282991), (330845, 382072), (386308, 424928)]),
            "sample__virus216": (617453, [(120834, 245138), (287071, 406485), (411139, 433498), (438859, 497873), (540843, 617453)]),
        }
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            with (directory / "candidates.fasta").open("w", encoding="utf-8") as handle:
                for parent_id, (length, _) in observed_patterns.items():
                    handle.write(f">{parent_id}\n{'A' * length}\n")

            common = {"sample_id": "sample", "classification": "virus"}
            ct3_rows = []
            genomad_rows = []
            virsorter2_rows = []
            vicat_rows = []
            for parent_id, (_, regions) in observed_patterns.items():
                for start, end in regions:
                    ct3_rows.append({
                        **common,
                        "sequence_id": f"{parent_id}|provirus_{start}_{end}",
                        "parent_sequence_id": parent_id,
                        "record_type": "provirus",
                        "coordinates": f"{start}-{end}",
                        "tool": "cenotetaker3",
                    })
                genomad_rows.append({
                    **common, "sequence_id": parent_id, "record_type": "input_contig",
                    "tool": "genomad", "evidence_strength": "strong",
                })
                virsorter2_rows.append({
                    **common, "sequence_id": f"{parent_id}||full",
                    "parent_sequence_id": parent_id, "record_type": "input_contig",
                    "tool": "virsorter2", "evidence_strength": "qualified",
                })
                vicat_rows.append({
                    **common, "sequence_id": parent_id, "record_type": "input_contig",
                    "tool": "vicat", "evidence_strength": "qualified",
                    "origin_pattern": "predominantly_viral",
                    "viral_supported_loci": "2", "cellular_supported_loci": "0",
                })

            evidence_paths = []
            for name, rows in (
                ("ct3", ct3_rows), ("genomad", genomad_rows),
                ("virsorter2", virsorter2_rows), ("vicat", vicat_rows),
            ):
                path = directory / f"{name}.tsv"
                write_evidence(path, rows)
                evidence_paths.append(path)

            paths = self.run_refiner(
                directory, evidence_paths, allow_ct3_only_refinement=True
            )
            fasta = read_fasta(paths["output"])
            self.assertEqual(set(fasta), set(observed_patterns))
            self.assertEqual(
                {identifier: len(sequence) for identifier, sequence in fasta.items()},
                {identifier: values[0] for identifier, values in observed_patterns.items()},
            )
            mapping = read_tsv(paths["map"])
            self.assertTrue(all(
                row["record_type"] == "input_contig"
                and row["boundary_status"]
                == "preserved_whole_contig_consensus_over_multi_ct3_fragments"
                for row in mapping
            ))
            audit = read_tsv(paths["audit"])
            self.assertTrue(all(row["selected"] == "false" for row in audit))
            self.assertTrue(all(
                row["boundary_status"]
                == "not_selected_multi_ct3_whole_contig_consensus"
                for row in audit
            ))
            summary = read_tsv(paths["summary"])[0]
            self.assertEqual(summary["whole_contig_consensus_veto_count"], "6")
            self.assertEqual(summary["refined_parent_count"], "0")

    def test_whole_contig_consensus_does_not_veto_one_ct3_provirus_region(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            (directory / "candidates.fasta").write_text(
                ">sample__c000001\n" + "A" * 200 + "\n", encoding="utf-8"
            )
            common = {"sample_id": "sample", "classification": "virus"}
            rows_by_name = {
                "ct3": [{
                    **common, "sequence_id": "sample__c000001|provirus_51_150",
                    "parent_sequence_id": "sample__c000001", "record_type": "provirus",
                    "coordinates": "51-150", "tool": "cenotetaker3",
                }],
                "genomad": [{
                    **common, "sequence_id": "sample__c000001", "record_type": "input_contig",
                    "tool": "genomad", "evidence_strength": "strong",
                }],
                "virsorter2": [{
                    **common, "sequence_id": "sample__c000001||full",
                    "parent_sequence_id": "sample__c000001", "record_type": "input_contig",
                    "tool": "virsorter2", "evidence_strength": "qualified",
                }],
                "vicat": [{
                    **common, "sequence_id": "sample__c000001", "record_type": "input_contig",
                    "tool": "vicat", "evidence_strength": "qualified",
                    "origin_pattern": "predominantly_viral", "viral_supported_loci": "5",
                    "cellular_supported_loci": "0",
                }],
            }
            evidence_paths = []
            for name, rows in rows_by_name.items():
                path = directory / f"{name}.tsv"
                write_evidence(path, rows)
                evidence_paths.append(path)
            paths = self.run_refiner(
                directory, evidence_paths, allow_ct3_only_refinement=True
            )
            self.assertEqual(
                read_fasta(paths["output"]),
                {"sample__c000001|viral_region_51_150": "A" * 100},
            )
            self.assertEqual(
                read_tsv(paths["summary"])[0]["whole_contig_consensus_veto_count"],
                "0",
            )

    def test_multi_ct3_regions_require_complete_clean_whole_contig_consensus(self) -> None:
        for case, include_virsorter2, cellular_loci in (
            ("missing_virsorter2", False, "0"),
            ("vicat_cellular_conflict", True, "1"),
        ):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temporary_directory:
                directory = Path(temporary_directory)
                (directory / "candidates.fasta").write_text(
                    ">sample__c000001\n" + "A" * 200 + "\n", encoding="utf-8"
                )
                common = {"sample_id": "sample", "classification": "virus"}
                rows_by_name = {
                    "ct3": [
                        {
                            **common,
                            "sequence_id": "sample__c000001|provirus_1_80",
                            "parent_sequence_id": "sample__c000001",
                            "record_type": "provirus", "coordinates": "1-80",
                            "tool": "cenotetaker3",
                        },
                        {
                            **common,
                            "sequence_id": "sample__c000001|provirus_121_200",
                            "parent_sequence_id": "sample__c000001",
                            "record_type": "provirus", "coordinates": "121-200",
                            "tool": "cenotetaker3",
                        },
                    ],
                    "genomad": [{
                        **common, "sequence_id": "sample__c000001",
                        "record_type": "input_contig", "tool": "genomad",
                        "evidence_strength": "strong",
                    }],
                    "vicat": [{
                        **common, "sequence_id": "sample__c000001",
                        "record_type": "input_contig", "tool": "vicat",
                        "evidence_strength": "qualified",
                        "origin_pattern": "predominantly_viral",
                        "viral_supported_loci": "5",
                        "cellular_supported_loci": cellular_loci,
                    }],
                }
                if include_virsorter2:
                    rows_by_name["virsorter2"] = [{
                        **common, "sequence_id": "sample__c000001||full",
                        "parent_sequence_id": "sample__c000001",
                        "record_type": "input_contig", "tool": "virsorter2",
                        "evidence_strength": "qualified",
                    }]
                evidence_paths = []
                for name, rows in rows_by_name.items():
                    path = directory / f"{name}.tsv"
                    write_evidence(path, rows)
                    evidence_paths.append(path)
                paths = self.run_refiner(
                    directory, evidence_paths, allow_ct3_only_refinement=True
                )
                self.assertEqual(
                    set(read_fasta(paths["output"])),
                    {
                        "sample__c000001|viral_region_1_80",
                        "sample__c000001|viral_region_121_200",
                    },
                )
                self.assertEqual(
                    read_tsv(paths["summary"])[0][
                        "whole_contig_consensus_veto_count"
                    ],
                    "0",
                )

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
