from __future__ import annotations

import csv
import gzip
import json
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import run_viharmony


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_refined_vicat_scope_supersedes_parent_without_counting_twice() -> None:
    parent = {
        "tool": "vicat", "sequence_id": "parent", "classification": "virus",
        "evidence_scope": "parent_discovery",
    }
    region = {
        "tool": "vicat", "sequence_id": "child", "parent_sequence_id": "parent",
        "classification": "cellular", "evidence_scope": "refined_region",
    }
    genomad = {"tool": "genomad", "sequence_id": "parent", "classification": "virus"}
    selected, region_rows, suppressed = run_viharmony.select_vicat_scope(
        [parent, region, genomad], "child", "parent"
    )
    assert suppressed
    assert region_rows == [region]
    assert parent not in selected
    assert region in selected
    assert genomad in selected


def test_parent_nonviral_evidence_becomes_context_only_for_refined_region() -> None:
    parent_cellular = {
        "tool": "deepmicroclass2", "sequence_id": "parent",
        "classification": "cellular",
    }
    parent_plasmid = {
        "tool": "genomad", "sequence_id": "parent",
        "classification": "plasmid",
    }
    region_cellular = {
        "tool": "vicat", "sequence_id": "child",
        "parent_sequence_id": "parent", "classification": "cellular",
        "evidence_scope": "refined_region",
    }
    region_viral = {
        "tool": "genomad", "sequence_id": "child",
        "parent_sequence_id": "parent", "classification": "virus",
    }
    decision, context = run_viharmony.partition_refined_region_context(
        [parent_cellular, parent_plasmid, region_cellular, region_viral],
        "child", "provirus",
    )
    assert context == [parent_cellular, parent_plasmid]
    assert decision == [region_cellular, region_viral]
    unchanged, context = run_viharmony.partition_refined_region_context(
        [parent_cellular], "parent", "input_contig"
    )
    assert unchanged == [parent_cellular]
    assert context == []


class ViharmonyTests(unittest.TestCase):
    def test_disputed_ct3_vicat_boundary_requires_regional_support(self) -> None:
        region = {
            "record_type": "provirus",
            "boundary_source": "cenotetaker3",
            "boundary_status": "selected_boundary_conflict",
            "supporting_boundary_tools": "cenotetaker3,vicat",
        }
        qualified_ct3 = [{
            "tool": "cenotetaker3", "classification": "virus",
            "evidence_strength": "qualified", "record_type": "provirus",
        }]
        passes, basis = run_viharmony.adjudicate_disputed_ct3_boundary(
            region,
            qualified_ct3,
            [{"viral_supported_loci": "3", "cellular_supported_loci": "2"}],
        )
        self.assertTrue(passes)
        self.assertEqual(
            basis, "vicat_regional_support_2_loci_fraction_ge_0.60"
        )

        passes, basis = run_viharmony.adjudicate_disputed_ct3_boundary(
            region,
            qualified_ct3,
            [{"viral_supported_loci": "2", "cellular_supported_loci": "2"}],
        )
        self.assertFalse(passes)
        self.assertEqual(basis, "provisional_disputed_ct3_vicat_boundary")

    def test_strong_ct3_rescues_disputed_boundary(self) -> None:
        region = {
            "record_type": "provirus",
            "boundary_source": "cenotetaker3",
            "boundary_status": "selected_boundary_conflict",
            "supporting_boundary_tools": "cenotetaker3,vicat",
        }
        passes, basis = run_viharmony.adjudicate_disputed_ct3_boundary(
            region,
            [{
                "tool": "cenotetaker3", "classification": "virus",
                "evidence_strength": "strong", "record_type": "provirus",
            }],
            [{"viral_supported_loci": "1", "cellular_supported_loci": "9"}],
        )
        self.assertTrue(passes)
        self.assertEqual(basis, "ct3_strong_multi_hallmark_rescue")

    def test_ct3_single_tool_boundary_is_not_changed_by_conflict_safeguard(self) -> None:
        passes, basis = run_viharmony.adjudicate_disputed_ct3_boundary(
            {
                "record_type": "provirus",
                "boundary_source": "cenotetaker3",
                "boundary_status": "selected_single_tool",
                "supporting_boundary_tools": "cenotetaker3",
            },
            [],
            [],
        )
        self.assertTrue(passes)
        self.assertEqual(basis, "not_disputed_ct3_vicat_boundary")

    def test_viral_entity_interpretations_preserve_plasmid_context(self) -> None:
        self.assertEqual(
            run_viharmony.viral_entity_interpretation(
                "retained_provirus", "provirus", ["genomad"], [], False,
                {"cenotetaker3"},
            ),
            "plasmid_associated_provirus",
        )
        self.assertEqual(
            run_viharmony.viral_entity_interpretation(
                "retained_viral", "input_contig", [], ["genomad"], False,
                {"genomad", "virsorter2"},
            ),
            "virus_plasmid_hybrid_candidate",
        )
        self.assertEqual(
            run_viharmony.viral_entity_interpretation(
                "retained_viral", "input_contig", [], [], False,
                {"genomad", "virsorter2"},
            ),
            "viral_contig",
        )

    def test_disputed_boundary_is_routed_to_provisional_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            normalized = directory / "sample.normalized.fasta"
            refined = directory / "sample.refined.fasta"
            normalized.write_text(
                ">sample__c000001\n" + "A" * 100 + "\n", encoding="utf-8"
            )
            final_id = "sample__c000001|viral_region_11_80"
            refined.write_text(
                f">{final_id}\n" + "A" * 70 + "\n", encoding="utf-8"
            )
            header = directory / "header.tsv"
            write_tsv(
                header,
                ["sample_id", "sequence_id", "original_id", "original_header", "length"],
                [{
                    "sample_id": "sample", "sequence_id": "sample__c000001",
                    "original_id": "original-one", "original_header": "original-one",
                    "length": "100",
                }],
            )
            gate = directory / "gate.tsv"
            write_tsv(
                gate,
                ["sequence_id", "advance_to_refinement", "discovery_status"],
                [{
                    "sequence_id": "sample__c000001",
                    "advance_to_refinement": "true", "discovery_status": "likely_viral",
                }],
            )
            regions = directory / "regions.tsv"
            write_tsv(
                regions,
                [
                    "sequence_id", "parent_sequence_id", "record_type", "coordinates",
                    "original_length", "refined_length", "boundary_source",
                    "supporting_boundary_tools", "boundary_status",
                ],
                [{
                    "sequence_id": final_id, "parent_sequence_id": "sample__c000001",
                    "record_type": "provirus", "coordinates": "11-80",
                    "original_length": "100", "refined_length": "70",
                    "boundary_source": "cenotetaker3",
                    "supporting_boundary_tools": "cenotetaker3,vicat",
                    "boundary_status": "selected_boundary_conflict",
                }],
            )
            evidence = directory / "evidence.tsv"
            evidence_columns = [
                "sample_id", "sequence_id", "parent_sequence_id", "record_type",
                "coordinates", "tool", "classification", "evidence_strength",
                "viral_supported_loci", "cellular_supported_loci",
            ]
            write_tsv(evidence, evidence_columns, [
                {
                    "sample_id": "sample", "sequence_id": final_id,
                    "parent_sequence_id": "sample__c000001", "record_type": "provirus",
                    "coordinates": "11-80", "tool": "cenotetaker3",
                    "classification": "virus", "evidence_strength": "qualified",
                    "viral_supported_loci": "", "cellular_supported_loci": "",
                },
                {
                    "sample_id": "sample", "sequence_id": final_id,
                    "parent_sequence_id": "sample__c000001", "record_type": "provirus",
                    "coordinates": "11-80", "tool": "vicat",
                    "classification": "virus", "evidence_strength": "qualified",
                    "viral_supported_loci": "2", "cellular_supported_loci": "2",
                },
            ])
            msl = directory / "MSL41.csv"
            msl.write_text(
                "Realm,Kingdom,Phylum,Class,Order,Family,Genus,Species\n"
                "Adnaviria,Zilligvirae,Taleaviricota,Tokiviricetes,"
                "Ligamenvirales,Chiyouviridae,Wargodvirus,Wargodvirus xiongnu\n",
                encoding="utf-8",
            )
            args = Namespace(
                sample_id="sample", input_type="dna", normalized_fasta=normalized,
                header_map=header, discovery_gate=gate, refined_fasta=refined,
                region_map=regions, ictv_msl=msl, evidence=[evidence],
                vcontact3_groups=[], audit_mode="full",
                vcontact3_min_taxonomy_length=1,
                output_prefix=directory / "sample",
            )
            run_viharmony.run(args)

            metadata = read_tsv(directory / "sample.final_metadata.tsv")
            self.assertEqual(metadata[0]["viral_decision"], "provisional_provirus")
            self.assertEqual(
                metadata[0]["provirus_retention_basis"],
                "provisional_disputed_ct3_vicat_boundary",
            )
            self.assertEqual(
                (directory / "sample.final.normalized.fasta").read_text(encoding="utf-8"),
                "",
            )
            self.assertIn(
                final_id,
                (directory / "sample.provisional.normalized.fasta").read_text(
                    encoding="utf-8"
                ),
            )

    def test_only_exact_full_span_region_evidence_follows_preserved_contig(self) -> None:
        common = {
            "exact_id": "sample__c000001|provirus_1_100",
            "parent_id": "sample__c000001",
            "row_record_type": "provirus",
            "final_id": "sample__c000001",
            "final_parent": "sample__c000001",
            "final_record_type": "input_contig",
            "final_interval": None,
            "final_length": 100,
            "known_final_ids": {"sample__c000001"},
        }
        self.assertTrue(
            run_viharmony.evidence_applies_to_final(
                row_interval=(1, 100), **common
            )
        )
        self.assertFalse(
            run_viharmony.evidence_applies_to_final(
                row_interval=(10, 90), **common
            )
        )

    def test_tesorter_changes_interpretation_without_voting_or_rejection(self) -> None:
        viral = [{"tool": "genomad", "classification": "virus", "evidence_strength": "qualified"}]

        strong_retroelement = viral + [{
            "tool": "tesorter", "classification": "retroelement",
            "evidence_strength": "strong", "tesorter_evidence_category": "retroelement",
            "tesorter_order": "LTR", "tesorter_superfamily": "Gypsy",
            "assignment_method": "direct_hmm",
        }]
        summary = run_viharmony.summarize_tesorter(strong_retroelement, True)
        self.assertEqual(summary["sequence_interpretation"], "viral_retroelement_conflict")
        self.assertEqual(summary["tesorter_status"], "strong_retroelement_conflict")
        self.assertEqual(summary["tesorter_orders"], "LTR")

        likely = run_viharmony.summarize_tesorter(strong_retroelement[1:], False)
        self.assertEqual(likely["sequence_interpretation"], "likely_retroelement")

        weak = run_viharmony.summarize_tesorter([{
            "tool": "tesorter", "classification": "retroelement",
            "evidence_strength": "weak", "tesorter_evidence_category": "retroelement",
            "tesorter_order": "LTR", "tesorter_superfamily": "Gypsy",
            "assignment_method": "second_pass",
        }], True)
        self.assertEqual(
            weak["sequence_interpretation"],
            "viral_candidate_with_weak_retroelement_signal",
        )

        viral_like = run_viharmony.summarize_tesorter([{
            "tool": "tesorter", "classification": "viral_like_mobile_element",
            "evidence_strength": "qualified",
            "tesorter_evidence_category": "viral_like_mobile_element",
            "tesorter_order": "LTR", "tesorter_superfamily": "Retrovirus",
            "assignment_method": "direct_hmm",
        }], True)
        self.assertEqual(viral_like["sequence_interpretation"], "retrovirus_compatible")

        mixed_viral_like = run_viharmony.summarize_tesorter([
            {
                "tool": "tesorter", "classification": "viral_like_mobile_element",
                "evidence_strength": "qualified",
                "tesorter_evidence_category": "viral_like_mobile_element",
                "tesorter_order": "LTR", "tesorter_superfamily": "Retrovirus",
                "assignment_method": "direct_hmm",
            },
            {
                "tool": "tesorter", "classification": "retroelement",
                "evidence_strength": "strong",
                "tesorter_evidence_category": "retroelement",
                "tesorter_order": "LTR", "tesorter_superfamily": "Gypsy",
                "assignment_method": "direct_hmm",
            },
        ], True)
        self.assertEqual(
            mixed_viral_like["sequence_interpretation"], "retrovirus_compatible"
        )
        self.assertEqual(
            mixed_viral_like["tesorter_status"],
            "viral_like_mobile_element_with_mixed_retroelement_evidence",
        )

        ambiguous = run_viharmony.summarize_tesorter([{
            "tool": "tesorter", "classification": "ambiguous_mobile_element",
            "evidence_strength": "strong",
            "tesorter_evidence_category": "ambiguous_mobile_element",
            "tesorter_order": "Maverick", "tesorter_superfamily": "Polinton",
            "assignment_method": "direct_hmm",
        }], True)
        self.assertEqual(ambiguous["sequence_interpretation"], "ambiguous_mobile_element")

        overridden = run_viharmony.summarize_tesorter(
            strong_retroelement, True, True
        )
        self.assertEqual(
            overridden["sequence_interpretation"],
            "viral_with_mobile_element_features",
        )
        self.assertEqual(
            overridden["tesorter_status"],
            "mobile_element_annotation_overridden_by_high_viral_consensus",
        )

    def test_checkv_complete_zero_host_pattern_is_strong_contig_support(self) -> None:
        rows = [{
            "tool": "checkv",
            "classification": "virus",
            "checkv_quality": "High-quality",
            "provirus": "No",
            "viral_genes": "4",
            "host_genes": "0",
        }]
        self.assertTrue(
            run_viharmony.checkv_supports_complete_viral_contig(
                rows, "input_contig"
            )
        )
        self.assertFalse(
            run_viharmony.checkv_supports_complete_viral_contig(rows, "provirus")
        )

    def test_combines_parent_evidence_region_taxonomy_and_original_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            normalized = directory / "sample.normalized.fasta"
            refined = directory / "sample.refined.fasta"
            normalized.write_text(">sample__c000001\n" + "A" * 100 + "\n>sample__c000002\n" + "C" * 100 + "\n", encoding="utf-8")
            refined.write_text(">sample__c000001\n" + "A" * 100 + "\n>sample__c000002|viral_region_11_80\n" + "C" * 70 + "\n", encoding="utf-8")

            header_map = directory / "header.tsv"
            header_columns = ["sample_id", "sequence_id", "parent_sequence_id", "record_type", "original_id", "original_header", "length", "coordinates", "extraction_tool"]
            write_tsv(header_map, header_columns, [
                {"sample_id": "sample", "sequence_id": "sample__c000001", "parent_sequence_id": "", "record_type": "input_contig", "original_id": "original-one", "original_header": "original-one note", "length": "100", "coordinates": "", "extraction_tool": ""},
                {"sample_id": "sample", "sequence_id": "sample__c000002", "parent_sequence_id": "", "record_type": "input_contig", "original_id": "original-two", "original_header": "original-two", "length": "100", "coordinates": "", "extraction_tool": ""},
            ])
            gate = directory / "gate.tsv"
            write_tsv(gate, ["sequence_id", "advance_to_refinement", "discovery_status"], [
                {"sequence_id": "sample__c000001", "advance_to_refinement": "true", "discovery_status": "viral"},
                {"sequence_id": "sample__c000002", "advance_to_refinement": "true", "discovery_status": "likely_viral"},
            ])
            region_map = directory / "regions.tsv"
            region_columns = ["sequence_id", "parent_sequence_id", "record_type", "coordinates", "original_length", "refined_length"]
            write_tsv(region_map, region_columns, [
                {"sequence_id": "sample__c000001", "parent_sequence_id": "", "record_type": "input_contig", "coordinates": "", "original_length": "100", "refined_length": "100"},
                {"sequence_id": "sample__c000002|viral_region_11_80", "parent_sequence_id": "sample__c000002", "record_type": "provirus", "coordinates": "11-80", "original_length": "100", "refined_length": "70"},
            ])

            msl = directory / "MSL41.csv"
            msl.write_text(
                "Realm,Kingdom,Phylum,Class,Order,Family,Genus,Species\n"
                "Adnaviria,Zilligvirae,Taleaviricota,Tokiviricetes,Ligamenvirales,Chiyouviridae,Wargodvirus,Wargodvirus xiongnu\n",
                encoding="utf-8",
            )
            evidence_columns = ["sample_id", "sequence_id", "parent_sequence_id", "tool", "classification", "evidence_strength", "score_type", *run_viharmony.RANK_COLUMNS.values(), "vitap_confidence_level", "vitap_assignment_method", "classification_rank", "tesorter_evidence_category", "tesorter_order", "tesorter_superfamily", "assignment_method"]
            evidence = directory / "evidence.tsv"
            base_taxonomy = {
                "d__Domain": "d__Viruses", "r__Realm": "r__Adnaviria", "k__Kingdom": "k__Zilligvirae",
                "p__Phylum": "p__Taleaviricota", "c__Class": "c__Tokiviricetes", "o__Order": "o__Ligamenvirales",
                "f__Family": "f__Chiyouviridae", "g__Genus": "g__unclassified", "s__Species": "s__unclassified",
            }
            rows = []
            for sequence_id, tool, strength in [
                ("sample__c000001", "genomad", "strong"),
                ("sample__c000001", "virsorter2", "strong"),
                ("sample__c000001", "deep6", ""),
                ("sample__c000001", "deepmicroclass2", "qualified"),
                ("sample__c000002", "genomad", "qualified"),
                ("sample__c000002", "deepmicroclass2", "qualified"),
                ("sample__c000002|viral_region_11_80", "vitap", ""),
            ]:
                row = {column: "" for column in evidence_columns}
                row.update({"sample_id": "sample", "sequence_id": sequence_id, "parent_sequence_id": "", "tool": tool, "classification": "cellular" if tool == "deepmicroclass2" else "virus", "evidence_strength": strength})
                if tool in {"genomad", "vitap"}:
                    row.update(base_taxonomy)
                    row["classification_rank"] = "family"
                if tool == "vitap":
                    row["vitap_confidence_level"] = "High-confidence"
                    row["vitap_assignment_method"] = "graph"
                rows.append(row)
            tesorter_row = {column: "" for column in evidence_columns}
            tesorter_row.update({
                "sample_id": "sample", "sequence_id": "sample__c000001",
                "tool": "tesorter", "classification": "retroelement",
                "evidence_strength": "strong",
                "tesorter_evidence_category": "retroelement",
                "tesorter_order": "LTR", "tesorter_superfamily": "Gypsy",
                "assignment_method": "direct_hmm",
            })
            rows.append(tesorter_row)
            write_tsv(evidence, evidence_columns, rows)

            groups = directory / "groups.tsv"
            write_tsv(groups, ["sequence_id", "realm_prediction", "family_prediction", "genus_prediction"], [
                {
                    "sequence_id": "sample__c000002|viral_region_11_80",
                    "realm_prediction": "Adnaviria",
                    "family_prediction": "Chiyouviridae",
                    "genus_prediction": "novel_genus_7_of_Chiyouviridae",
                }
            ])
            args = Namespace(
                sample_id="sample", input_type="dna", normalized_fasta=normalized,
                header_map=header_map, discovery_gate=gate, refined_fasta=refined,
                region_map=region_map, ictv_msl=msl, evidence=[evidence],
                vcontact3_groups=[groups], audit_mode="full",
                vcontact3_min_taxonomy_length=1,
                output_prefix=directory / "sample",
            )
            run_viharmony.run(args)

            metadata = read_tsv(directory / "sample.final_metadata.tsv")
            self.assertEqual(len(metadata), 2)
            self.assertEqual(metadata[0]["viral_confidence"], "high")
            self.assertEqual(metadata[0]["viral_decision"], "retained_viral")
            self.assertEqual(metadata[0]["cellular_conflict_tools"], "deepmicroclass2")
            self.assertEqual(metadata[0]["generic_cellular_conflict_overridden"], "true")
            self.assertEqual(metadata[0]["sequence_interpretation"], "viral_with_mobile_element_features")
            self.assertEqual(metadata[0]["tesorter_status"], "mobile_element_annotation_overridden_by_high_viral_consensus")
            self.assertEqual(metadata[0]["strict_taxonomy_rank"], "family")
            self.assertEqual(metadata[0]["provirus_coordinates"], "NA")
            self.assertEqual(metadata[1]["provirus_coordinates"], "11-80")
            self.assertEqual(metadata[1]["cellular_conflict_tools"], "NA")
            self.assertEqual(
                metadata[1]["parent_cellular_context_tools"],
                "deepmicroclass2",
            )
            self.assertEqual(metadata[1]["strict_taxonomy_rank"], "family")
            self.assertIn("g__viharmony_sample_novel_genus_7_of_Chiyouviridae", metadata[1]["analysis_taxonomy"])
            primary_fasta = (directory / "sample.final.original_ids.fasta").read_text(
                encoding="utf-8"
            )
            self.assertIn(">original-one", primary_fasta)
            self.assertNotIn(">original-two|provirus_11_80", primary_fasta)
            provisional_fasta = (
                directory / "sample.provisional.original_ids.fasta"
            ).read_text(encoding="utf-8")
            self.assertIn(">original-two|provirus_11_80", provisional_fasta)
            provisional_metadata = read_tsv(
                directory / "sample.provisional_metadata.tsv"
            )
            self.assertEqual(
                provisional_metadata[0]["viral_decision"],
                "provisional_provirus",
            )
            review_fasta = (directory / "sample.review_candidates.fasta").read_text(
                encoding="utf-8"
            )
            self.assertIn(">sample__c000001", review_fasta)
            self.assertTrue((directory / "sample.evidence_audit.tsv.gz").exists())
            with gzip.open(directory / "sample.evidence_audit.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
                audit = list(csv.DictReader(handle, delimiter="\t"))
            deep6 = next(row for row in audit if row["tool"] == "deep6")
            self.assertEqual(deep6["evidence_strength"], "qualified")
            parent_cellular = next(
                row for row in audit
                if row["tool"] == "deepmicroclass2"
                and row["final_sequence_id"] == "sample__c000002|viral_region_11_80"
            )
            self.assertEqual(
                parent_cellular["evidence_application"], "parent_context_only"
            )
            manifest = json.loads((directory / "sample.harmonizer_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema_version"], "viharmony-0.9")
            self.assertIn(
                "parent-only cellular/plasmid calls are retained as context",
                manifest["policy"]["provirus_adjudication"],
            )

    def test_adjudication_separates_primary_review_and_nonviral_calls(self) -> None:
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "viral", "input_contig", {"genomad"}, set(), [], [], "viral_candidate"
            ),
            "retained_viral",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "ambiguous", "input_contig", set(), {"vicat"}, ["checkv"], [],
                "viral_candidate",
            ),
            "likely_nonviral",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "ambiguous", "input_contig", {"genomad"}, set(), ["checkv"], [],
                "viral_candidate",
            ),
            "ambiguous_review",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "viral", "input_contig", {"genomad", "virsorter2"}, set(),
                [], ["genomad"], "viral_candidate",
            ),
            "retained_viral",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "ambiguous", "provirus", set(), {"checkv"}, ["checkv"], [],
                "viral_candidate",
            ),
            "likely_nonviral",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "likely_viral", "provirus", set(), {"checkv"}, [], [],
                "viral_candidate",
            ),
            "provisional_provirus",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "viral", "provirus", {"genomad"}, set(), ["checkv"], [],
                "viral_candidate",
            ),
            "ambiguous_review",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "viral", "provirus", {"genomad"}, set(), [], [],
                "likely_retroelement",
            ),
            "likely_retroelement",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "ambiguous", "input_contig", {"genomad", "virsorter2"}, set(),
                ["deepmicroclass2"], [], "viral_candidate",
            ),
            "retained_viral",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "ambiguous", "input_contig", {"genomad", "virsorter2"}, set(),
                ["deepmicroclass2", "checkv"], [], "viral_candidate",
            ),
            "ambiguous_review",
        )
        self.assertEqual(
            run_viharmony.adjudicate_viral_decision(
                "ambiguous", "input_contig", {"genomad", "virsorter2"}, set(),
                ["deepmicroclass2"], ["genomad"], "viral_candidate",
            ),
            "ambiguous_review",
        )

    def test_vcontact3_groups_require_length_context_and_unique_rank(self) -> None:
        strict = {rank: "" for rank in run_viharmony.RANKS}
        strict.update({"domain": "Viruses", "realm": "Floreoviria", "family": "Polyomaviridae"})

        short = [{
            "realm_prediction": "Floreoviria",
            "genus_prediction": "novel_genus_1_of_Polyomaviridae",
        }]
        analysis, labels, statuses, audit = run_viharmony.apply_vcontact3_groups(
            short, strict, 556, 1000, "sample"
        )
        self.assertFalse(analysis["genus"])
        self.assertIn("below_vcontact3_minimum_length", statuses)
        self.assertEqual(audit[0]["reason"], "below_vcontact3_minimum_length")
        self.assertTrue(labels)

        compatible_analysis, _, statuses, _ = run_viharmony.apply_vcontact3_groups(
            short, strict, 1500, 1000, "sample"
        )
        self.assertEqual(
            compatible_analysis["genus"],
            "viharmony_sample_novel_genus_1_of_Polyomaviridae",
        )
        self.assertIn("selected_compatible_project_group", statuses)

        incompatible = [{
            "realm_prediction": "Monodnaviria",
            "genus_prediction": "novel_genus_1_within_Monodnaviria",
        }]
        incompatible_analysis, _, statuses, _ = run_viharmony.apply_vcontact3_groups(
            incompatible, strict, 1500, 1000, "sample"
        )
        self.assertFalse(incompatible_analysis["genus"])
        self.assertIn("incompatible_or_unknown_parent_context", statuses)

        ambiguous = [
            {"realm_prediction": "Floreoviria", "genus_prediction": "novel_genus_1"},
            {"realm_prediction": "Floreoviria", "genus_prediction": "novel_genus_2"},
        ]
        ambiguous_analysis, _, statuses, _ = run_viharmony.apply_vcontact3_groups(
            ambiguous, strict, 1500, 1000, "sample"
        )
        self.assertFalse(ambiguous_analysis["genus"])
        self.assertIn("ambiguous_vcontact3_groups", statuses)

    def test_short_vcontact3_formal_vote_is_audit_only(self) -> None:
        msl = {("Floreoviria", "Shotokuvirae", "Cossaviricota", "Papovaviricetes", "Sepolyvirales", "Polyomaviridae", "", "")}
        row = {
            "tool": "vcontact3", "classification": "virus", "length": "556",
            "vcontact3_assignment_method": "realm_only",
            "d__Domain": "d__Viruses", "r__Realm": "r__Floreoviria",
        }
        strict, rank, confidence, tools, votes = run_viharmony.decide_taxonomy(
            [row], msl, 1000
        )
        self.assertEqual(rank, "domain")
        self.assertEqual(confidence, "unclassified")
        self.assertEqual(tools, [])
        self.assertEqual(votes[0]["reason"], "below_vcontact3_minimum_length")

    def test_zero_candidates_still_writes_disposition_and_all_final_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            normalized = directory / "sample.normalized.fasta"
            normalized.write_text(">sample__c000001\nAAAA\n", encoding="utf-8")
            refined = directory / "sample.refined.fasta"
            refined.write_text("", encoding="utf-8")
            header = directory / "header.tsv"
            write_tsv(header, ["sample_id", "sequence_id", "original_id", "original_header", "length"], [{"sample_id": "sample", "sequence_id": "sample__c000001", "original_id": "one", "original_header": "one", "length": "4"}])
            gate = directory / "gate.tsv"
            write_tsv(gate, ["sequence_id", "advance_to_refinement", "discovery_status"], [{"sequence_id": "sample__c000001", "advance_to_refinement": "false", "discovery_status": "unresolved"}])
            regions = directory / "regions.tsv"
            write_tsv(regions, ["sequence_id", "parent_sequence_id", "record_type", "coordinates", "original_length", "refined_length"], [])
            msl = directory / "MSL41.csv"
            msl.write_text("Realm,Kingdom,Phylum,Class,Order,Family,Genus,Species\nAdnaviria,Zilligvirae,Taleaviricota,Tokiviricetes,Ligamenvirales,Chiyouviridae,Wargodvirus,Wargodvirus xiongnu\n", encoding="utf-8")
            args = Namespace(sample_id="sample", input_type="dna", normalized_fasta=normalized, header_map=header, discovery_gate=gate, refined_fasta=refined, region_map=regions, ictv_msl=msl, evidence=[], vcontact3_groups=[], audit_mode="compact", output_prefix=directory / "sample")
            run_viharmony.run(args)
            self.assertEqual(read_tsv(directory / "sample.final_metadata.tsv"), [])
            disposition = read_tsv(directory / "sample.sequence_disposition.tsv")
            self.assertEqual(disposition[0]["disposition"], "discovery_noncandidate_insufficient_support")
            self.assertFalse((directory / "sample.evidence_audit.tsv.gz").exists())


if __name__ == "__main__":
    unittest.main()
