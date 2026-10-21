"""Reporting ceilings must not become new inference or retention rules."""
from argparse import Namespace
import copy
import csv
from dataclasses import replace
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import run_viharmony as harmony
from reporting_test_paths import artifact
from test_viharmony_reconciliation import fixture, BUNDLE


def table(path):
    with artifact(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compressed_table(path):
    with gzip.open(artifact(path), "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def example():
    lineage = dict(zip(harmony.RANKS[1:], ("R", "K", "P", "C", "O", "F", "G", "G one")))
    catalog = harmony.make_msl_catalog([tuple(lineage.values())])
    rows = [dict(sequence_id="c1", tool="vitap", classification="virus", length="2000",
                 evidence_strength="qualified", vitap_confidence_level="High-confidence",
                 **{harmony.RANK_COLUMNS[r]: value for r, value in lineage.items()}),
            dict(sequence_id="c1", tool="genomad", classification="virus", length="2000",
                 evidence_strength="strong", **{harmony.RANK_COLUMNS[r]: value for r, value in lineage.items()})]
    task = harmony.CandidateTask("c1", 2000,
        dict(record_type="input_contig", coordinates="", original_length="2000"), "c1",
        dict(original_id="original", original_header="original comment"), "viral", rows, [], ())
    context = harmony.CandidateContext("sample", "dna", 1000, False, False, False, True)
    return catalog, task, context


class ReportingCapTests(unittest.TestCase):
    def test_cli_defaults_to_family_and_accepts_every_supported_rank(self):
        required = ["--sample-id", "x", "--input-type", "dna", "--normalized-fasta", "x",
                    "--header-map", "x", "--discovery-gate", "x", "--refined-fasta", "x",
                    "--region-map", "x", "--ictv-msl", "x", "--output-prefix", "x"]
        from unittest.mock import patch
        with patch.object(sys, "argv", ["run_viharmony", *required]):
            self.assertEqual(harmony.parse_args().max_rank, "family")
        for rank in harmony.RANKS:
            with patch.object(sys, "argv", ["run_viharmony", *required, "--max-rank", rank]):
                self.assertEqual(harmony.parse_args().max_rank, rank)

    def test_invalid_rank_rejected_before_any_input_is_read(self):
        for rank in ("strain", "", None, False, 1):
            with self.subTest(rank=rank), self.assertRaisesRegex(ValueError, "maximum reporting rank"):
                harmony.run(Namespace(workers=1, max_rank=rank))

    def test_all_ceiling_ranks_mask_only_deeper_taxa(self):
        catalog, task, context = example()
        uncapped = harmony.evaluate_candidate(task, replace(context, max_rank="species"), catalog)
        original = copy.deepcopy(task)
        for rank in harmony.RANKS:
            with self.subTest(rank=rank):
                result = harmony.evaluate_candidate(task, replace(context, max_rank=rank), catalog)
                self.assertEqual(result.metadata["viral_decision"], uncapped.metadata["viral_decision"])
                for mode in ("strict", "exploratory"):
                    self.assertEqual(result.metadata[mode + "_taxonomy_rank"], rank)
                    taxa = result.metadata[mode + "_taxonomy"].split(";")
                    self.assertEqual(len(taxa), len(harmony.RANKS))
                    for index, name in enumerate(harmony.RANKS):
                        if index > harmony.RANKS.index(rank):
                            self.assertEqual(taxa[index], harmony.PREFIXES[name] + "unclassified")
                for name in harmony.REPORTING_TAXONOMY_COLUMNS:
                    self.assertEqual(result.reporting_audit["uncapped_" + name], uncapped.metadata[name])
                self.assertEqual(result.evidence_audit, uncapped.evidence_audit)
                self.assertEqual(result.taxonomy_audit, uncapped.taxonomy_audit)
        self.assertEqual(task, original)

    def test_family_ceiling_never_promotes_order_only_evidence(self):
        catalog, task, context = example()
        rows = [{key: value for key, value in row.items()
                 if key not in {harmony.RANK_COLUMNS[r] for r in ("family", "genus", "species")}}
                for row in task.evidence]
        result = harmony.evaluate_candidate(replace(task, evidence=rows), context, catalog)
        self.assertEqual(result.metadata["strict_taxonomy_rank"], "order")
        self.assertIn("f__unclassified;g__unclassified;s__unclassified", result.metadata["strict_taxonomy"])
        self.assertEqual(result.reporting_audit["reporting_status"], "unchanged")

    def test_species_only_call_does_not_invent_direct_family_support(self):
        catalog, task, context = example()
        row = dict(sequence_id="c1", tool="vitap", classification="virus", length="2000",
                   evidence_strength="qualified", vitap_confidence_level="High-confidence", s__Species="G one")
        result = harmony.evaluate_candidate(replace(task, evidence=[row]), context, catalog)
        self.assertEqual(result.metadata["strict_taxonomy_rank"], "family")
        self.assertEqual(result.metadata["strict_taxonomy_confidence"], "lineage_inferred")
        self.assertEqual(result.metadata["strict_taxonomy_supporting_tools"], "NA")
        self.assertEqual(result.reporting_audit["uncapped_strict_taxonomy_confidence"], "single_method_high")

    def test_cap_does_not_modify_inference_engine_or_decision_object(self):
        catalog, task, context = example()
        taxonomy = harmony.decide_taxonomy(task.evidence, catalog)
        before = copy.deepcopy(taxonomy)
        metadata = harmony.evaluate_candidate(task, replace(context, max_rank="species"), catalog).metadata
        harmony.apply_taxonomy_reporting_cap(metadata, taxonomy, "family", "sample")
        self.assertEqual(taxonomy, before)
        self.assertEqual(taxonomy.strict_rank, "species")

    def test_conflicts_and_review_reasons_survive_ceiling(self):
        catalog, task, context = example()
        other = ("R", "K", "P", "C", "O", "F", "H", "H two")
        catalog = harmony.make_msl_catalog([catalog.lineages[0], other])
        alt = {harmony.RANK_COLUMNS[r]: value for r, value in zip(harmony.RANKS[1:], other)}
        rows = [task.evidence[0],
                dict(sequence_id="c1", tool="genomad", classification="virus", evidence_strength="strong", **alt),
                dict(sequence_id="c1", tool="vicat", classification="virus", evidence_strength="qualified", **alt)]
        task = replace(task, evidence=rows)
        uncapped = harmony.evaluate_candidate(task, replace(context, max_rank="species"), catalog)
        capped = harmony.evaluate_candidate(task, context, catalog)
        self.assertEqual(capped.metadata["taxonomy_conflict"], "true")
        self.assertEqual(capped.metadata["strict_taxonomy_conflict_rank"], "genus")
        self.assertEqual(capped.metadata["strict_taxonomy_conflicting_taxa"], "NA")
        self.assertEqual(capped.review["review_reasons"], uncapped.review["review_reasons"])
        self.assertEqual(capped.reporting_audit["uncapped_strict_taxonomy_conflicting_taxa"],
                         uncapped.metadata["strict_taxonomy_conflicting_taxa"])

    def test_every_audit_mode_preserves_uncapped_consensus_and_retained_outputs(self):
        for mode in ("none", "compact", "full"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                args = fixture(root, mode)
                args.ictv_reconciliation = BUNDLE
                outputs = {}
                for rank in ("species", "family", "genus"):
                    dest = root / rank
                    dest.mkdir()
                    args.max_rank = rank
                    args.output_prefix = dest / "sample"
                    harmony.run(args)
                    outputs[rank] = dest
                baseline = table(outputs["species"] / "sample.final_metadata.tsv")
                for rank in ("family", "genus"):
                    dest = outputs[rank]
                    audits = compressed_table(dest / "sample.taxonomy_reporting_audit.tsv.gz")
                    self.assertEqual(len(audits), len(baseline))
                    for original, audit, capped in zip(baseline, audits, table(dest / "sample.final_metadata.tsv"), strict=True):
                        for key, value in original.items():
                            if key not in harmony.REPORTING_TAXONOMY_COLUMNS:
                                self.assertEqual(capped[key], value, key)
                        for key in harmony.REPORTING_TAXONOMY_COLUMNS:
                            self.assertEqual(audit["uncapped_" + key], original[key], key)
                    self.assertTrue(any(r["reporting_status"] == "capped" for r in audits))
                    for path in outputs["species"].rglob("*.fasta"):
                        self.assertEqual(path.read_bytes(), (dest / path.relative_to(outputs["species"])).read_bytes(), path.name)
                    for suffix in ("sequence_map.tsv", "sequence_disposition.tsv"):
                        self.assertEqual(artifact(outputs["species"] / ("sample." + suffix)).read_bytes(),
                                         artifact(dest / ("sample." + suffix)).read_bytes())
                    manifest = json.loads(artifact(dest / "sample.harmonizer_manifest.json").read_text())
                    self.assertEqual(manifest["policy"]["harmonizer_max_rank"], rank)
                    self.assertIn("taxonomy_reporting_audit", manifest["outputs"])
                    self.assertTrue(all(harmony.sha256(dest / rec["file"]) == rec["sha256"]
                                        for rec in manifest["outputs"].values()))
                    if mode == "full":
                        for suffix in ("evidence_audit.tsv.gz", "taxonomy_vote_audit.tsv.gz"):
                            self.assertEqual(gzip.decompress(artifact(outputs["species"] / ("sample." + suffix)).read_bytes()),
                                             gzip.decompress(artifact(dest / ("sample." + suffix)).read_bytes()))

    def test_family_cap_serial_and_eight_worker_outputs_match(self):
        catalog, task, context = example()
        tasks = [replace(task, final_id="c" + str(i)) for i in range(32)]
        with harmony.candidate_results(tasks, context, catalog, 1) as results:
            expected = list(results)
        with harmony.candidate_results(tasks, context, catalog, 8, 1) as results:
            self.assertEqual(list(results), expected)


if __name__ == "__main__":
    unittest.main()
