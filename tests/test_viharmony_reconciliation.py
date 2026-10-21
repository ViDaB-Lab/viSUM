"""Taxonomy-only integration: retention, source evidence, and spawn safety."""
from argparse import Namespace
import copy
import csv
import gzip
import json
import shutil
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import run_viharmony as harmony
from reporting_test_paths import artifact

BUNDLE = ROOT / "assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz"
MSL = ROOT / "assets/ICTV_VMR_MSL41.csv"


def write_table(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fixture(directory, audit="full"):
    sample, ids = "sample", ["c1", "c2"]
    normalized = directory / "normalized.fasta"
    normalized.write_text("".join(">" + sid + "\n" + "A" * 1200 + "\n" for sid in ids))
    header, gate, regions, evidence = [directory / (name + ".tsv") for name in ("header", "gate", "regions", "evidence")]
    write_table(header, ["sample_id", "sequence_id", "original_id", "original_header", "length"],
                [dict(sample_id=sample, sequence_id=sid, original_id=sid, original_header=sid, length="1200") for sid in ids])
    write_table(gate, ["sequence_id", "advance_to_refinement", "discovery_status"],
                [dict(sequence_id=sid, advance_to_refinement="true", discovery_status="ambiguous") for sid in ids])
    write_table(regions, ["sequence_id", "parent_sequence_id", "record_type", "coordinates", "original_length", "refined_length"],
                [dict(sequence_id=sid, parent_sequence_id="", record_type="input_contig", coordinates="", original_length="1200", refined_length="1200") for sid in ids])
    rows = [
        dict(sequence_id="c1", tool="gianthunter", classification="virus", evidence_strength="qualified",
             r__Realm="r__Varidnaviria", c__Class="c__Megaviricetes"),
        dict(sequence_id="c1", tool="vcontact3", classification="virus", evidence_strength="",
             r__Realm="r__Duplodnaviria", c__Class="c__Naldaviricetes"),
        dict(sequence_id="c1", tool="deepmicroclass2", classification="cellular", evidence_strength="qualified"),
        dict(sequence_id="c1", tool="vicat", classification="cellular", evidence_strength="qualified"),
        dict(sequence_id="c2", tool="genomad", classification="virus", evidence_strength="strong", r__Realm="r__Riboviria"),
        dict(sequence_id="c2", tool="vitap", classification="virus", evidence_strength="qualified",
             vitap_confidence_level="High-confidence", r__Realm="r__Riboviria", g__Genus="g__Gammacystovirus", s__Species="s__Gammacystovirus phiYY"),
        dict(sequence_id="c2", tool="vicat", classification="virus", evidence_strength="qualified",
             r__Realm="r__Riboviria", g__Genus="g__Cystovirus", s__Species="s__Cystovirus phiYY"),
    ]
    rows = [dict(sample_id=sample, parent_sequence_id="", length="1200", **row) for row in rows]
    columns = sorted({name for row in rows for name in row})
    write_table(evidence, columns, rows)
    return Namespace(sample_id=sample, input_type="dna", normalized_fasta=normalized,
        refined_fasta=normalized, header_map=header, discovery_gate=gate, region_map=regions,
        ictv_msl=MSL, evidence=[evidence], vcontact3_groups=[], audit_mode=audit, workers=1,
        vcontact3_min_taxonomy_length=1000, rna_pair_homology_floor=False,
        disable_rna_checkv_exception=False, output_prefix=directory / sample)


class ReconciliationIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog, _ = harmony.load_msl(MSL)

    def context(self, enabled=True, full=True):
        return harmony.CandidateContext("sample", "dna", 1000, False, False, False, full,
            str(BUNDLE) if enabled else None, harmony.sha256(BUNDLE) if enabled else None)

    def test_taxonomy_copies_preserve_scores_strength_and_originals(self):
        rows = [dict(tool="vcontact3", classification="virus", length="2000", score="1.0",
                     evidence_strength="", r__Realm="r__Varidnaviria", c__Class="c__Naldaviricetes")]
        original = copy.deepcopy(rows)
        audit = []
        result = harmony.decide_candidate_taxonomy(rows, self.catalog, 1000, self.context(), audit)
        self.assertEqual(rows, original)
        self.assertEqual(result.strict["class"], "Naldaviricetes")
        self.assertEqual(result.strict["realm"], "")
        self.assertEqual(json.loads(audit[0]["original_taxonomy"])["realm"], "Varidnaviria")
        self.assertEqual(json.loads(audit[0]["unresolved"])[0]["reason"], "target_higher_rank_unassigned")

    def test_below_length_vcontact3_cannot_be_rescued_by_reconciliation(self):
        rows = [dict(tool="vcontact3", classification="virus", length="999",
                     r__Realm="r__Varidnaviria", c__Class="c__Naldaviricetes")]
        audit = []
        result = harmony.decide_candidate_taxonomy(rows, self.catalog, 1000, self.context(), audit)
        self.assertEqual(audit, [])
        self.assertEqual(result.strict["class"], "")

    def test_actual_assigned_conflict_stays_rejected_and_audited(self):
        rows = [dict(tool="vcontact3", classification="virus", length="2000",
                     r__Realm="r__Varidnaviria", c__Class="c__Caudoviricetes")]
        audit = []
        result = harmony.decide_candidate_taxonomy(rows, self.catalog, 1000, self.context(), audit)
        self.assertEqual(audit[0]["status"], "incompatible")
        self.assertEqual(result.strict["class"], "")

    def test_disabled_adapter_is_identical_to_legacy_taxonomy_engine(self):
        rows = [dict(tool="genomad", classification="virus", evidence_strength="strong", r__Realm="r__Riboviria")]
        audit = []
        self.assertEqual(harmony.decide_candidate_taxonomy(rows, self.catalog, 1000, self.context(False), audit),
                         harmony.decide_taxonomy(rows, self.catalog, 1000))
        self.assertEqual(audit, [])

    def test_loader_rejects_changed_pinned_bundle(self):
        with self.assertRaisesRegex(ValueError, "changed before loading"):
            harmony.load_reconciliation(str(BUNDLE), "0" * 64)

    def test_complete_run_retention_and_original_audit_identical(self):
        from replay_integrated_ictv_reconciliation import compare
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = fixture(root)
            baseline, integrated = root / "baseline", root / "integrated"
            baseline.mkdir(); integrated.mkdir()
            before = {str(path): harmony.sha256(path) for path in (args.normalized_fasta, args.header_map, args.region_map, *args.evidence)}
            args.output_prefix = baseline / "sample"
            harmony.run(args)
            args.ictv_reconciliation = BUNDLE
            args.output_prefix = integrated / "sample"
            harmony.run(args)
            count, changes = compare(baseline, integrated, "sample")
            self.assertEqual(count, 2)
            self.assertTrue(changes)
            self.assertTrue(all(harmony.sha256(Path(p)) == h for p, h in before.items()))
            manifest = json.loads(artifact(integrated / "sample.harmonizer_manifest.json").read_text())
            self.assertEqual(manifest["taxonomy_reconciliation"]["sha256"], harmony.sha256(BUNDLE))
            self.assertIn("taxonomy_reconciliation_audit", manifest["outputs"])

    def test_compact_and_none_still_preserve_changed_call_provenance(self):
        for mode in ("compact", "none"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                args = fixture(root, mode)
                args.ictv_reconciliation = BUNDLE
                harmony.run(args)
                with gzip.open(artifact(root / "sample.taxonomy_reconciliation_audit.tsv.gz"), "rt") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertTrue(rows)
                self.assertTrue(all(r["status"] != "unchanged" for r in rows))
                self.assertFalse((root / "sample.evidence_audit.tsv.gz").exists())

    def test_bundle_target_must_match_configured_msl(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = fixture(root)
            args.ictv_reconciliation = BUNDLE
            args.ictv_msl = root / "other_msl.csv"
            args.ictv_msl.write_text("Realm,Kingdom,Phylum,Class,Order,Family,Genus,Species\nRiboviria,,,,,,,Fake species\n")
            with self.assertRaisesRegex(ValueError, "target differs"):
                harmony.run(args)

    def test_each_run_revalidates_manifest_even_after_cache_is_populated(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = fixture(root)
            copied = root / BUNDLE.name
            manifest_path = copied.with_name(copied.name + ".manifest.json")
            shutil.copyfile(BUNDLE, copied)
            shutil.copyfile(BUNDLE.with_name(BUNDLE.name + ".manifest.json"), manifest_path)
            args.ictv_reconciliation = copied
            harmony.run(args)
            manifest = json.loads(manifest_path.read_text())
            manifest["bundle_sha256"] = "0" * 64
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                harmony.run(args)

    def test_reconciled_spawn_results_match_serial(self):
        rows = [dict(tool="vcontact3", classification="virus", length="2000", sequence_id="c1",
                     r__Realm="r__Varidnaviria", c__Class="c__Naldaviricetes")]
        task = harmony.CandidateTask("c1", 2000,
            dict(record_type="input_contig", coordinates="", original_length="2000"), "c1",
            dict(original_id="c1", original_header="c1"), "ambiguous", rows, [], ())
        context = self.context()
        with harmony.candidate_results([task, task], context, self.catalog, 1) as results:
            serial = list(results)
        with harmony.candidate_results([task, task], context, self.catalog, 2, 1) as results:
            parallel = list(results)
        self.assertEqual(serial, parallel)

    def test_cli_two_workers_with_bundle_matches_serial_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = fixture(root)
            serial, parallel = root / "serial", root / "parallel"
            serial.mkdir(); parallel.mkdir()
            args.ictv_reconciliation = BUNDLE
            args.output_prefix = serial / "sample"
            harmony.run(args)
            args.output_prefix = parallel / "sample"
            args.workers = 2
            command = [sys.executable, "-B", str(ROOT / "bin/run_viharmony.py")]
            for name, value in vars(args).items():
                option = "--" + name.replace("_", "-")
                if isinstance(value, bool):
                    if value:
                        command.append(option)
                elif isinstance(value, list):
                    command.extend([option, *map(str, value)])
                else:
                    command.extend([option, str(value)])
            process = subprocess.run(command, text=True, capture_output=True, timeout=120)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            self.assertIn("workers=2", process.stdout)
            left = json.loads(artifact(serial / "sample.harmonizer_manifest.json").read_text())
            right = json.loads(artifact(parallel / "sample.harmonizer_manifest.json").read_text())
            self.assertEqual(left["taxonomy_reconciliation"], right["taxonomy_reconciliation"])
            self.assertEqual(left["outputs"].keys(), right["outputs"].keys())
            for record in left["outputs"].values():
                filename = record["file"]
                opener = gzip.decompress if filename.endswith(".gz") else lambda b: b
                self.assertEqual(opener((serial / filename).read_bytes()), opener((parallel / filename).read_bytes()), filename)


if __name__ == "__main__":
    unittest.main()
