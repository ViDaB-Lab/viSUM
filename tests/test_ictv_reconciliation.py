import copy
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from ictv_reconciliation import RANKS, SCHEMA, Reconciler
from build_ictv_reconciliation_bundle import read_export, verified_sources, forward_history_link

BUNDLE = ROOT / "assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz"


def toy_bundle():
    """Explicitly synthetic tree: two families, with two species in genus G."""
    rows = [("0", "0", "tree", "root"), ("1", "0", "realm", "R"),
            ("2", "1", "order", "O"), ("3", "2", "family", "F"),
            ("4", "3", "genus", "G"), ("5", "4", "species", "G a"),
            ("6", "4", "species", "G b"), ("7", "2", "family", "OtherF"),
            ("8", "7", "genus", "H"), ("9", "8", "species", "H c")]
    aliases = {}
    for rank, name, targets in [("genus", "OldGenus", ["4", "8"]),
                                ("species", "OldSplitSpecies", ["5", "6"]),
                                ("family", "OldPromotedFamily", ["2"]),
                                ("species", "OldMergedSpecies", ["5"])]:
        aliases[rank + "\t" + name] = [{"source_node_id": "old-" + name,
            "ictv_id": "old-" + name, "release": 39, "targets": targets,
            "proposals": ["synthetic-test-only.zip"]}]
    return {"schema": SCHEMA, "target_release": "MSL41.v1", "source_releases": [39, 41],
            "nodes": [{"id": i, "parent": p, "rank": r, "name": n} for i, p, r, n in rows],
            "aliases": aliases}


class ResolverUnitTests(unittest.TestCase):
    def setUp(self):
        self.resolver = Reconciler(toy_bundle())

    def test_current_name_is_not_changed_and_does_not_infer_descendants(self):
        r = self.resolver.resolve({"family": "F"})
        self.assertEqual(r.status, "unchanged")
        self.assertEqual(r.normalized, {"realm": "R", "order": "O", "family": "F"})

    def test_split_is_not_arbitrarily_resolved(self):
        r = self.resolver.resolve({"genus": "OldGenus"})
        self.assertEqual(r.status, "ambiguous")
        self.assertEqual(r.normalized, {"realm": "R", "order": "O"})

    def test_context_can_resolve_split(self):
        r = self.resolver.resolve({"family": "F", "genus": "OldGenus"})
        self.assertEqual(r.status, "translated")
        self.assertEqual(r.normalized["genus"], "G")
        self.assertNotIn("species", r.normalized)

    def test_species_split_reports_only_shared_ancestor(self):
        r = self.resolver.resolve({"species": "OldSplitSpecies"})
        self.assertEqual(r.status, "ambiguous")
        self.assertEqual(r.normalized["genus"], "G")
        self.assertNotIn("species", r.normalized)

    def test_promotion_does_not_invent_new_family(self):
        r = self.resolver.resolve({"family": "OldPromotedFamily"})
        self.assertEqual(r.normalized, {"realm": "R", "order": "O"})
        self.assertEqual(r.status, "translated")

    def test_documented_merge_can_have_one_destination(self):
        r = self.resolver.resolve({"species": "OldMergedSpecies"})
        self.assertEqual(r.normalized["species"], "G a")
        self.assertEqual(r.status, "translated")

    def test_demotion_to_unreported_subrank_stops_at_formal_ancestor(self):
        b = toy_bundle()
        b["nodes"].append({"id": "10", "parent": "3", "rank": "subfamily", "name": "NewSubfamily"})
        b["nodes"][4]["parent"] = "10"
        b["aliases"]["family\tOldDemotedFamily"] = [{"release": 39, "targets": ["10"]}]
        r = Reconciler(b).resolve({"family": "OldDemotedFamily"})
        self.assertEqual(r.normalized["family"], "F")
        self.assertNotIn("genus", r.normalized)
        self.assertEqual(r.status, "translated")

    def test_unknown_suffix_preserves_only_verified_prefix(self):
        r = self.resolver.resolve({"family": "F", "genus": "NotInHistory", "species": "G a"})
        self.assertEqual(r.status, "partial")
        self.assertEqual(r.normalized["family"], "F")
        self.assertNotIn("genus", r.normalized)
        self.assertNotIn("species", r.normalized)

    def test_unknown_ancestor_does_not_get_guessed_from_species(self):
        r = self.resolver.resolve({"realm": "Unknown", "species": "G a"})
        self.assertEqual(r.normalized, {})
        self.assertEqual(r.status, "partial")

    def test_incompatible_current_claims_are_not_silently_repaired(self):
        r = self.resolver.resolve({"family": "OtherF", "genus": "G"})
        self.assertEqual(r.status, "incompatible")
        self.assertEqual(r.normalized, {})

    def orphan_resolver(self):
        b = toy_bundle()
        b["nodes"].extend([
            {"id": "10", "parent": "0", "rank": "class", "name": "OrphanClass"},
            {"id": "11", "parent": "10", "rank": "family", "name": "OrphanFamily"},
            {"id": "12", "parent": "11", "rank": "genus", "name": "OrphanGenus"},
            {"id": "13", "parent": "0", "rank": "realm", "name": "OtherR"},
            {"id": "14", "parent": "0", "rank": "phylum", "name": "OtherP"},
            {"id": "15", "parent": "0", "rank": "kingdom", "name": "OtherK"},
        ])
        return Reconciler(b)

    def test_unassigned_higher_rank_preserves_lower_taxonomy_and_original(self):
        resolver = self.orphan_resolver()
        row = {"realm": "R", "class": "OrphanClass", "family": "OrphanFamily"}
        original = dict(row)
        bundle = copy.deepcopy(resolver.bundle)
        r = resolver.resolve(row)
        self.assertEqual(r.status, "partial")
        self.assertEqual(r.normalized, {"class": "OrphanClass", "family": "OrphanFamily"})
        self.assertEqual(r.unresolved[0]["reason"], "target_higher_rank_unassigned")
        self.assertEqual(r.unresolved[0]["rank"], "realm")
        self.assertEqual(r.original, original)
        self.assertEqual(row, original)
        self.assertEqual(resolver.bundle, bundle)
        self.assertEqual(r.mappings[0]["disposition"], "not_used_target_higher_rank_unassigned")
        self.assertEqual(r.mappings[0]["target_node_ids"], ["1"])

    def test_multiple_unassigned_higher_ranks_are_flagged_individually(self):
        r = self.orphan_resolver().resolve({"realm": "R", "kingdom": "OtherK",
            "phylum": "OtherP", "class": "OrphanClass", "family": "OrphanFamily"})
        self.assertEqual(r.normalized, {"class": "OrphanClass", "family": "OrphanFamily"})
        self.assertEqual([x["rank"] for x in r.unresolved], ["realm", "kingdom", "phylum"])
        self.assertTrue(all(x["reason"] == "target_higher_rank_unassigned" for x in r.unresolved))

    def test_unknown_higher_label_is_not_treated_as_verified_unassigned(self):
        r = self.orphan_resolver().resolve({"realm": "UnknownR", "class": "OrphanClass"})
        self.assertEqual(r.normalized, {})
        self.assertEqual(r.status, "partial")
        self.assertEqual(r.unresolved[0]["reason"], "no_verified_current_mapping")

    def test_lower_branch_conflict_is_not_masked_by_unassigned_realm(self):
        b = self.orphan_resolver().bundle
        b["nodes"].append({"id": "17", "parent": "2", "rank": "class", "name": "AssignedClass"})
        next(n for n in b["nodes"] if n["name"] == "OtherF")["parent"] = "17"
        r = Reconciler(b).resolve({"realm": "R", "class": "OrphanClass", "family": "OtherF"})
        self.assertEqual(r.status, "incompatible")
        self.assertEqual(r.normalized, {})

    def test_actual_realm_conflict_survives_an_unassigned_intermediate_rank(self):
        r = self.orphan_resolver().resolve({"realm": "OtherR", "class": "OrphanClass", "family": "F"})
        self.assertEqual(r.status, "incompatible")
        self.assertEqual(r.normalized, {})

    def test_lower_unassigned_taxon_cannot_hide_conflict_between_upper_claims(self):
        r = self.orphan_resolver().resolve({"realm": "OtherR", "order": "O", "family": "OrphanFamily"})
        self.assertEqual(r.status, "incompatible")
        self.assertEqual(r.normalized, {})

    def test_mixed_assigned_and_unassigned_split_is_not_used_to_drop_ancestor(self):
        b = self.orphan_resolver().bundle
        b["aliases"]["genus\tMixedOldGenus"] = [{"release": 39, "targets": ["4", "12"]}]
        r = Reconciler(b).resolve({"realm": "OtherR", "genus": "MixedOldGenus"})
        self.assertEqual(r.status, "incompatible")
        self.assertEqual(r.normalized, {})

    def test_unknown_suffix_barrier_remains_after_unassigned_rank_recovery(self):
        r = self.orphan_resolver().resolve({"realm": "R", "class": "OrphanClass",
                                          "genus": "UnknownG", "species": "G a"})
        self.assertEqual(r.normalized, {"class": "OrphanClass"})
        self.assertEqual([x["reason"] for x in r.unresolved],
                         ["target_higher_rank_unassigned", "no_verified_current_mapping"])

    def test_ambiguous_unassigned_lower_options_remain_ambiguous_in_audit(self):
        b = self.orphan_resolver().bundle
        b["nodes"].append({"id": "16", "parent": "11", "rank": "genus", "name": "SecondOrphanGenus"})
        b["aliases"]["genus\tOldOrphanGenus"] = [{"release": 39, "targets": ["12", "16"]}]
        r = Reconciler(b).resolve({"realm": "R", "genus": "OldOrphanGenus"})
        self.assertEqual(r.status, "partial")
        self.assertEqual(r.normalized, {"class": "OrphanClass", "family": "OrphanFamily"})
        self.assertEqual(r.mappings[-1]["compatible_target_node_ids"], ["12", "16"])

    def test_evidence_wrapper_only_removes_unsupported_taxonomic_rank(self):
        from validate_ictv_reconciliation_replay import normalize_row
        import run_viharmony
        row = {"tool": "vcontact3", "classification": "virus", "score": "0.95",
               "evidence_strength": "qualified", "r__Realm": "r__R",
               "c__Class": "c__OrphanClass", "f__Family": "f__OrphanFamily"}
        before = dict(row)
        copied, resolution = normalize_row(row, self.orphan_resolver())
        self.assertEqual(row, before)
        self.assertEqual(copied["r__Realm"], "r__unclassified")
        self.assertEqual(copied["f__Family"], "f__OrphanFamily")
        for key in ("tool", "classification", "score", "evidence_strength"):
            self.assertEqual(copied[key], row[key])
        self.assertEqual(resolution.status, "partial")

    def test_no_named_taxonomy_is_distinct_from_invalid_taxonomy(self):
        self.assertEqual(self.resolver.resolve({"family": "unclassified"}).status, "no_named_taxonomy")

    def test_no_fuzzy_name_matching(self):
        self.assertEqual(self.resolver.resolve({"genus": "g"}).status, "partial")

    def test_source_release_must_exist(self):
        with self.assertRaisesRegex(ValueError, "Source release"):
            self.resolver.resolve({"family": "F"}, 999)

    def test_known_source_release_restricts_alias_observations(self):
        r = self.resolver.resolve({"genus": "OldGenus"}, 41)
        self.assertEqual(r.status, "partial")
        self.assertEqual(self.resolver.resolve({"genus": "OldGenus"}, 39).status, "ambiguous")

    def test_original_input_and_bundle_are_not_mutated(self):
        original = {"genus": "OldGenus"}
        bundle_before = copy.deepcopy(self.resolver.bundle)
        r = self.resolver.resolve(original)
        self.assertEqual(original, {"genus": "OldGenus"})
        self.assertEqual(r.original, original)
        self.assertEqual(self.resolver.bundle, bundle_before)

    def test_deterministic(self):
        self.assertEqual(self.resolver.resolve({"genus": "OldGenus"}), self.resolver.resolve({"genus": "OldGenus"}))

    def test_cycle_rejected(self):
        b = toy_bundle(); b["nodes"][1]["parent"] = "3"
        with self.assertRaisesRegex(ValueError, "Cycle"):
            Reconciler(b)

    def test_missing_parent_rejected(self):
        b = toy_bundle(); b["nodes"][1]["parent"] = "missing"
        with self.assertRaisesRegex(ValueError, "Missing target parent"):
            Reconciler(b)

    def test_foreign_mapping_target_rejected(self):
        b = toy_bundle(); b["aliases"]["genus\tOldGenus"][0]["targets"] = ["missing"]
        with self.assertRaisesRegex(ValueError, "outside target"):
            Reconciler(b)

    def test_unsupported_schema_rejected(self):
        b = toy_bundle(); b["schema"] = "future-incompatible-schema"
        with self.assertRaisesRegex(ValueError, "schema"):
            Reconciler(b)

    def test_empty_target_rejected(self):
        b = toy_bundle(); b["nodes"] = []
        with self.assertRaisesRegex(ValueError, "Empty target"):
            Reconciler(b)

    def test_corrupt_bundle_checksum_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "test.json.gz"
            p.write_bytes(gzip.compress(json.dumps(toy_bundle()).encode()))
            p.with_name(p.name + ".manifest.json").write_text(json.dumps({"bundle_sha256": "wrong"}))
            with self.assertRaisesRegex(ValueError, "checksum"):
                Reconciler.load(p)

    def test_export_null_and_backslash_escaping(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "source.tsv"
            p.write_text('name\tnotes\n"Example"\t\\N\n"Quoted \\"value\\""\t"ok"\n')
            rows = list(read_export(p))
            self.assertEqual(rows[0]["notes"], "")
            self.assertEqual(rows[1]["name"], 'Quoted "value"')

    def test_malformed_export_not_silently_truncated(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "source.tsv"
            p.write_text("a\tb\n1\t2\t3\n")
            with self.assertRaisesRegex(ValueError, "Malformed"):
                list(read_export(p))

    def test_source_snapshot_checksum_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "tiny.txt"
            p.write_bytes(b"corrupt")
            with patch.dict("build_ictv_reconciliation_bundle.SOURCE_BLOBS", {"tiny.txt": "wrong"}, clear=True):
                with self.assertRaisesRegex(ValueError, "pinned official Git blob"):
                    verified_sources(Path(directory))

    def test_manifest_target_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "test.json.gz"
            data = gzip.compress(json.dumps(toy_bundle()).encode())
            p.write_bytes(data)
            p.with_name(p.name + ".manifest.json").write_text(json.dumps({
                "bundle_sha256": hashlib.sha256(data).hexdigest(), "target_release": "MSL999"}))
            with self.assertRaisesRegex(ValueError, "target release mismatch"):
                Reconciler.load(p)

    def test_reverse_or_recreated_taxa_are_not_equivalence_links(self):
        self.assertTrue(forward_history_link({"rev_count": "0", "is_recreated": "0"}))
        self.assertFalse(forward_history_link({"rev_count": "1", "is_recreated": "0"}))
        self.assertFalse(forward_history_link({"rev_count": "0", "is_recreated": "1"}))

    def test_evidence_wrapper_preserves_scores_classification_and_original(self):
        from validate_ictv_reconciliation_replay import normalize_row
        import run_viharmony
        row = {"tool": "vitap", "classification": "virus", "score": "0.95",
               "vitap_confidence_level": "High-confidence", "evidence_strength": "qualified",
               run_viharmony.RANK_COLUMNS["genus"]: "g__OldGenus"}
        before = dict(row)
        copied, resolution = normalize_row(row, self.resolver)
        self.assertEqual(row, before)
        for key in ("tool", "classification", "score", "vitap_confidence_level", "evidence_strength"):
            self.assertEqual(copied[key], row[key])
        self.assertEqual(resolution.status, "ambiguous")
        self.assertEqual(copied[run_viharmony.RANK_COLUMNS["genus"]], "g__unclassified")


class OfficialSnapshotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.resolver = Reconciler.load(BUNDLE)

    def test_documented_species_rename_from_official_identity(self):
        r = self.resolver.resolve({"realm": "Riboviria", "family": "Cystoviridae", "genus": "Cystovirus", "species": "Cystovirus phiYY"}, 39)
        self.assertEqual(r.normalized["species"], "Gammacystovirus phiYY")
        self.assertEqual(r.normalized["genus"], "Gammacystovirus")
        self.assertEqual(r.status, "translated")
        self.assertTrue(any("2024.043B" in p for x in r.mappings for h in x["history"] for p in h["proposals"]))

    def test_old_realm_alone_does_not_choose_surviving_renamed_branch(self):
        r = self.resolver.resolve({"realm": "Monodnaviria"}, 40)
        self.assertEqual(r.status, "ambiguous")
        self.assertNotIn("realm", r.normalized)
        targets = r.mappings[0]["compatible_target_node_ids"]
        self.assertEqual({self.resolver.nodes[n]["name"] for n in targets}, {"Efunaviria", "Floreoviria", "Pleomoviria", "Volvereviria"})

    def test_old_realm_with_kingdom_context(self):
        for kingdom, realm in [("Shotokuvirae", "Floreoviria"), ("Loebvirae", "Efunaviria"), ("Sangervirae", "Volvereviria"), ("Trapavirae", "Pleomoviria")]:
            with self.subTest(kingdom=kingdom):
                r = self.resolver.resolve({"realm": "Monodnaviria", "kingdom": kingdom}, 40)
                self.assertEqual(r.normalized["realm"], realm)
                self.assertEqual(r.status, "translated")

    def test_official_rank_promotion_preserves_supported_rank(self):
        r = self.resolver.resolve({"family": "Autographiviridae"}, 39)
        self.assertEqual(r.normalized["order"], "Autographivirales")
        self.assertNotIn("family", r.normalized)

    def test_complete_current_lineage_unchanged(self):
        row = {"realm": "Riboviria", "kingdom": "Orthornavirae", "phylum": "Duplornaviricota", "class": "Vidaverviricetes", "order": "Mindivirales", "family": "Cystoviridae", "genus": "Gammacystovirus", "species": "Gammacystovirus phiYY"}
        r = self.resolver.resolve(row)
        self.assertEqual(r.normalized, row)
        self.assertEqual(r.status, "unchanged")

    def test_current_ancestry_matches_bundled_msl_exactly(self):
        import csv
        with (ROOT / "assets/ICTV_VMR_MSL41.csv").open(encoding="utf-8-sig", newline="") as handle:
            expected = {tuple(row[r.title()] for r in RANKS) for row in csv.DictReader(handle)}
        actual = {tuple(self.resolver.lineages[nid].get(r, "") for r in RANKS) for nid, n in self.resolver.nodes.items() if n["rank"] == "species"}
        self.assertEqual(actual, expected)
        self.assertEqual(len(actual), 17554)

    def test_naldaviricetes_unassigned_ancestors_do_not_discard_lower_calls(self):
        lower = {"class": "Naldaviricetes", "order": "Lefavirales",
                 "family": "Baculoviridae", "genus": "Gammabaculovirus"}
        for realm in ("Varidnaviria", "Duplodnaviria"):
            with self.subTest(realm=realm):
                original = {"realm": realm, **lower}
                r = self.resolver.resolve(original)
                self.assertEqual(r.original, original)
                self.assertEqual(r.normalized, lower)
                self.assertEqual(r.status, "partial")
                self.assertEqual(r.unresolved[0]["reason"], "target_higher_rank_unassigned")
        self.assertEqual(self.resolver.resolve(lower).status, "unchanged")

    def test_assigned_caudoviricetes_realm_conflict_remains_incompatible(self):
        r = self.resolver.resolve({"realm": "Varidnaviria", "class": "Caudoviricetes"})
        self.assertEqual(r.status, "incompatible")
        self.assertEqual(r.normalized, {})

    def test_known_source_release_retained_when_recovering_unassigned_ancestors(self):
        r = self.resolver.resolve({"realm": "Varidnaviria", "class": "Naldaviricetes"}, 40)
        self.assertEqual(r.normalized, {"class": "Naldaviricetes"})
        self.assertEqual(r.source_release, 40)
        self.assertTrue(any(h["release"] == 40 for m in r.mappings for h in m["history"]))

    def test_entire_target_release_handles_unassigned_formal_ranks(self):
        from audit_ictv_unassigned_ranks import audit_release
        report = audit_release(self.resolver)
        self.assertGreater(report["synthetic_recovery_cases_verified"], 0)
        self.assertEqual(report["unchanged_valid_lineages_verified"], report["taxa_with_missing_higher_formal_ranks"])
        self.assertEqual(report["synthetic_recovery_cases_verified"] + report["synthetic_assigned_conflicts_verified"],
                         sum(report["by_missing_rank"].values()))


if __name__ == "__main__":
    unittest.main()
