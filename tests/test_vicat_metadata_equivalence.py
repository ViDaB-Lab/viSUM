"""Small offline vOTU partition-equivalence regression; no acquisition/build."""
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import prepare_vicat_metavr_sources as adapter
from audit_vicat_metadata_equivalence import audit_metadata_equivalence


class PartitionAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.base = [
            ('U1', 'old-A', 'd__Viruses;o__Ortervirales;f__Retroviridae', 'method'),
            ('U2', 'old-A', 'd__Viruses;o__Ortervirales;f__Retroviridae', 'method'),
            ('U3', 'old-B', 'd__Viruses;o__Ortervirales;f__Retroviridae', 'method'),
            ('U4', 'old-B', 'd__Viruses;o__Ortervirales;f__Retroviridae', 'method')]

    def write(self, name, rows):
        path = self.root / (name + '.tsv.gz')
        with gzip.open(path, 'wt', encoding='utf-8') as stream:
            stream.write('uvig\tvotu\tictv_taxonomy\tictv_taxonomy_method\n')
            for row in rows:
                stream.write('\t'.join(row) + '\n')
        return path

    def audit(self, new, original=None):
        return audit_metadata_equivalence(self.write('old', self.base if original is None else original),
            self.write('new', new), self.root / 'work', 4, adapter.connect)

    def test_exact_match_and_bijective_renaming(self):
        self.assertEqual(self.audit(self.base)['status'], 'EXACT_MATCH')
        renamed = [(row[0], {'old-A': 'new-Z', 'old-B': 'new-Y'}[row[1]], *row[2:]) for row in reversed(self.base)]
        result = self.audit(renamed)
        self.assertEqual(result['status'], 'EQUIVALENT_VOTU_RELABELING')
        self.assertTrue(result['builder_semantics_equivalent'])
        self.assertFalse(result['exact_builder_metadata_match'])
        self.assertEqual(result['votu_label_changed_rows'], 4)
        self.assertEqual(result['split_original_groups'], 0)
        self.assertEqual(result['merged_portal_groups'], 0)
        self.assertEqual(result['mapping_pairs'], 2)

    def test_split_merge_and_balanced_regrouping_never_pass(self):
        cases = (
            [('U1', 'X', *self.base[0][2:]), ('U2', 'Y', *self.base[1][2:]), self.base[2], self.base[3]],
            [(row[0], 'merged', *row[2:]) for row in self.base],
            [(row[0], ('X' if row[0] in ('U1', 'U3') else 'Y'), *row[2:]) for row in self.base])
        for rows in cases:
            with self.subTest(rows=rows):
                result = self.audit(rows)
                self.assertFalse(result['builder_semantics_equivalent'])
                self.assertGreater(result['split_original_groups'] + result['merged_portal_groups'], 0)

    def test_missing_duplicate_and_changed_taxonomy_method_do_not_pass(self):
        for rows in (
            [('another-UViG', *self.base[0][1:])] + self.base[1:],
            [self.base[0], self.base[0], self.base[2], self.base[3]],
            [(self.base[0][0], self.base[0][1], 'd__Viruses;f__Changed', self.base[0][3])] + self.base[1:],
            [(*self.base[0][:3], 'changed-method')] + self.base[1:]):
            with self.subTest(rows=rows):
                self.assertFalse(self.audit(rows)['builder_semantics_equivalent'])

    def test_exact_uploaded_twenty_examples_are_renaming_not_rank_changes(self):
        source = json.loads((ROOT / 'tests/fixtures/metavr_linux_votu_examples.json').read_text())
        prefixes = [('realm', 'r__'), ('kingdom', 'k__'), ('phylum', 'p__'), ('tax_class', 'c__'),
                    ('tax_order', 'o__'), ('family', 'f__'), ('genus', 'g__'), ('species', 's__')]
        def rows(side):
            return [(entry['uvig'], entry[side + '_votu_key'],
                     'd__Viruses;' + ';'.join(prefix + entry[side + '_' + key] for key, prefix in prefixes if entry[side + '_' + key] is not None),
                     entry[side + '_taxonomy_method'] or '\\N') for entry in source['examples']]
        old, new = self.write('reported-old', rows('original')), self.write('reported-new', rows('portal'))
        result = audit_metadata_equivalence(old, new, self.root / 'work', 20, adapter.connect)
        self.assertEqual(result['status'], 'EQUIVALENT_VOTU_RELABELING')
        self.assertEqual(result['rank_changed_rows'], 0)
        self.assertEqual(result['method_changed_rows'], 0)

    def test_actual_lookup_and_audit_outputs_unchanged_after_bijective_relabeling(self):
        import build_vicat_taxonomy_lookup as builder
        members = self.root / 'members.tsv.gz'
        with gzip.open(members, 'wt', encoding='utf-8') as stream:
            stream.write('U1|p1\tU1|p1\nU1|p1\tU2|p2\nU1|p1\tU3|p3\nU4|p4\tU4|p4\n')
        renamed = [(row[0], {'old-A': 'new-Z', 'old-B': 'new-Y'}[row[1]], *row[2:]) for row in reversed(self.base)]
        results = []
        for label, rows in (('old', self.base), ('new', renamed)):
            output = self.root / (label + '-output')
            args = ['builder', '--cluster-members', str(members), '--metadata', str(self.write(label, rows)),
                '--output-dir', str(output), '--work-dir', str(self.root / (label + '-build-work')),
                '--threads', '1', '--memory-limit', '1GB', '--expected-member-count', '4', '--expected-representative-count', '2']
            with patch.object(sys, 'argv', args):
                builder.main()
            with adapter.connect(self.root) as connection:
                results.append([connection.execute('SELECT * FROM read_parquet(?) ORDER BY 1', [str(path)]).fetchall()
                                for path in sorted(output.glob('*.parquet'))])
        self.assertEqual(results[0], results[1])

if __name__ == '__main__':
    unittest.main()

