"""Tiny offline patch regression tests. No networks, proteins or database build."""
import argparse
import contextlib
import csv
import gzip
import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'bin'), str(ROOT / 'tests')]
import patch_vicat_metavr_metadata as p
from audit_vicat_metadata_equivalence import audit_metadata_equivalence
try:
    import metadata_adapter as adapter
except ImportError:
    import prepare_vicat_metavr_sources as adapter


class PatchTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.work = self.root / 'work'
        self.work.mkdir()
        self.old_rows = [('U1', 'A'), ('U2', 'A'), ('U3', 'B'), ('U4', 'B')]
        # Both original groups split; portal Y merges members of both originals.
        self.new_rows = [('U1', 'X'), ('U2', 'Y'), ('U3', 'Y'), ('U4', 'B')]
        self.old = self.write('original.tsv.gz', self.old_rows)
        self.new = self.write('converted.tsv.gz', self.new_rows)

    def write(self, name, rows, family='Retroviridae', method='tool', optional='kept'):
        path = self.root / name
        with gzip.open(path, 'wt', encoding='utf-8', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(('uvig', 'votu', 'ictv_taxonomy', 'ictv_taxonomy_method', 'source'))
            for uvig, votu in rows:
                writer.writerow((uvig, votu, 'd__Viruses;o__Ortervirales;f__' + family, method, optional))
        return path

    def export(self, expected=3, name='exports'):
        output = self.root / name
        manifest = p.export_patch(self.old, self.new, output, self.work, 4, expected,
                                  dict(portal_sha256='tiny-test-not-production'), adapter.connect)
        return output, manifest

    def apply(self, output):
        result = self.root / 'patched.tsv.gz'
        p.apply_patch_copy(self.new, output / 'votu_compatibility_patch.tsv',
            output / 'patch_manifest.json', result, self.work, adapter.connect)
        return result

    def test_split_merge_restoration_is_exact_and_inputs_unchanged(self):
        before = self.old.read_bytes(), self.new.read_bytes()
        output, manifest = self.export()
        self.assertEqual(manifest['before']['votu_changed_rows'], 3)
        self.assertEqual(manifest['before']['split_original_groups'], 2)
        self.assertEqual(manifest['before']['merged_portal_groups'], 1)
        result = self.apply(output)
        audit = audit_metadata_equivalence(self.old, result, self.work, 4, adapter.connect)
        self.assertEqual(audit['status'], 'EXACT_MATCH')
        self.assertEqual(p.verify_non_votu_columns(self.new, result, self.work, 4, adapter.connect)['changed_rows'], 0)
        self.assertEqual(before, (self.old.read_bytes(), self.new.read_bytes()))
        self.assertFalse(manifest['full_source_equivalence'])

    def test_exactly_1409_changes_exported_and_restored(self):
        self.old = self.write('many-original.tsv.gz', [('U' + str(i), 'old' + str(i)) for i in range(1410)])
        self.new = self.write('many-portal.tsv.gz', [('U' + str(i), ('new' if i < 1409 else 'old') + str(i)) for i in range(1410)])
        output = self.root / 'many-export'
        result = p.export_patch(self.old, self.new, output, self.work, 1410, 1409, {}, adapter.connect)
        self.assertEqual(result['expected_changes'], 1409)
        self.assertEqual(len((output / 'votu_compatibility_patch.tsv').read_text().splitlines()), 1410)
        patched = self.apply(output)
        self.assertEqual(audit_metadata_equivalence(self.old, patched, self.work, 1410, adapter.connect)['status'], 'EXACT_MATCH')

    def test_rank_and_method_changes_refused(self):
        for key, value in (('family', 'Changed'), ('method', 'Changed')):
            self.new = self.write(key + '.tsv.gz', self.new_rows, **{key: value})
            with self.assertRaisesRegex(ValueError, 'taxonomy/method changes'):
                self.export(name=key + '-exports')

    def test_missing_and_duplicate_uvigs_refused(self):
        for label, rows in (('missing', [('different', 'X')] + self.new_rows[1:]),
                            ('duplicate', [self.new_rows[0], self.new_rows[0], *self.new_rows[2:]])):
            self.new = self.write(label + '.tsv.gz', rows)
            with self.assertRaises(ValueError):
                self.export(name=label + '-exports')

    def test_unexpected_change_count_refused(self):
        with self.assertRaisesRegex(ValueError, 'Expected exactly 1409'):
            self.export(expected=1409)

    def test_tampered_patch_and_changed_source_fail_checksums(self):
        output, _ = self.export()
        patch_file = output / 'votu_compatibility_patch.tsv'
        original = patch_file.read_bytes()
        patch_file.write_bytes(original + b'bad')
        with self.assertRaisesRegex(ValueError, 'patch checksum mismatch'):
            self.apply(output)
        patch_file.write_bytes(original)
        self.new.write_bytes(self.new.read_bytes() + b'bad')
        with self.assertRaisesRegex(ValueError, 'metadata checksum differs'):
            self.apply(output)

    def reseal(self, output):
        path = output / 'patch_manifest.json'
        manifest = json.loads(path.read_text())
        manifest['patch_sha256'] = p.digest(output / 'votu_compatibility_patch.tsv')
        p.json_write(path, manifest)

    def test_wrong_portal_precondition_and_absent_id_rejected(self):
        for label, column, replacement in (('wrong', 1, 'unexpected-group'), ('absent', 0, 'absent-UViG')):
            output, _ = self.export(name=label)
            path = output / 'votu_compatibility_patch.tsv'
            with path.open(encoding='utf-8', newline='') as stream:
                rows = list(csv.reader(stream, delimiter='\t'))
            rows[1][column] = replacement
            with path.open('w', encoding='utf-8', newline='') as stream:
                csv.writer(stream, delimiter='\t', lineterminator='\n').writerows(rows)
            self.reseal(output)
            with self.assertRaisesRegex(ValueError, 'precondition failed'):
                self.apply(output)

    def test_duplicate_patch_id_rejected_even_with_resealed_test_manifest(self):
        output, _ = self.export()
        path = output / 'votu_compatibility_patch.tsv'
        with path.open(encoding='utf-8', newline='') as stream:
            rows = list(csv.reader(stream, delimiter='\t'))
        rows[2][0] = rows[1][0]
        with path.open('w', encoding='utf-8', newline='') as stream:
            csv.writer(stream, delimiter='\t', lineterminator='\n').writerows(rows)
        self.reseal(output)
        with self.assertRaisesRegex(ValueError, 'Duplicate UViGs'):
            self.apply(output)

    def test_original_null_group_restored_as_equivalent_normalized_key(self):
        self.old = self.write('null-original.tsv.gz', [('U1', '\\N'), *self.old_rows[1:]])
        output, _ = self.export()
        result = self.apply(output)
        self.assertEqual(audit_metadata_equivalence(self.old, result, self.work, 4, adapter.connect)['status'], 'EXACT_MATCH')

    def test_overwrites_refused(self):
        output, _ = self.export()
        with self.assertRaises(FileExistsError):
            self.export()
        self.apply(output)
        with self.assertRaises(FileExistsError):
            self.apply(output)

    def test_non_votu_optional_column_change_rejected(self):
        changed = self.write('changed-other.tsv.gz', self.new_rows, optional='changed')
        with self.assertRaisesRegex(ValueError, 'Non-vOTU fields changed'):
            p.verify_non_votu_columns(self.new, changed, self.work, 4, adapter.connect)

    def test_optional_nulls_and_empty_strings_survive_patching(self):
        for label, value in (('null', '\\N'), ('empty', '')):
            self.new = self.write(label + '-portal.tsv.gz', self.new_rows, optional=value)
            output, _ = self.export(name=label + '-exports')
            patched = self.root / (label + '-patched.tsv.gz')
            p.apply_patch_copy(self.new, output / 'votu_compatibility_patch.tsv',
                output / 'patch_manifest.json', patched, self.work, adapter.connect)
            self.assertEqual(p.verify_non_votu_columns(self.new, patched, self.work, 4, adapter.connect)['changed_rows'], 0)


@unittest.skipUnless((ROOT / 'PATCH_PACKAGE.json').is_file(), 'Package-only integration tests')
class PackageTests(unittest.TestCase):
    def test_complete_offline_runner_on_tiny_fixture(self):
        import run_metavr_patch_validation as runner
        helper = PatchTests()
        helper.setUp()
        self.addCleanup(helper.temporary.cleanup)
        original = helper.root / 'IMGVR5_UViG.tsv.gz'
        original.write_bytes(helper.old.read_bytes())
        portal = helper.root / 'METAVR_main_table.parquet'
        with adapter.connect(helper.work) as connection:
            connection.execute(f'COPY (SELECT * FROM {p.raw_reader(helper.new)}) TO {p.sql(portal)} (FORMAT PARQUET)')
        with tempfile.TemporaryDirectory(dir=ROOT.parent) as temporary, \
             patch.object(runner, 'FROZEN_SHA', p.digest(original)), \
             patch.object(runner, 'PORTAL_SHA', p.digest(portal)), \
             patch.object(runner, 'PORTAL_BYTES', portal.stat().st_size), \
             patch.object(runner, 'ROWS', 4), patch.object(runner, 'CHANGES', 3), \
             patch.object(runner, 'REVIEWED', dict(original_votu_groups=2, portal_votu_groups=3,
                        split_original_groups=2, merged_portal_groups=1)), \
             patch.object(runner.shutil, 'disk_usage', return_value=type('Disk', (), {'free': 100 * 1024**3})()):
            output = Path(temporary) / 'results'
            status = runner.run(argparse.Namespace(original_metadata=original, portal_metadata=portal, output_dir=output))
            self.assertEqual(status, 0)
            report = json.loads((output / 'report.json').read_text())
            self.assertEqual(report['after']['status'], 'EXACT_MATCH')
            self.assertTrue(all(report['inputs_unchanged'].values()))
            self.assertFalse(report['network_requests_performed'])
            self.assertFalse(report['database_build_performed'])

    def test_runner_wrong_pin_no_patch_and_report(self):
        import run_metavr_patch_validation as runner
        with tempfile.TemporaryDirectory(dir=ROOT.parent) as temporary:
            directory = Path(temporary)
            source = directory / 'source'
            source.mkdir()
            original, portal = source / 'IMGVR5_UViG.tsv.gz', source / 'METAVR_main_table.parquet'
            original.write_bytes(b'not matching metadata')
            portal.write_bytes(b'tiny')
            output = directory / 'results'
            with patch.object(runner, 'PORTAL_BYTES', 4), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(runner.run(argparse.Namespace(original_metadata=original, portal_metadata=portal, output_dir=output)), 1)
            report = json.loads((output / 'report.json').read_text())
            self.assertIn('checksum differs', report['error'])
            self.assertFalse((output / 'patch_manifest.json').exists())

    def test_review_excludes_all_large_metadata_and_environment_files(self):
        import package_metavr_patch_review as reviews
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'results/private_metadata').mkdir(parents=True)
            (root / 'results/report.json').write_text('{"status":"PASS"}')
            (root / 'results/votu_differences.tsv').write_text('small-export')
            (root / 'results/private_metadata/secret.tsv.gz').write_bytes(b'private-data')
            output = root / 'review.tar.gz'
            reviews.package_review(root, output)
            with tarfile.open(output) as archive:
                self.assertEqual(set(archive.getnames()), {'results/report.json', 'results/votu_differences.tsv', 'REVIEW_MANIFEST.json'})


if __name__ == '__main__':
    unittest.main(verbosity=2, buffer=True)
