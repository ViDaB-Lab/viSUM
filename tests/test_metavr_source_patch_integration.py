"""Offline source-preparation wiring tests; synthetic sources, never live data."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'bin'), str(ROOT / 'tests')]
import prepare_vicat_metavr_sources as sources
import patch_vicat_metavr_metadata as patches
import test_vicat_metavr_sources as fixtures


class IntegratedPatchTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.MetavrTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.temporary.cleanup)
        f = self.fixture
        changed = f.root / 'changed.parquet'
        with sources.connect(f.root) as connection:
            connection.execute(f'''COPY (SELECT * REPLACE (
                CASE WHEN uvig='U1' THEN 'portal-new-group' ELSE votu END AS votu)
                FROM read_parquet({sources.sql(f.new_metadata)})) TO {sources.sql(changed)} (FORMAT PARQUET)''')
        f.new_metadata = changed
        raw = f.root / 'raw.tsv.gz'
        sources.convert_metadata(changed, raw, 2, f.root)
        f.manifest['files']['metadata'].update(bytes=changed.stat().st_size, sha256=sources.digest(changed))
        self.bundle = f.root / 'bundle'
        patches.export_patch(f.old_metadata, raw, self.bundle, f.root, 2, 1,
            dict(original_metadata_sha256=sources.digest(f.old_metadata), portal_sha256=sources.digest(changed)), sources.connect)
        f.manifest['metadata_compatibility_patch'] = dict(manifest='bundle/patch_manifest.json',
            patch='bundle/votu_compatibility_patch.tsv', manifest_sha256=sources.digest(self.bundle / 'patch_manifest.json'))
        f.save_template()

    def test_metadata_path_uses_patch_without_protein_access_or_network(self):
        f = self.fixture
        before = f.old_metadata.read_bytes(), f.new_metadata.read_bytes()
        destination = f.root / 'prepared.tsv.gz'
        with (patch.object(sources, 'protein_fingerprint', side_effect=AssertionError('No proteins')),
             patch.object(sources, 'fetch', side_effect=AssertionError('No fetch')),
             patch.object(sources, 'urlopen', side_effect=AssertionError('No network'))):
            receipt = sources.prepare_metadata(f.template, f.manifest, f.new_metadata, destination, f.root / 'work')
        self.assertEqual(receipt['applied_rows'], 1)
        self.assertEqual(receipt['non_votu_changed_rows'], 0)
        self.assertEqual(sources.taxonomy_fingerprint(destination, f.root), sources.taxonomy_fingerprint(f.old_metadata, f.root))
        self.assertEqual(before, (f.old_metadata.read_bytes(), f.new_metadata.read_bytes()))

    def test_compare_prepare_and_offline_reuse_share_patched_metadata(self):
        f = self.fixture
        result = f.compare()
        self.assertEqual(result['comparison']['metadata_compatibility_patch']['applied_rows'], 1)
        provenance, calls = f.prepare()
        self.assertEqual(calls, 4)  # HEAD/GET on each tiny in-memory fixture; no network.
        self.assertEqual(provenance['metadata_compatibility_patch']['applied_rows'], 1)
        prepared = Path(provenance['metadata'])
        self.assertEqual(sources.taxonomy_fingerprint(prepared, f.root), sources.taxonomy_fingerprint(f.old_metadata, f.root))
        before = prepared.read_bytes(), prepared.stat().st_mtime_ns
        with (patch.object(sources, 'urlopen', side_effect=AssertionError('No network')),
             patch.object(sources, 'convert_metadata', side_effect=AssertionError('No reconversion')),
             patch.object(sources, 'prepare_metadata', side_effect=AssertionError('No repatching'))):
            reused = sources.prepare_sources(f.verified, f.root / 'cache', f.manifest['frozen'])
        self.assertEqual(provenance, reused)
        self.assertEqual(before, (prepared.read_bytes(), prepared.stat().st_mtime_ns))

    def test_native_prepared_checksum_approval_avoids_repeat_protein_scan(self):
        f = self.fixture
        result = f.compare()
        expected = f.root / 'expected.tsv.gz'
        sources.prepare_metadata(f.template, f.manifest, f.new_metadata, expected, f.root / 'expected-work')
        result['comparison']['metadata'] = dict(records=2, prepared_metadata_sha256=sources.digest(expected))
        sources.atomic_json(f.verified, result)
        with patch.object(sources,'protein_fingerprint',side_effect=AssertionError('No repeat protein scan')):
            provenance, calls = f.prepare()
        self.assertEqual(provenance['metadata_sha256'], sources.digest(expected))
        receipt = f.root / 'cache/download_provenance.json'
        before = receipt.read_bytes(), receipt.stat().st_mtime_ns
        with (patch.object(sources,'urlopen',side_effect=AssertionError('No network')),
              patch.object(sources,'prepare_metadata',side_effect=AssertionError('No repatching'))):
            sources.prepare_sources(f.verified,f.root/'cache',f.manifest['frozen'])
        self.assertEqual(before,(receipt.read_bytes(),receipt.stat().st_mtime_ns))

    def test_recipe_tampering_and_path_escape_rejected_before_fetch(self):
        f = self.fixture
        for field, value in (('manifest_sha256', '0' * 64), ('patch', '../escaped.tsv')):
            with self.subTest(field=field):
                altered = copy.deepcopy(f.manifest)
                altered['metadata_compatibility_patch'][field] = value
                sources.atomic_json(f.template, altered)
                with patch.object(sources, 'urlopen', side_effect=AssertionError('No network')):
                    with self.assertRaises(ValueError):
                        sources.load_manifest(f.template, verified=False)
        f.save_template()
        (self.bundle / 'votu_compatibility_patch.tsv').write_bytes(b'changed-patch')
        with self.assertRaisesRegex(ValueError, 'patch checksum mismatch'):
            sources.load_manifest(f.template, verified=False)

    def test_changed_portal_source_rejected_before_conversion(self):
        f = self.fixture
        f.new_metadata.write_bytes(f.new_metadata.read_bytes() + b'changed')
        with patch.object(sources, 'convert_metadata', side_effect=AssertionError('No conversion')):
            with self.assertRaisesRegex(ValueError, 'recipe source'):
                sources.prepare_metadata(f.template, f.manifest, f.new_metadata, f.root / 'out.tsv.gz', f.root / 'work')

    def test_wrong_recipe_frozen_identity_rejected(self):
        f = self.fixture
        f.manifest['frozen']['metadata_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'pinned source identities'):
            sources.prepare_metadata(f.template, f.manifest, f.new_metadata, f.root / 'out.tsv.gz', f.root / 'work')

    def test_default_native_bundle_and_approved_protein_identity(self):
        path = ROOT / 'assets/vicat/metavr_portal_v1.json'
        manifest = sources.load_manifest(path, verified=False)
        bundle = sources.patch_bundle(path, manifest)
        self.assertEqual(bundle['identity']['expected_changes'], 1409)
        with (bundle['paths']['patch']).open(encoding='utf-8') as stream:
            self.assertEqual(sum(1 for _ in stream), 1410)
        self.assertEqual(manifest['files']['proteins']['sha256'],
            '30ac6f74eceae3490ddd26e169d63efe66c1ae410775890ab029653afad49a8d')
        self.assertEqual(sources.load_manifest(path)['compatibility'], 'verified-full-source-equivalence')


if __name__ == '__main__':
    unittest.main(verbosity=2)
