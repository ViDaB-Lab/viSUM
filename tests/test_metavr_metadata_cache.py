"""Tiny metadata-only cache tests; no protein fixtures or network required."""
import contextlib
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import prepare_vicat_metavr_sources as sources
import patch_vicat_metavr_metadata as patches
from audit_vicat_metadata_equivalence import audit_metadata_equivalence


class MetadataCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.original = self.root / 'original.tsv.gz'
        self.portal = self.root / 'portal.parquet'
        with gzip.open(self.original, 'wt', encoding='utf-8') as stream:
            stream.write('uvig\tvotu\tictv_taxonomy\tictv_taxonomy_method\n'
                'U1\told\td__Viruses;f__Retroviridae\tmethod\n'
                'U2\told\td__Viruses;f__Retroviridae\tmethod\n')
        with sources.connect(self.root) as connection:
            connection.execute(f'''COPY (SELECT * REPLACE ('new' AS votu) FROM {patches.raw_reader(self.original)})
                TO {sources.sql(self.portal)} (FORMAT PARQUET)''')
        raw = self.root / 'raw.tsv.gz'
        sources.convert_metadata(self.portal, raw, 2, self.root)
        self.manifest_path = self.root / 'manifest.json'
        self.manifest = dict(schema='visum-metavr-sources-v1', compatibility='pending-full-source-comparison',
            frozen=dict(metadata_sha256=sources.digest(self.original), protein_sha256='0'*64, uvigs=2, proteins=1),
            files=dict(metadata=dict(name='METAVR_main_table.parquet', bytes=self.portal.stat().st_size,
                sha256=sources.digest(self.portal), url='https://example.invalid/metadata'),
                proteins=dict(name='METAVR.faa.bgz', bytes=1, sha256=None, url='https://example.invalid/proteins')))
        patches.export_patch(self.original, raw, self.root / 'patch', self.root, 2, 2,
            dict(original_metadata_sha256=sources.digest(self.original), portal_sha256=sources.digest(self.portal)), sources.connect)
        self.manifest['metadata_compatibility_patch'] = dict(manifest='patch/patch_manifest.json',
            patch='patch/votu_compatibility_patch.tsv', manifest_sha256=sources.digest(self.root / 'patch/patch_manifest.json'))
        sources.atomic_json(self.manifest_path, self.manifest)
        self.destination = self.root / 'cache/metadata.tsv.gz'

    def prepare(self, previous=None):
        with (patch.object(sources, 'protein_fingerprint', side_effect=AssertionError('No protein I/O')),
              patch.object(sources, 'fetch', side_effect=AssertionError('No source acquisition')),
              patch.object(sources, 'urlopen', side_effect=AssertionError('No network'))):
            return sources.prepare_metadata_cached(self.manifest_path, self.manifest, self.portal,
                                                   self.destination, self.root / 'cache', previous)

    def receipt(self, result):
        return dict(result, manifest_sha256=sources.digest(self.manifest_path))

    def test_prepare_and_unchanged_reuse_use_real_shared_helper(self):
        result = self.prepare()
        self.assertEqual(result['cache_status'], 'prepared')
        self.assertEqual(result['metadata_compatibility_patch']['applied_rows'], 2)
        self.assertEqual(audit_metadata_equivalence(self.original, self.destination, self.root, 2, sources.connect)['status'], 'EXACT_MATCH')
        before = self.destination.read_bytes(), self.destination.stat().st_mtime_ns
        with (patch.object(sources, 'prepare_metadata', side_effect=AssertionError('No repatching')),
              patch.object(sources, 'convert_metadata', side_effect=AssertionError('No conversion'))):
            reused = self.prepare(self.receipt(result))
        self.assertEqual(reused['cache_status'], 'reused')
        self.assertEqual(reused['metadata_compatibility_patch'], result['metadata_compatibility_patch'])
        self.assertEqual(before, (self.destination.read_bytes(), self.destination.stat().st_mtime_ns))

    def test_unreceipted_or_changed_metadata_rejected(self):
        result = self.prepare()
        with self.assertRaisesRegex(ValueError, 'changed/unreceipted'):
            self.prepare()
        self.destination.write_bytes(b'changed-private-fixture')
        with self.assertRaisesRegex(ValueError, 'changed/unreceipted'):
            self.prepare(self.receipt(result))

    def test_changed_manifest_or_missing_patch_provenance_rejected(self):
        result = self.prepare()
        wrong = dict(self.receipt(result), manifest_sha256='0'*64)
        with self.assertRaisesRegex(ValueError, 'manifest changed'):
            self.prepare(wrong)
        wrong = self.receipt(result)
        wrong.pop('metadata_compatibility_patch')
        with self.assertRaisesRegex(ValueError, 'patch provenance missing'):
            self.prepare(wrong)

    def test_pending_gate_remains_closed_before_protein_or_network(self):
        with (patch.object(sources, 'fetch', side_effect=AssertionError('No fetch')),
              patch.object(sources, 'protein_fingerprint', side_effect=AssertionError('No proteins'))):
            with self.assertRaisesRegex(ValueError, 'NOT yet verified'):
                sources.prepare_sources(self.manifest_path, self.root / 'blocked', self.manifest['frozen'])
        self.assertFalse((self.root / 'blocked').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2, buffer=True)
