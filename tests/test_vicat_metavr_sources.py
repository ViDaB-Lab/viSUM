"""Tiny offline tests. No public downloads, database installations or searches."""
import copy
import csv
import gzip
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
sys.path.insert(0, str(ROOT / 'tests'))
import prepare_vicat_metavr_sources as sources
import probe_metavr_portal as probe


class Response(io.BytesIO):
    def __init__(self, payload, status=200, headers=None):
        super().__init__(payload)
        self.status = status
        self.headers = {'Content-Length': str(len(payload)), 'ETag': 'stable'} if headers is None else headers

    def getcode(self):
        return self.status


class MetavrTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.old_proteins = self.root / 'frozen.faa.gz'
        self.new_proteins = self.root / 'METAVR.faa.bgz'
        # Different gzip bytes/descriptions/wrapping, same exact IDs/order/sequences.
        self.old_proteins.write_bytes(gzip.compress(b'>U1|p1 old\nMAA\nK\n>U2|p2\nMCD\n', mtime=1))
        self.new_proteins.write_bytes(gzip.compress(b'>U1|p1 portal\nMAAK\n>U2|p2 new\nMC\nD\n', mtime=2))
        self.old_metadata = self.root / 'frozen.tsv.gz'
        self.new_metadata = self.root / 'METAVR_main_table.parquet'
        with gzip.open(self.old_metadata, 'wt', newline='', encoding='utf-8') as target:
            writer = csv.writer(target, delimiter='\t', lineterminator='\n')
            writer.writerow(sources.REQUIRED_COLUMNS)
            writer.writerows([['U1', 'V1', 'd__Viruses;o__Ortervirales;f__Retroviridae', 'method1'],
                              ['U2', '\\N', 'd__Viruses;c__Naldaviricetes', '\\N']])
        with sources.connect(self.root) as connection:
            connection.execute('CREATE TABLE m (uvig VARCHAR, votu VARCHAR, ictv_taxonomy VARCHAR, ictv_taxonomy_method VARCHAR, irrelevant VARCHAR)')
            connection.execute("INSERT INTO m VALUES ('U2', NULL, 'd__Viruses;c__Naldaviricetes', NULL, 'unused'), ('U1', 'V1', 'd__Viruses;o__Ortervirales;f__Retroviridae', 'method1', 'unused')")
            connection.execute('COPY m TO ' + sources.sql(self.new_metadata) + ' (FORMAT PARQUET)')
        self.manifest = json.loads((ROOT / 'assets/vicat/metavr_portal_v1.json').read_text())
        # Legacy/unpatched tiny fixtures are separate from the packaged real patch.
        self.manifest.pop('metadata_compatibility_patch', None)
        self.manifest.pop('validation_evidence', None)
        self.manifest.pop('comparison', None)
        self.manifest['compatibility'] = 'pending-full-source-comparison'
        self.manifest['frozen'] = {'protein_sha256': sources.digest(self.old_proteins),
            'metadata_sha256': sources.digest(self.old_metadata), 'uvigs': 2, 'proteins': 2}
        self.template = self.root / 'template.json'
        self.save_template()
        self.verified = self.root / 'verified.json'

    def save_template(self):
        sources.atomic_json(self.template, self.manifest)

    def compare(self):
        return sources.compare_sources(self.template, self.old_proteins, self.old_metadata,
            self.new_proteins, self.new_metadata, self.verified, self.root / 'compare')

    def test_full_tiny_equivalence_approves_different_compression_and_metadata_order(self):
        before = {str(p): p.read_bytes() for p in (self.old_proteins, self.old_metadata, self.new_proteins, self.new_metadata)}
        result = self.compare()
        self.assertEqual(result['compatibility'], 'verified-full-source-equivalence')
        self.assertEqual(result['files']['proteins']['sha256'], sources.digest(self.new_proteins))
        self.assertNotEqual(result['files']['proteins']['sha256'], self.manifest['frozen']['protein_sha256'])
        self.assertEqual(result['comparison']['proteins']['records'], 2)
        self.assertFalse(result['comparison']['database_build_performed'])
        self.assertEqual(before, {str(p): p.read_bytes() for p in (self.old_proteins, self.old_metadata, self.new_proteins, self.new_metadata)})

    def test_changed_sequence_identifier_order_or_case_never_approves_manifest(self):
        for text in (b'>U1|p1\nMAAT\n>U2|p2\nMCD\n', b'>U1|wrong\nMAAK\n>U2|p2\nMCD\n',
                     b'>U2|p2\nMCD\n>U1|p1\nMAAK\n', b'>U1|p1\nMaak\n>U2|p2\nMCD\n'):
            with self.subTest(text=text):
                self.new_proteins.write_bytes(gzip.compress(text))
                with self.assertRaisesRegex(ValueError, 'protein IDs/sequences/order differ'):
                    self.compare()
                self.assertFalse(self.verified.exists())

    def test_changed_taxonomy_votu_and_method_rejected(self):
        for column in ('ictv_taxonomy', 'votu', 'ictv_taxonomy_method'):
            with self.subTest(column=column):
                changed = self.root / (column + '.parquet')
                with sources.connect(self.root) as connection:
                    connection.execute(f"COPY (SELECT * REPLACE ('changed' AS {column}) FROM read_parquet({sources.sql(self.new_metadata)})) TO {sources.sql(changed)} (FORMAT PARQUET)")
                with self.assertRaisesRegex(ValueError, 'taxonomy/vOTU/method metadata differ'):
                    sources.compare_sources(self.template, self.old_proteins, self.old_metadata, self.new_proteins, changed, self.verified, self.root / 'compare')
                self.assertFalse(self.verified.exists())

    def test_frozen_original_checksum_is_required(self):
        self.old_proteins.write_bytes(gzip.compress(b'>not-frozen\nMAA\n'))
        with self.assertRaisesRegex(ValueError, 'Frozen original source checksum mismatch'):
            self.compare()

    def test_incorrect_expected_protein_count_rejected(self):
        self.manifest['frozen']['proteins'] = 3
        self.save_template()
        with self.assertRaisesRegex(ValueError, 'protein IDs/sequences/order differ'):
            self.compare()

    def test_pending_manifest_fails_before_any_network_or_source_directory(self):
        with patch.object(sources, 'urlopen', side_effect=AssertionError('Network must not run')):
            with self.assertRaisesRegex(ValueError, 'NOT yet verified'):
                sources.prepare_sources(self.template, self.root / 'cache', self.manifest['frozen'])
        self.assertFalse((self.root / 'cache').exists())

    def test_manifest_frozen_identity_and_source_checksums_enforced(self):
        self.compare()
        wrong = dict(self.manifest['frozen'], proteins=3)
        with self.assertRaisesRegex(ValueError, 'configured frozen source'):
            sources.load_manifest(self.verified, wrong)
        value = json.loads(self.verified.read_text())
        value['files']['metadata']['sha256'] = None
        sources.atomic_json(self.verified, value)
        with self.assertRaisesRegex(ValueError, 'full SHA256'):
            sources.load_manifest(self.verified)

    def test_unsafe_source_names_and_non_https_urls_rejected(self):
        for field, value in (('name', '../escape'), ('url', 'http://example.invalid/data')):
            with self.subTest(field=field):
                changed = copy.deepcopy(self.manifest)
                changed['files']['proteins'][field] = value
                sources.atomic_json(self.template, changed)
                with self.assertRaisesRegex(ValueError, 'source name/HTTPS URL'):
                    sources.load_manifest(self.template, verified=False)

    def test_offline_probe_inspects_only_footer_and_prefix_not_full_data(self):
        import struct
        capture = self.root / 'capture'
        capture.mkdir()
        (capture / 'protein-prefix.bgz.part').write_bytes(self.new_proteins.read_bytes())
        parquet = self.new_metadata.read_bytes()
        length = struct.unpack('<I', parquet[-8:-4])[0]
        (capture / 'parquet-footer-only.bin').write_bytes(b'PAR1' + parquet[-8-length:])
        sources.atomic_json(capture / 'capture.json', {'file_headers': {}, 'expected_uvigs': 2, 'binary_response_bytes': length + 8 + self.new_proteins.stat().st_size})
        with patch.object(probe, 'urlopen', side_effect=AssertionError('Offline footer inspection must not use network')):
            result = probe.inspect_capture(capture, self.root / 'probe.json')
        self.assertEqual(result['metadata_rows'], 2)
        self.assertEqual(result['full_source_equivalence'], 'NOT VERIFIED')
        self.assertFalse(result['full_database_downloaded'])

    def test_conversion_preserves_required_values_and_nulls_for_actual_builder(self):
        target = self.root / 'converted.tsv.gz'
        sources.convert_metadata(self.new_metadata, target, 2, self.root)
        self.assertEqual(sources.taxonomy_fingerprint(self.old_metadata, self.root), sources.taxonomy_fingerprint(target, self.root))
        with gzip.open(target, 'rt', newline='', encoding='utf-8') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual({row['uvig']: row['ictv_taxonomy'] for row in rows}, {'U1': 'd__Viruses;o__Ortervirales;f__Retroviridae', 'U2': 'd__Viruses;c__Naldaviricetes'})
        self.assertEqual(rows[0]['votu'], '\\N')
        self.assertNotIn('irrelevant', rows[0])

    def test_actual_taxonomy_builder_outputs_match_for_frozen_and_converted_metadata(self):
        import build_vicat_taxonomy_lookup as builder
        converted = self.root / 'converted.tsv.gz'
        sources.convert_metadata(self.new_metadata, converted, 2, self.root)
        members = self.root / 'members.tsv.gz'
        with gzip.open(members, 'wt', encoding='utf-8') as stream:
            stream.write('U1|p1\tU1|p1\nU1|p1\tU2|p2\n')
        outputs = []
        for label, metadata in (('frozen', self.old_metadata), ('portal', converted)):
            output = self.root / (label + '-taxonomy')
            arguments = ['builder', '--cluster-members', str(members), '--metadata', str(metadata),
                '--output-dir', str(output), '--work-dir', str(self.root / (label + '-work')),
                '--threads', '1', '--memory-limit', '1GB', '--expected-member-count', '2',
                '--expected-representative-count', '1']
            with patch.object(sys, 'argv', arguments):
                builder.main()
            with sources.connect(self.root) as connection:
                outputs.append(connection.execute('SELECT * FROM read_parquet(?) ORDER BY 1',
                    [str(output / 'IMGVR5_UViG.vicat_taxonomy_lookup.parquet')]).fetchall())
        self.assertEqual(outputs[0], outputs[1])

    def test_missing_schema_duplicate_uvigs_and_wrong_count_rejected(self):
        for label, query, expected in (
            ('missing', "SELECT 'U1' AS uvig", 'missing required columns'),
            ('duplicate', f'SELECT * FROM read_parquet({sources.sql(self.new_metadata)}) UNION ALL SELECT * FROM read_parquet({sources.sql(self.new_metadata)})', 'row/unique-ID mismatch'),
            ('count', f'SELECT * FROM read_parquet({sources.sql(self.new_metadata)}) LIMIT 1', 'row/unique-ID mismatch')):
            path = self.root / (label + '.parquet')
            with sources.connect(self.root) as connection:
                connection.execute(f'COPY ({query}) TO {sources.sql(path)} (FORMAT PARQUET)')
            with self.assertRaisesRegex(ValueError, expected):
                sources.convert_metadata(path, self.root / (label + '.tsv.gz'), 2, self.root)

    def prepare(self):
        value = json.loads(self.verified.read_text())
        payloads = {value['files']['proteins']['url']: self.new_proteins.read_bytes(),
                    value['files']['metadata']['url']: self.new_metadata.read_bytes()}
        def open_fixture(request, **kwargs):
            return Response(payloads[request.full_url])
        with patch.object(sources, 'urlopen', side_effect=open_fixture) as opened:
            result = sources.prepare_sources(self.verified, self.root / 'cache', self.manifest['frozen'])
        return result, opened.call_count

    def test_automatic_download_and_offline_reuse(self):
        self.compare()
        result, calls = self.prepare()
        self.assertEqual(calls, 4)
        cache = self.root / 'cache'
        before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in cache.iterdir() if p.is_file() and p.name != 'download_provenance.json'}
        with patch.object(sources, 'urlopen', side_effect=AssertionError('Offline reuse must not use network')):
            reused = sources.prepare_sources(self.verified, cache, self.manifest['frozen'])
        self.assertEqual(result, reused)
        self.assertEqual(before, {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in cache.iterdir() if p.is_file() and p.name != 'download_provenance.json'})

    def test_changed_completed_or_converted_sources_are_not_silently_repaired(self):
        self.compare()
        self.prepare()
        converted = self.root / 'cache/METAVR.vicat_metadata.tsv.gz'
        converted.write_bytes(b'corrupt')
        with patch.object(sources, 'urlopen', side_effect=AssertionError('No repair network')):
            with self.assertRaisesRegex(ValueError, 'metadata changed'):
                sources.prepare_sources(self.verified, self.root / 'cache', self.manifest['frozen'])
        protein = self.root / 'cache/METAVR.faa.bgz'
        protein.write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'Completed MetaVR source changed'):
            sources.prepare_sources(self.verified, self.root / 'cache', self.manifest['frozen'])

    def test_changed_manifest_rejected_before_network(self):
        self.compare()
        self.prepare()
        value = json.loads(self.verified.read_text())
        value['additional_note'] = 'changed identity'
        sources.atomic_json(self.verified, value)
        with patch.object(sources, 'urlopen', side_effect=AssertionError('No network')):
            with self.assertRaisesRegex(ValueError, 'manifest changed'):
                sources.prepare_sources(self.verified, self.root / 'cache', self.manifest['frozen'])


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.path = Path(self.temporary.name) / 'source'
        self.payload = b'verified-source-bytes'
        self.entry = dict(bytes=len(self.payload), sha256=sources.hashlib.sha256(self.payload).hexdigest(),
                          url='https://example.invalid/source')

    def head(self):
        return Response(b'', headers={'Content-Length':str(len(self.payload)), 'ETag':'stable'})

    def seed(self, prefix):
        self.path.with_name('source.part').write_bytes(prefix)
        sources.atomic_json(self.path.with_name('source.http_identity.json'),
            dict(url=self.entry['url'], bytes=len(self.payload), etag='stable', last_modified=None))

    def server(self, request, **kwargs):
        if request.get_method() == 'HEAD':
            return self.head()
        start,end = map(int, request.get_header('Range').removeprefix('bytes=').split('-'))
        return Response(self.payload[start:end+1],206,{'Content-Length':str(end-start+1),
            'Content-Range':f'bytes {start}-{end}/{len(self.payload)}', 'ETag':'stable'})

    def test_valid_partial_resumes_with_range_and_if_range(self):
        self.seed(self.payload[:5])
        with patch.object(sources,'urlopen',side_effect=self.server) as opened:
            sources.fetch(self.entry,self.path,retries=1)
        self.assertEqual(opened.call_args.args[0].get_header('Range'),f'bytes=5-{len(self.payload)-1}')
        self.assertEqual(opened.call_args.args[0].get_header('If-range'),'stable')
        self.assertEqual(self.path.read_bytes(),self.payload)

    def test_server_ignoring_range_preserves_prefix(self):
        self.seed(self.payload[:5])
        with patch.object(sources,'urlopen',side_effect=[self.head(),Response(self.payload)]):
            with self.assertRaisesRegex(ValueError,'NOT truncated'):
                sources.fetch(self.entry,self.path,retries=1)
        self.assertEqual(self.path.with_name('source.part').read_bytes(),self.payload[:5])

    def test_exact_broken_pipe_retries_from_written_prefix(self):
        class Broken(Response):
            def read(inner,size=-1):
                if inner.tell()==0:
                    return super(Broken,inner).read(5)
                raise BrokenPipeError(32,'Broken pipe')
        calls=[]
        def server(request,**kwargs):
            calls.append(request)
            if request.get_method()=='HEAD':
                return self.head()
            if len(calls)==2:
                return Broken(self.payload)
            return self.server(request)
        with patch.object(sources,'urlopen',side_effect=server),patch.object(sources.time,'sleep'):
            sources.fetch(self.entry,self.path,retries=2)
        self.assertEqual(calls[-1].get_header('Range'),f'bytes=5-{len(self.payload)-1}')
        self.assertEqual(self.path.read_bytes(),self.payload)

    def test_consecutive_retry_budget_resets_after_bounded_chunks(self):
        class Broken(Response):
            def read(inner,size=-1):
                raise BrokenPipeError(32,'Broken pipe')
        failed=set()
        def server(request,**kwargs):
            response=self.server(request)
            if request.get_method()=='HEAD':
                return response
            start=int(request.get_header('Range').split('=')[1].split('-')[0])
            if start not in failed:
                failed.add(start)
                return Broken(b'',response.status,response.headers)
            return response
        with patch.object(sources,'urlopen',side_effect=server),patch.object(sources.time,'sleep'):
            sources.fetch(self.entry,self.path,retries=2,chunk_bytes=4)
        self.assertGreater(len(failed),3)
        self.assertEqual(self.path.read_bytes(),self.payload)

    def test_complete_partial_is_checksum_verified_without_get(self):
        self.seed(self.payload)
        with patch.object(sources,'urlopen',return_value=self.head()) as opened:
            sources.fetch(self.entry,self.path,retries=1)
        self.assertEqual(opened.call_count,1)
        self.assertEqual(self.path.read_bytes(),self.payload)

    def test_wrong_range_length_checksum_and_oversized_responses_never_publish(self):
        cases=(('range',Response(self.payload,206,{'Content-Range':'bytes 1-20/21'})),
               ('length',Response(self.payload,headers={'Content-Length':'999'})),
               ('hash',Response(b'x'*len(self.payload))),
               ('oversized',Response(self.payload+b'extra',headers={'ETag':'stable'})),
               ('short',Response(b'short',headers={'ETag':'stable'})))
        for label,response in cases:
            with self.subTest(label=label),tempfile.TemporaryDirectory() as temporary:
                destination=Path(temporary)/'out'
                with patch.object(sources,'urlopen',side_effect=[self.head(),response]):
                    with self.assertRaises((OSError,ValueError)):
                        sources.fetch(self.entry,destination,retries=1)
                self.assertFalse(destination.exists())

    def test_unreceipted_partial_fails_before_network(self):
        self.path.with_name('source.part').write_bytes(self.payload[:5])
        with patch.object(sources,'urlopen',side_effect=AssertionError('No network')):
            with self.assertRaisesRegex(ValueError,'lacks HTTP identity'):
                sources.fetch(self.entry,self.path)

    def test_bounded_probe_refuses_unbounded_server_without_reading_body(self):
        response=Response(b'ignored')
        with patch.object(probe,'urlopen',return_value=response),patch.object(response,'read',side_effect=AssertionError('Must not read')):
            with self.assertRaisesRegex(ValueError,'refused to read'):
                probe.bounded_range('https://example.invalid/proteins',0,5,99)

class WiringTests(unittest.TestCase):
    def test_setup_only_defaults_lock_reuse_and_unverified_preflight(self):
        config = (ROOT / 'visum.config').read_text()
        workflow = (ROOT / 'workflows/database_setup.nf').read_text()
        module = (ROOT / 'modules/local/vicat_database.nf').read_text()
        self.assertIn('vicat_metavr_auto_download = true', config)
        self.assertIn('verified-full-source-equivalence', workflow)
        self.assertLess(workflow.index('portalIdentity.compatibility'), workflow.index('ready = Channel.value(true)'))
        self.assertLess(module.index('flock 9'), module.index('prepare_vicat_metavr_sources.py'))
        self.assertLess(module.index("emit_database 'skipped-existing-managed-database'"), module.index('prepare_vicat_metavr_sources.py'))
        self.assertIn('METAVR_SOURCE_DIR="\\${DB_DEST}.sources"', module)
        self.assertIn('overwrite: true', module)
        self.assertIn('failOnError: true', module)
        self.assertNotIn('prepare_vicat_metavr_sources.py', (ROOT / 'visum_nextflow.nf').read_text())

    def test_frozen_source_pins_remain_unchanged(self):
        config = (ROOT / 'visum.config').read_text()
        manifest = json.loads((ROOT / 'assets/vicat/metavr_portal_v1.json').read_text())
        for key in ('protein', 'metadata'):
            parameter = 'vicat_metavr_' + ('proteins' if key == 'protein' else 'metadata') + '_sha256'
            self.assertEqual(re.search(parameter + r"\s*=\s*'([^']+)'", config).group(1), manifest['frozen'][key + '_sha256'])
        self.assertEqual(manifest['compatibility'], 'verified-full-source-equivalence')
        self.assertEqual(manifest['files']['proteins']['sha256'],
            '30ac6f74eceae3490ddd26e169d63efe66c1ae410775890ab029653afad49a8d')
        self.assertEqual(manifest['files']['metadata']['sha256'],
            'b1260bd415f0d09503e6e8c36c9f773b0c4411f53aecec449d6ed7db8e87dfb7')


if __name__ == '__main__':
    unittest.main()
