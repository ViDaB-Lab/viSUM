"""Small offline acquisition/build-wiring tests. No public downloads or searches."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import download_vicat_nonviral_sources as acquisition
import prepare_vicat_nonviral_build as preparation


class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.remote = self.root / 'remote'
        self.cache = self.root / 'cache'
        self.accession = 'GCF_000000001.1'
        self.manifest = self.root / 'manifest.tsv'
        self.manifest.write_text('assembly_accession\tcellular_group\torganism_name\ttaxid\n'
                                 + self.accession + '\tbacteria\tExample\t1\n')
        parent = self.remote / 'GCF/000/000/001'
        parent.mkdir(parents=True)
        self.basename = self.accession + '_ASM1v1'
        self.directory = parent / self.basename
        self.directory.mkdir()
        # file: URLs serve this offline directory-index fixture.
        self.listing = parent / 'index.html'
        self.listing.write_text(f'<a href="{self.basename}/">selected</a>'
                               '<a href="GCF_000000001.2_ASM1v2/">newer</a>')
        self.original_fetch = acquisition.fetch

        def fixture_fetch(url, destination, retries=3, timeout=120):
            if url.endswith('/'):
                url += 'index.html'
            return self.original_fetch(url, destination, retries, timeout)

        self.fetch_patch = patch.object(acquisition, 'fetch', side_effect=fixture_fetch)
        self.fetch = self.fetch_patch.start()
        self.addCleanup(self.fetch_patch.stop)
        self.addCleanup(patch.stopall)
        patch.object(acquisition, 'NCBI_ROOT', self.remote.as_uri() + '/').start()
        self.protein = self.directory / (self.basename + '_protein.faa.gz')
        self.feature = self.directory / (self.basename + '_feature_table.txt.gz')
        self.sequences = ''.join(f'>P{i}\nMKKLLAA\n' for i in range(7))
        self.protein.write_bytes(gzip.compress(self.sequences.encode()))
        types = ['chromosome', 'unplaced scaffold', 'plasmid', 'plastid', 'mitochondrion']
        table = '# feature\tseq_type\tgenomic_accession\tproduct_accession\n'
        table += ''.join(f'CDS\t{kind}\tNC_{i}\tP{i}\n' for i, kind in enumerate(types))
        table += 'CDS\tchromosome\tNC_C\tP5\nCDS\tplasmid\tNC_P\tP5\n'
        # P6 has no feature mapping: the existing unplaced fallback is preserved.
        self.feature.write_bytes(gzip.compress(table.encode()))
        self.write_checksums()

    def write_checksums(self):
        (self.directory / 'md5checksums.txt').write_text(''.join(
            hashlib.md5(path.read_bytes()).hexdigest() + '  ./' + path.name + '\n'
            for path in (self.protein, self.feature)))

    def download(self):
        return acquisition.download_sources(self.manifest, self.cache, workers=1, retries=1)

    def test_exact_version_download_layout_provenance_and_offline_reuse(self):
        proteins, features, provenance = self.download()
        self.assertEqual((proteins / 'ncbi_dataset/data' / self.accession / 'protein.faa').read_text(), self.sequences)
        self.assertTrue((features / 'bacteria' / self.accession / (self.accession + '.feature_table.txt.gz')).is_file())
        recorded = provenance.read_bytes()
        record = json.loads(recorded)['assemblies'][0]
        self.assertIn(self.basename, record['protein_url'])
        self.assertEqual(record['protein_sha256'], hashlib.sha256(self.sequences.encode()).hexdigest())
        with patch.object(acquisition, 'fetch', side_effect=AssertionError('Network on reuse')):
            self.download()
        self.assertEqual(provenance.read_bytes(), recorded)

    def test_missing_exact_version_never_uses_newer(self):
        self.listing.write_text('<a href="GCF_000000001.2_ASM1v2/">newer</a>')
        with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
            self.download()
        summary = json.loads((self.cache / 'download_summary.json').read_text())
        self.assertIn('No version substitution', summary['errors'][0]['error'])

    def test_ambiguous_directory_is_rejected(self):
        self.listing.write_text(f'<a href="{self.basename}/">a</a><a href="{self.accession}_OTHER/">b</a>')
        with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
            self.download()

    def test_wrong_checksum_stops_before_classification(self):
        self.protein.write_bytes(gzip.compress(b'>P0\nMUTATED\n'))
        with patch.object(preparation.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
                self.download()
            run.assert_not_called()
        self.assertIn('MD5 mismatch', (self.cache / 'download_summary.json').read_text())

    def test_missing_feature_table_checksum_fails(self):
        (self.directory / 'md5checksums.txt').write_text('')
        with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
            self.download()

    def test_corrupt_gzip_rejected_even_when_download_checksum_matches(self):
        self.feature.write_bytes(self.feature.read_bytes()[:-5])
        self.write_checksums()
        with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
            self.download()
        self.assertFalse((self.cache / 'downloads' / self.accession / 'receipt.json').exists())

    def test_invalid_fasta_rejected(self):
        self.protein.write_bytes(gzip.compress(b'not a FASTA\n'))
        self.write_checksums()
        with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
            self.download()

    def test_changed_completed_source_is_not_silently_replaced(self):
        proteins, _, _ = self.download()
        (proteins / 'ncbi_dataset/data' / self.accession / 'protein.faa').write_text('>x\nCHANGED\n')
        with patch.object(acquisition, 'fetch', side_effect=AssertionError('Unexpected repair')):
            with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
                self.download()
        self.assertIn('changed/missing', (self.cache / 'download_summary.json').read_text())

    def test_changed_manifest_rejected_before_downloads(self):
        self.download()
        self.manifest.write_text(self.manifest.read_text().replace('Example', 'Different'))
        with patch.object(acquisition, 'fetch', side_effect=AssertionError('Unexpected network')):
            with self.assertRaisesRegex(ValueError, 'manifest changed'):
                self.download()

    def test_partial_assembly_reuses_checksum_valid_protein_download(self):
        original = self.fetch.side_effect

        def fail_feature(url, destination, retries=3, timeout=120):
            if url.endswith('_feature_table.txt.gz'):
                raise OSError('interrupted')
            return original(url, destination, retries, timeout)

        self.fetch.side_effect = fail_feature
        with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
            self.download()
        self.fetch.side_effect = original
        self.fetch.reset_mock()
        self.download()
        self.assertFalse(any(call.args[0].endswith('_protein.faa.gz') for call in self.fetch.call_args_list))

    def test_downloaded_and_manual_inputs_have_identical_classification(self):
        proteins, features, provenance = self.download()
        args = argparse.Namespace(classified_dir=None, manifest=str(self.manifest),
            package_root=str(proteins), metadata_root=str(features), work_dir=str(self.root / 'manual'))
        manual = preparation.classified_inputs(args)
        args.work_dir = str(self.root / 'automatic')
        args.download_provenance = str(provenance)
        automatic = preparation.classified_inputs(args)
        for path in manual.iterdir():
            if path.name.endswith('.gz'):
                self.assertEqual(gzip.decompress(path.read_bytes()), gzip.decompress((automatic / path.name).read_bytes()), path.name)
            elif path.name.endswith('.tsv'):
                self.assertEqual(path.read_bytes(), (automatic / path.name).read_bytes(), path.name)
        evidence = json.loads((automatic / 'vicat_nonviral_source_provenance.json').read_text())
        self.assertEqual(evidence['download_sources']['assembly_count'], 1)
        self.assertEqual(json.loads((self.cache / 'download_summary.json').read_text())['status'], 'PASS')

    def test_prepare_cli_downloads_then_uses_unchanged_builder_arguments(self):
        arguments = ['prepare', '--manifest', str(self.manifest), '--download-sources',
            '--destination', str(self.root / 'db'), '--work-dir', str(self.root / 'build'), '--threads', '2']
        actual_run = subprocess.run
        builder_calls = []

        def run(command, **kwargs):
            if command[0] == 'bash':
                builder_calls.append(command)
                return subprocess.CompletedProcess(command, 0)
            return actual_run(command, **kwargs)

        with patch.object(sys, 'argv', arguments), patch.object(preparation.subprocess, 'run', side_effect=run):
            preparation.main()
        self.assertEqual(len(builder_calls), 1)
        self.assertNotIn('--min-seq-id', builder_calls[0])
        self.assertTrue((self.root / 'build/classified/vicat_nonviral_source_provenance.json').is_file())

    def test_download_failure_never_reaches_classifier_or_builder(self):
        self.listing.write_text('')
        arguments = ['prepare', '--manifest', str(self.manifest), '--download-sources',
            '--destination', str(self.root / 'db'), '--work-dir', str(self.root / 'build'), '--threads', '1']
        with patch.object(sys, 'argv', arguments), patch.object(preparation.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'required assemblies failed'):
                preparation.main()
            run.assert_not_called()

    def test_manifest_validation_and_frozen_inventory(self):
        frozen = ROOT / 'assets/vicat/nonviral_classaware_v1_20260904/assembly_manifest.tsv'
        self.assertEqual(len(acquisition.read_inventory(frozen)), 480)
        original = self.manifest.read_text()
        for invalid in (original.replace(self.accession, 'GCF_000000001'),
                        original.replace('bacteria', '../outside'), original + original.splitlines()[-1] + '\n'):
            self.manifest.write_text(invalid)
            with self.assertRaises(ValueError):
                acquisition.read_inventory(self.manifest)


class HttpIntegrityTests(unittest.TestCase):
    def test_incomplete_response_leaves_no_completed_file(self):
        class Response(io.BytesIO):
            headers = {'Content-Length': '100'}
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'file'
            with patch.object(acquisition, 'urlopen', return_value=Response(b'short')):
                with self.assertRaisesRegex(ValueError, 'Incomplete HTTP'):
                    acquisition.fetch('https://example.invalid/data', path, retries=1)
            self.assertFalse(path.exists())
            self.assertFalse(path.with_name('file.part').exists())

    def test_transient_error_is_retried(self):
        class Response(io.BytesIO):
            headers = {'Content-Length': '2'}
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'file'
            with patch.object(acquisition, 'urlopen', side_effect=[OSError('temporary'), Response(b'ok')]) as opened, patch.object(acquisition.time, 'sleep'):
                acquisition.fetch('https://example.invalid/data', path, retries=2)
            self.assertEqual(path.read_bytes(), b'ok')
            self.assertEqual(opened.call_count, 2)


class WorkflowWiringTests(unittest.TestCase):
    def test_setup_only_defaults_lock_and_helper_staging(self):
        setup = (ROOT / 'workflows/database_setup.nf').read_text()
        module = (ROOT / 'modules/local/vicat_nonviral_database.nf').read_text()
        self.assertIn('vicat_nonviral_auto_download = true', (ROOT / 'visum.config').read_text())
        self.assertIn('assets/vicat/nonviral_classaware_v1_20260904/assembly_manifest.tsv', setup)
        self.assertIn('download_vicat_nonviral_sources.py', setup)
        self.assertIn('download_vicat_refseq_metadata.py', setup)
        self.assertIn('hasLocalRoots && !nonviralRaw.every', setup)
        self.assertIn('Manifest-only nonviral setup requires', setup)
        self.assertLess(module.index('flock 9'), module.index('--download-sources'))
        self.assertIn('downloaded-and-built-from-assembly-manifest', module)
        self.assertNotIn('--download-sources', (ROOT / 'visum_nextflow.nf').read_text())


if __name__ == '__main__':
    unittest.main()
