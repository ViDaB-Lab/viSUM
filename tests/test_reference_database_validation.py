import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('reference_validation', ROOT / 'bin/validate_reference_database.py')
validation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)


class ReferenceValidationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def asset(self, relative, data=b'fixture\n'):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def db(self, relative, headers=True):
        self.asset(relative)
        self.asset(relative + '.index', b'0\t0\t8\n')
        self.asset(relative + '.dbtype', b'\x00\x00\x00\x00')
        if headers:
            self.db(relative + '_h', False)

    def genomad_fixture(self):
        self.asset('version.txt', b'1.9\n')
        self.asset('genomad_marker_metadata.tsv', b'id\tvalue\nmarker\t1\n')
        self.asset('nodes.dmp')
        self.asset('names.dmp')
        for name in ('genomad_db', 'genomad_mini_db', 'genomad_integrase_db'):
            self.db(name)

    def vc_fixture(self):
        domains = {}
        for domain in ('prokaryotes', 'eukaryotes'):
            block = {}
            for key in ('Genomes Report', 'gene2genome', 'proteins', 'nucleotides'):
                path = f'v232/{domain}/{key.replace(" ", "_")}'
                self.asset(path)
                block[key] = path
            block['identities'] = {}
            for identity in ('0.3', '0.4', '0.5', '0.6', '0.7'):
                prefix = f'v232/{domain}/{identity}'
                self.db(prefix + '/db')
                self.db(prefix + '/cluster', False)
                self.asset(prefix + '/cluster2genome')
                block['identities'][identity] = {key: prefix + '/' + key for key in ('db', 'cluster', 'cluster2genome')}
            domains[domain] = block
        self.db('v232/vog/profiles')
        data = {'232': {'RefSeq': {'domains': domains}, 'VOGDB': {'mmseqs_db': 'v232/vog/profiles'}}}
        manifest = self.root / '232.json'
        manifest.write_text(json.dumps(data))
        return manifest, data

    def test_complete_genomad(self):
        self.genomad_fixture()
        validation.genomad(self.root)

    def test_old_minimal_genomad_check_would_accept_incomplete(self):
        self.asset('junk')
        self.asset('junk.dbtype')
        with self.assertRaises(ValueError):
            validation.genomad(self.root)

    def test_genomad_missing_header_index(self):
        self.genomad_fixture()
        (self.root / 'genomad_db_h.index').unlink()
        with self.assertRaisesRegex(ValueError, 'genomad_db_h.index'):
            validation.genomad(self.root)

    def test_genomad_empty_metadata_and_bad_version(self):
        self.genomad_fixture()
        self.asset('version.txt', b'nan')
        with self.assertRaises(ValueError):
            validation.genomad(self.root)
        self.asset('version.txt', b'1.9')
        self.asset('genomad_marker_metadata.tsv', b'header\tonly\n')
        with self.assertRaises(ValueError):
            validation.genomad(self.root)

    def test_split_mmseqs_data(self):
        self.db('split', False)
        (self.root / 'split').unlink()
        self.asset('split.0')
        self.asset('split.1')
        validation.mmseqs(self.root / 'split')
        (self.root / 'split.0').unlink()
        with self.assertRaises(ValueError):
            validation.mmseqs(self.root / 'split')

    def test_complete_vcontact3(self):
        manifest, _ = self.vc_fixture()
        validation.vcontact3(manifest)

    def test_old_minimal_vcontact3_check_would_accept_incomplete(self):
        self.asset('v232/junk')
        manifest = self.asset('232.json', b'{}')
        with self.assertRaises(KeyError):
            validation.vcontact3(manifest)

    def test_missing_domain_file_and_selected_domain(self):
        manifest, _ = self.vc_fixture()
        (self.root / 'v232/eukaryotes/proteins').unlink()
        with self.assertRaises(ValueError):
            validation.vcontact3(manifest)
        validation.vcontact3(manifest, 'prokaryotes')

    def test_missing_identity_and_vog(self):
        manifest, data = self.vc_fixture()
        del data['232']['RefSeq']['domains']['prokaryotes']['identities']['0.5']
        manifest.write_text(json.dumps(data))
        with self.assertRaises(KeyError):
            validation.vcontact3(manifest)
        manifest, _ = self.vc_fixture()
        (self.root / 'v232/vog/profiles.dbtype').unlink()
        with self.assertRaises(ValueError):
            validation.vcontact3(manifest)

    def test_manifest_wrong_version_and_path_escape(self):
        manifest, data = self.vc_fixture()
        data['233'] = data.pop('232')
        manifest.write_text(json.dumps(data))
        with self.assertRaises(KeyError):
            validation.vcontact3(manifest)
        data['232'] = data.pop('233')
        data['232']['RefSeq']['domains']['prokaryotes']['proteins'] = '../outside.faa'
        manifest.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'Unsafe'):
            validation.vcontact3(manifest)


if __name__ == '__main__':
    unittest.main()
