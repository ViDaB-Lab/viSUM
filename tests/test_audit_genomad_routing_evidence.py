import hashlib
import gzip
import io
import importlib.util
import json
from pathlib import Path
import tempfile
import tarfile
import unittest

spec = importlib.util.spec_from_file_location('routing_audit',Path(__file__).resolve().parents[1]/'bin/audit_genomad_routing_evidence.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class RoutingEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.parents = {'AA1.1':dict(length_nt='20',canonical_sequence_sha256='abc')}
        self.rna = dict(matches=[],source_url='https://example.org/rna',source_sha256='sourcehash')

    def record(self,accession='AA1.1',length=20,molecule='dna'):
        path = Path(self.tmp.name)/'cache.json'
        path.write_text(json.dumps({'result':{'uids':['123'],'123':dict(accessionversion=accession,slen=length,moltype=molecule)}}))
        return path

    def test_exact_accession_and_length(self):
        path = self.record()
        evidence,rejected = mod.collect_evidence(self.parents,self.rna,[path])
        self.assertEqual(evidence['AA1.1'][0]['route'],'dna')
        self.assertIn(hashlib.sha256(path.read_bytes()).hexdigest(),evidence['AA1.1'][0]['record'])
        self.assertFalse(rejected)

    def test_version_mismatch_not_joined(self):
        evidence,_ = mod.collect_evidence(self.parents,self.rna,[self.record(accession='AA1.2')])
        self.assertFalse(evidence)

    def test_length_mismatch_rejected(self):
        evidence,rejected = mod.collect_evidence(self.parents,self.rna,[self.record(length=19)])
        self.assertFalse(evidence)
        self.assertEqual(rejected[0]['reason'],'NCBI_length_mismatch')

    def test_unknown_molecule_not_guessed(self):
        evidence,rejected = mod.collect_evidence(self.parents,self.rna,[self.record(molecule='unknown')])
        self.assertFalse(evidence)
        self.assertEqual(len(rejected),1)

    def test_conflicting_sources_preserved(self):
        self.rna['matches'] = [dict(parent='AA1.1',length=20,sha256='abc',source_id='ND_1')]
        evidence,_ = mod.collect_evidence(self.parents,self.rna,[self.record()])
        self.assertEqual({v['route'] for v in evidence['AA1.1']},{'dna','rna'})

    def test_false_rna_match_rejected(self):
        self.rna['matches'] = [dict(parent='AA1.1',length=20,sha256='bad',source_id='ND_1')]
        with self.assertRaises(ValueError): mod.collect_evidence(self.parents,self.rna,[])

    def test_missing_parent_not_defaulted(self):
        evidence,_ = mod.collect_evidence(self.parents,self.rna,[])
        self.assertFalse(evidence)

    def package(self, bad_checksum=False, duplicate=False):
        row = 'AA1.1\t2\t20\tabc\tabc\tAA1.1 original header\n'
        table = 'parent_sequence_id\tfragment_count\tlength_nt\tsequence_sha256\tcanonical_sequence_sha256\toriginal_header\n' + row
        if duplicate: table += row
        packed = gzip.compress(table.encode())
        summary = dict(status='PASS',required_parents=1,viral_fragments=2,
            fingerprints_sha256='bad' if bad_checksum else hashlib.sha256(packed).hexdigest(),truth_sha256='truthhash')
        path = Path(self.tmp.name)/'fingerprints.tar.gz'
        with tarfile.open(path,'w:gz') as archive:
            for name,data in [('summary.json',json.dumps(summary).encode()),('viral_parent_fingerprints.tsv.gz',packed),('missing_parents.txt',b'')]:
                info = tarfile.TarInfo(name)
                info.size = len(data)
                archive.addfile(info,io.BytesIO(data))
        return path

    def test_fingerprint_checksum_required(self):
        with self.assertRaisesRegex(ValueError,'checksum'): mod.load_parents(self.package(bad_checksum=True))

    def test_duplicate_parents_rejected(self):
        with self.assertRaisesRegex(ValueError,'Duplicate'): mod.load_parents(self.package(duplicate=True))

    def test_audit_never_authorizes_launch_and_keeps_conflict(self):
        package = self.package()
        self.rna['matches'] = [dict(parent='AA1.1',length=20,sha256='abc',source_id='ND_1')]
        rna_path = Path(self.tmp.name)/'rna.json'
        rna_path.write_text(json.dumps(self.rna))
        cache = Path(self.tmp.name)/'cache'
        cache.mkdir()
        self.record().rename(cache/'record.json')
        out = Path(self.tmp.name)/'audit'
        report = mod.audit(package,rna_path,cache,out)
        self.assertFalse(report['launch_ready'])
        self.assertEqual(report['fragment_counts'],{'unresolved':2})
        self.assertEqual(report['conflict_count'],1)
        self.assertEqual(len((out/'routing.candidate.tsv').read_text().splitlines()),1)
        self.assertIn('conflicting_source_annotations',(out/'unresolved.tsv').read_text())
        self.assertTrue((out/'checksums.sha256').is_file())
        with self.assertRaises(FileExistsError): mod.audit(package,rna_path,cache,out)


if __name__ == '__main__': unittest.main()
