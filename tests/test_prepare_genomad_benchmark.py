import argparse
import contextlib
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

BIN = Path(__file__).resolve().parents[1] / 'bin'
sys.path.insert(0, str(BIN))
import prepare_genomad_benchmark as prep
import build_genomad_routing as routing
import verify_genomad_preparation as verify


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'source'
        for name in [prep.TRUTH, prep.CLUSTERS, prep.SUMMARY, prep.EXCLUDED]:
            (self.source / name).parent.mkdir(parents=True, exist_ok=True)
        self.rows = []
        for parent, cls, cluster, split in [('C.1', 'chromosome', 'c', 1),
                                           ('P.1', 'plasmid', 'p', 2),
                                           ('V.1', 'virus', 'v', 3),
                                           ('R.1', 'virus', 'r', 4)]:
            self.rows.append(dict(zip(prep.REQUIRED, [parent + '|fragment_1', cluster, cls,
                                                      '1', '0', '1', '1', '1', '1', '4'])))
        (self.source / prep.CLUSTERS).write_text('c\t1\t1\tC.1\np\t1\t2\tP.1\nv\t1\t3\tV.1\nr\t1\t4\tR.1\n')
        (self.source / prep.EXCLUDED).write_text('EXCLUDED|fragment_1\n')
        (self.source / prep.SUMMARY).write_text(json.dumps(dict(retained_sequences=4,
            retained_by_source_class={'chromosome': 1, 'plasmid': 1, 'virus': 2},
            total_nucleotides=16, author_exclusions_removed=1)))
        self.write_truth()
        self.routes = self.root / 'routes.tsv'
        self.routes.write_text('parent_sequence_id\tinput_type\tevidence_source\tevidence_record\n'
                               'V.1\tdna\tfixture\tV.1 genomic DNA\nR.1\trna\tfixture\tR.1 genomic RNA\n')
        self.fasta = self.root / 'input.fna'
        self.fasta.write_text(''.join('>' + row['sequence_id'] + '\nACGT\n' for row in self.rows))

    def write_truth(self):
        with (self.source / prep.TRUTH).open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=prep.REQUIRED, delimiter='\t')
            w.writeheader()
            w.writerows(self.rows)

    def run_prep(self, command='build', name='out', routing_path=True, unresolved_policy='error'):
        args = argparse.Namespace(command=command, out=self.root / name, metadata=self.source,
                                  routing=self.routes if routing_path else None,
                                  fasta=self.fasta, seed='test', shards=2,unresolved_policy=unresolved_policy)
        with contextlib.redirect_stdout(io.StringIO()):
            return prep.prepare(args)

    def test_inspect_marks_unresolved(self):
        report = self.run_prep('inspect', routing_path=False)
        self.assertEqual(report['unresolved_fragments'], 2)
        self.assertFalse((self.root / 'out/manifest.csv').exists())

    def test_unresolved_build_fails(self):
        with self.assertRaisesRegex(ValueError, 'source-supported'):
            self.run_prep(routing_path=False)
        self.assertFalse((self.root / 'out/manifest.csv').exists())

    def test_build_preserves_truth_and_sequences(self):
        report = self.run_prep()
        self.assertEqual(sum(report['batches'].values()), 4)
        with (self.root / 'out/truth.routed.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual([{k: r[k] for k in prep.REQUIRED} for r in rows], self.rows)
        self.assertEqual([r['author_cluster_split'] for r in rows], ['1', '2', '3', '4'])
        self.assertEqual([r['input_type'] for r in rows], ['dna', 'dna', 'dna', 'rna'])

    def test_deterministic_sharding(self):
        self.run_prep(name='a')
        self.fasta.write_text(''.join('>' + r['sequence_id'] + '\nACGT\n' for r in reversed(self.rows)))
        self.run_prep(name='b')
        def assignments(name):
            with (self.root / name / 'truth.routed.tsv').open() as f:
                return {r['sequence_id']: r['benchmark_batch'] for r in csv.DictReader(f, delimiter='\t')}
        self.assertEqual(assignments('a'), assignments('b'))

    def test_refuses_existing_output(self):
        self.run_prep()
        with self.assertRaises(FileExistsError):
            self.run_prep()

    def test_duplicate_fasta(self):
        self.fasta.write_text(self.fasta.read_text() + '>V.1|fragment_1\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.run_prep()
        self.assertFalse((self.root / 'out/manifest.csv').exists())

    def test_missing_fasta(self):
        self.fasta.write_text('>V.1|fragment_1\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'missing from FASTA'):
            self.run_prep()

    def test_length_mismatch(self):
        self.fasta.write_text('>V.1|fragment_1\nACG\n')
        with self.assertRaisesRegex(ValueError, 'length/alphabet'):
            self.run_prep()

    def test_nan_weight(self):
        self.rows[0]['author_weight_split_0'] = 'nan'
        self.write_truth()
        with self.assertRaisesRegex(ValueError, 'Invalid weight'):
            self.run_prep()

    def test_cluster_mismatch(self):
        self.rows[0]['reference_cluster'] = 'v'
        self.write_truth()
        with self.assertRaisesRegex(ValueError, 'Conflicting parent'):
            self.run_prep()

    def test_blank_route_provenance(self):
        self.routes.write_text(self.routes.read_text().replace('fixture', ''))
        with self.assertRaisesRegex(ValueError, 'provenance'):
            self.run_prep()

    def test_fragment_parser_keeps_embedded_pipes(self):
        self.assertEqual(prep.parent_id('IMGVR|123|fragment_2'), 'IMGVR|123')

    def test_gzip_fasta(self):
        path = self.root / 'input.fna.gz'
        with gzip.open(path, 'wt') as f:
            f.write(self.fasta.read_text())
        self.fasta = path
        self.assertEqual(self.run_prep()['status'], 'INPUTS_PREPARED_NOT_LAUNCHED')

    def test_excluded_fragment_rejected(self):
        (self.source / prep.EXCLUDED).write_text(self.rows[0]['sequence_id'] + '\n')
        with self.assertRaisesRegex(ValueError, 'excluded fragment'):
            self.run_prep()

    def test_duplicate_truth_rejected(self):
        self.rows.append(dict(self.rows[0]))
        self.write_truth()
        with self.assertRaises(prep.sqlite3.IntegrityError):
            self.run_prep()

    def test_unexpected_fasta_rejected(self):
        self.fasta.write_text(self.fasta.read_text() + '>UNKNOWN|fragment_1\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'Unexpected'):
            self.run_prep()

    def test_extra_routes_rejected(self):
        self.routes.write_text(self.routes.read_text() + 'missing\tdna\tfixture\trecord\n')
        with self.assertRaisesRegex(ValueError, 'absent from viral truth'):
            self.run_prep()

    def test_summary_mismatch_rejected(self):
        path = self.source / prep.SUMMARY
        data = json.loads(path.read_text())
        data['total_nucleotides'] = 17
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'disagree'):
            self.run_prep()

    def test_exact_xml_mapping_only(self):
        required = self.root / 'required.tsv'
        required.write_text('parent_sequence_id\nV.1\nR.1\nIMGVR_unknown\nV.2\n')
        xml = self.root / 'records.xml'
        xml.write_text('<GBSet><GBSeq><GBSeq_accession-version>V.1</GBSeq_accession-version>'
                       '<GBSeq_primary-accession>V</GBSeq_primary-accession><GBSeq_moltype>DNA</GBSeq_moltype>'
                       '</GBSeq><GBSeq><GBSeq_accession-version>R.1</GBSeq_accession-version>'
                       '<GBSeq_moltype>ss-RNA</GBSeq_moltype></GBSeq></GBSet>')
        with contextlib.redirect_stdout(io.StringIO()):
            report = routing.resolve(required, [xml], self.root / 'mapped', 'fixture')
        self.assertEqual(report['resolved_parents'], 2)
        self.assertEqual(report['unresolved_parents'], 2)
        self.assertIsNone(routing.molecule_route('DNA/RNA'))

    def test_dual_routing_keeps_unique_truth(self):
        report=self.run_prep(routing_path=False,unresolved_policy='both')
        self.assertEqual(report['unique_input_sequences'],4)
        self.assertEqual(report['planned_execution_sequences'],6)
        with (self.root/'out/truth.routed.tsv').open() as handle:
            truth=list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(len(truth),4)
        with (self.root/'out/execution_map.tsv').open() as handle:
            executions=list(csv.DictReader(handle,delimiter='\t'))
        for row in self.rows:
            matches=[r for r in executions if r['sequence_id']==row['sequence_id']]
            self.assertEqual(len(matches),2 if row['source_class']=='virus' else 1)
            if len(matches)==2:
                self.assertEqual({r['evaluation_stratum'] for r in matches},{'unresolved_dna','unresolved_rna'})
                self.assertEqual(len({r['sequence_sha256'] for r in matches}),1)
                self.assertEqual(len({r['benchmark_batch'].rsplit('_',1)[1] for r in matches}),1)
        with contextlib.redirect_stdout(io.StringIO()): result=verify.verify(self.root/'out')
        self.assertEqual(result['execution_sequences'],6)

    def test_dual_policy_preserves_known_routes(self):
        report=self.run_prep(unresolved_policy='both')
        self.assertEqual(report['planned_execution_sequences'],4)
        with contextlib.redirect_stdout(io.StringIO()): result=verify.verify(self.root/'out')
        self.assertEqual(result['strata'],{'source_routed':4})

    def test_partial_known_and_dual_routes(self):
        self.routes.write_text('parent_sequence_id\tinput_type\tevidence_source\tevidence_record\nR.1\trna\tfixture\tRNA genome\n')
        report=self.run_prep(unresolved_policy='both')
        self.assertEqual(report['planned_execution_sequences'],5)
        with contextlib.redirect_stdout(io.StringIO()): result=verify.verify(self.root/'out')
        self.assertEqual(result['strata'],{'source_routed':3,'unresolved_dna':1,'unresolved_rna':1})

    def test_verifier_rejects_modified_fasta(self):
        self.run_prep(unresolved_policy='both')
        path=next((self.root/'out').glob('*.fasta'))
        path.write_text(path.read_text().replace('ACGT','AAAA'))
        with self.assertRaisesRegex(ValueError,'Checksum mismatch'): verify.verify(self.root/'out')

    def test_dual_build_still_rejects_duplicate_source(self):
        self.fasta.write_text(self.fasta.read_text()+'>V.1|fragment_1\nACGT\n')
        with self.assertRaisesRegex(ValueError,'duplicate'):
            self.run_prep(routing_path=False,unresolved_policy='both')
        self.assertFalse((self.root/'out/manifest.csv').exists())

    def rehash_outputs(self):
        root=self.root/'out'
        (root/'checksums.sha256').write_text(''.join(f'{prep.sha256(path)}  {path.name}\n'
            for path in sorted(root.iterdir()) if path.name!='checksums.sha256'))

    def test_verifier_checks_pairing_beyond_checksums(self):
        self.run_prep(routing_path=False,unresolved_policy='both')
        path=self.root/'out/execution_map.tsv'
        lines=path.read_text().splitlines()
        path.write_text('\n'.join(line for line in lines if 'unresolved_rna' not in line)+'\n')
        self.rehash_outputs()
        with self.assertRaisesRegex(ValueError,'Missing planned executions'):
            verify.verify(self.root/'out')

    def test_failed_dual_build_never_exports_manifest(self):
        self.fasta.write_text('>C.1|fragment_1\nACGT\n')
        with self.assertRaisesRegex(ValueError,'missing from FASTA'):
            self.run_prep(routing_path=False,unresolved_policy='both')
        self.assertFalse((self.root/'out/manifest.csv').exists())


if __name__ == '__main__':
    unittest.main()
