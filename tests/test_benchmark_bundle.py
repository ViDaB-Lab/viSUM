"""Evaluation contracts only; no classifier or upstream-tool execution."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from benchmark_bundle import Bundle, sha
import evaluate_prepared_benchmark as prepared_eval


def write_table(path, values, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(values[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(values)


def fixture(root, sample='S', nested=True, originals=('a','b','c','d'), calls=None, legacy=False):
    folder = Path(root) / (sample + '_results') / 'viharmony'
    folder.mkdir(parents=True)
    disp = [dict(normalized_name=sample+'__c'+str(i),original_contig_name=name) for i,name in enumerate(originals)]
    calls = calls or [('retained_viral','viral_candidates','high',0),
                      ('provisional_provirus','viral_candidates','supported',1),
                      ('TE_candidate','TE_candidates','supported',2),
                      ('likely_nonviral','excluded','provisional',3)]
    meta=[]
    for i,(decision,destination,confidence,parent) in enumerate(calls):
        row=dict(final_sequence_id=sample+'__record'+str(i),normalized_name=disp[parent]['normalized_name'],
                 original_contig_name=originals[parent],viral_decision=decision,viral_confidence=confidence,
                 candidate_destination=destination)
        if legacy: row.pop('candidate_destination')
        meta.append(row)
    primary=[r for r in meta if r['viral_decision'].startswith('retained_')]
    provisional=[r for r in meta if r['viral_decision'].startswith('provisional_')]
    outputs={}
    for label,suffix,values in [('disposition','sequence_disposition.tsv',disp),('metadata','final_metadata.tsv',meta),
                                ('database_tsv','database_candidates.tsv',primary),
                                ('provisional_metadata','provisional_metadata.tsv',provisional)]:
        relative = ('audit/legacy/' if label in ('database_tsv','provisional_metadata') else 'audit/') if nested else ''
        name=relative+sample+'.'+suffix
        fields=list(disp[0]) if label=='disposition' else list(meta[0])
        write_table(folder/name,values,fields)
        outputs[label]=dict(file=name,sha256=sha(folder/name))
    manifest=dict(sample_id=sample,outputs=outputs,policy=dict(rna_pair_homology_floor=True))
    if nested:manifest['reporting_bundle']=dict(output_root='..')
    path=folder/('audit/' if nested else '')/(sample+'.harmonizer_manifest.json')
    path.write_text(json.dumps(manifest))
    return folder


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def test_nested_endpoints_supported_confidence_is_not_primary(self):
        bundle=Bundle(fixture(self.root),'S');mapping,sets=bundle.endpoints()
        self.assertEqual(sets['primary'],{'S__c0'})
        self.assertEqual(sets['viral_candidates'],{'S__c0','S__c1'})
        self.assertEqual(sets['TE_candidates'],{'S__c2'})
        self.assertEqual(sets['TE_only'],{'S__c2'})
        self.assertEqual(len(mapping),4)

    def test_flat_and_nested_have_identical_sets(self):
        old=Bundle(fixture(self.root/'flat',nested=False),'S')
        new=Bundle(fixture(self.root/'nested'),'S')
        self.assertEqual(old.endpoints(),new.endpoints())

    def test_historical_without_destination_column(self):
        b=Bundle(fixture(self.root,nested=False,legacy=True),'S')
        self.assertEqual(b.endpoints()[1]['primary'],{'S__c0'})
        self.assertEqual(b.endpoints()[1]['viral_candidates'],{'S__c0','S__c1'})
        self.assertFalse(b.endpoints()[1]['TE_candidates'])

    def test_regions_count_once_and_mixed_parent_is_not_TE_only(self):
        b=Bundle(fixture(self.root,calls=[('retained_provirus','viral_candidates','supported',0),
            ('retained_provirus','viral_candidates','supported',0),('TE_candidate','TE_candidates','high',0)]),'S')
        sets=b.endpoints()[1]
        self.assertEqual(sets['primary'],{'S__c0'})
        self.assertEqual(sets['TE_candidates'],{'S__c0'})
        self.assertFalse(sets['TE_only'])

    def test_ambiguous_supported_is_not_a_viral_export(self):
        b=Bundle(fixture(self.root,calls=[('ambiguous_review','review_only','supported',0)]),'S')
        self.assertFalse(b.endpoints()[1]['viral_candidates'])

    def test_corrupt_artifact_fails(self):
        path=fixture(self.root)
        (path/'audit/S.final_metadata.tsv').write_text('changed')
        with self.assertRaisesRegex(ValueError,'checksum'):Bundle(path,'S')

    def test_two_manifest_layouts_fail(self):
        path=fixture(self.root)
        (path/'S.harmonizer_manifest.json').write_bytes((path/'audit/S.harmonizer_manifest.json').read_bytes())
        with self.assertRaisesRegex(ValueError,'exactly one'):Bundle(path,'S')

    def test_missing_artifact_requires_explicit_archived_mode(self):
        path=fixture(self.root);(path/'audit/legacy/S.provisional_metadata.tsv').unlink()
        with self.assertRaises(FileNotFoundError):Bundle(path,'S')
        b=Bundle(path,'S',allow_missing=True)
        self.assertTrue(any(r['status']=='not_extracted_or_not_packaged' for r in b.integrity))
        with self.assertRaises(FileNotFoundError):b.table('provisional_metadata')

    def test_manifest_path_traversal_fails(self):
        path=fixture(self.root);m=path/'audit/S.harmonizer_manifest.json'
        value=json.loads(m.read_text());value['outputs']['metadata']['file']='../S.final_metadata.tsv'
        m.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError,'Unsafe'):Bundle(path,'S')


class PreparedEvaluationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.prepared=self.root/'prepared';self.prepared.mkdir()
        self.results=self.root/'results';self.output=self.root/'evaluation'
        truth=[]
        for name,cls,route in [('a','virus','dna'),('b','chromosome','dna'),('c','virus','unresolved'),('d','virus','rna')]:
            truth.append(dict(sequence_id=name,parent_sequence_id='parent_'+name,source_class=cls,input_type=route,
                              sequence_sha256='f'*64,length_nt='4'))
        write_table(self.prepared/'truth.routed.tsv',truth)
        execution=[]
        for batch,route,ids in [('DNA','dna','abc'),('RNA','rna','cd')]:
            for name in ids:
                r=next(r for r in truth if r['sequence_id']==name)
                execution.append(dict(benchmark_batch=batch,sequence_id=name,input_type=route,
                    evaluation_stratum='unresolved_'+route if name=='c' else 'source_routed',
                    parent_sequence_id=r['parent_sequence_id'],source_class=r['source_class'],sequence_sha256='f'*64,length_nt='4'))
        write_table(self.prepared/'execution_map.tsv',execution)
        (self.prepared/'manifest.csv').write_text('prefix,type,fasta\nDNA,dna,DNA.fa\nRNA,rna,RNA.fa\n')
        (self.prepared/'preparation_summary.json').write_text(json.dumps(dict(status='INPUTS_PREPARED_NOT_LAUNCHED',
            unresolved_policy='both',unique_input_sequences=4,planned_execution_sequences=5,batches=dict(DNA=3,RNA=2))))
        self.rehash()
        fixture(self.results,'DNA',originals=tuple('abc'),calls=[('retained_viral','viral_candidates','high',0),
            ('TE_candidate','TE_candidates','supported',1),('provisional_viral','viral_candidates','provisional',2)])
        fixture(self.results,'RNA',originals=tuple('cd'),calls=[('likely_nonviral','excluded','provisional',0),
            ('retained_viral','viral_candidates','supported',1)])

    def rehash(self):
        names=('preparation_summary.json','truth.routed.tsv','execution_map.tsv','manifest.csv')
        (self.prepared/'checksums.sha256').write_text(''.join(sha(self.prepared/name)+'  '+name+'\n' for name in names))

    def test_strata_separate_truth_unchanged_and_TE_not_viral(self):
        before=sha(self.prepared/'truth.routed.tsv')
        result=prepared_eval.evaluate(self.prepared,self.results,self.output)
        self.assertEqual(result['unique_inputs'],4);self.assertEqual(result['execution_sequences'],5)
        self.assertFalse(result['paired_arms_pooled']);self.assertEqual(sha(self.prepared/'truth.routed.tsv'),before)
        with (self.output/'metrics.tsv').open() as f:metrics=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual({r['evaluation_stratum'] for r in metrics},{'source_routed','unresolved_dna','unresolved_rna'})
        row=next(r for r in metrics if r['evaluation_stratum']=='unresolved_dna' and r['endpoint']=='primary')
        self.assertEqual(row['called'],'0')
        row=next(r for r in metrics if r['evaluation_stratum']=='unresolved_dna' and r['endpoint']=='viral_candidates')
        self.assertEqual(row['called'],'1')
        row=next(r for r in metrics if r['source_class']=='chromosome' and r['endpoint']=='TE_only')
        self.assertEqual(row['called'],'1');self.assertEqual(row['interpretation'],'TE_export_rate_not_viral_detection')

    def test_changed_truth_without_updated_checksum_fails(self):
        with (self.prepared/'truth.routed.tsv').open('a') as f:f.write('changed\n')
        with self.assertRaisesRegex(ValueError,'checksum'):prepared_eval.evaluate(self.prepared,self.results,self.output)

    def test_missing_result_batch_fails(self):
        with self.assertRaisesRegex(ValueError,'sample directories'):
            prepared_eval.evaluate(self.prepared,self.root/'missing',self.output)

    def test_wrong_paired_stratum_fails_even_if_rehashed(self):
        p=self.prepared/'execution_map.tsv';p.write_text(p.read_text().replace('unresolved_rna','source_routed'))
        self.rehash()
        with self.assertRaisesRegex(ValueError,'frozen truth/routing'):
            prepared_eval.evaluate(self.prepared,self.results,self.output)

    def test_missing_paired_arm_fails(self):
        p=self.prepared/'execution_map.tsv';lines=p.read_text().splitlines()
        p.write_text('\n'.join(line for line in lines if 'unresolved_rna' not in line)+'\n');self.rehash()
        with self.assertRaisesRegex(ValueError,'paired execution coverage'):
            prepared_eval.evaluate(self.prepared,self.results,self.output)

    def test_existing_output_refused(self):
        self.output.mkdir()
        with self.assertRaises(FileExistsError):prepared_eval.evaluate(self.prepared,self.results,self.output)


if __name__=='__main__':unittest.main()
