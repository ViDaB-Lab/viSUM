"""Presentation contracts; no classification or sequence inference."""
from copy import deepcopy
import csv
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import reporting_bundle as r


def metadata(**updates):
    result = dict.fromkeys(r.SUMMARY_COLUMNS, 'NA')
    result.update(final_sequence_id='s__c1', normalized_name='s__c1', original_contig_name='original',
                  record_type='input_contig', provirus_coordinates='', refined_length=100,
                  viral_decision='provisional_viral', candidate_destination='viral_candidates',
                  strict_taxonomy='d__Viruses;' + ';'.join(p+'__unclassified' for p in ('r','k','p','c','o','f','g','s')),
                  strict_taxonomy_rank='domain', strict_taxonomy_confidence='unclassified',
                  strict_taxonomy_supporting_tools='NA')
    result.update(updates)
    return result


class ReportingTests(unittest.TestCase):
    def test_exact_schema_sizes(self):
        self.assertEqual(len(r.SUMMARY_COLUMNS),26)
        self.assertEqual(len(r.TAXONOMY_COLUMNS),10)
        self.assertEqual(len(r.QUALITY_COLUMNS),24)
        self.assertEqual(len(r.GROUP_COLUMNS),14)
        self.assertEqual(r.REVIEW_COLUMNS, r.SUMMARY_COLUMNS+['review_category','review_reason'])

    def test_missing_not_zero(self):
        self.assertEqual(r.cell(0),'0')
        self.assertEqual(r.cell(''),'NA')
        self.assertEqual(r.cell(None),'NA')
        self.assertEqual(r.cell('Not-determined'),'Not-determined')

    def test_complete_category_not_viral_identity(self):
        m=metadata(viral_decision='TE_candidate',candidate_destination='TE_candidates')
        evidence=[dict(tool='checkv',sequence_id='s__c1',checkv_quality='Complete',contamination='0.0')]
        q=list(r.quality_rows('s',m,evidence))
        self.assertEqual(q[0]['checkv_complete_call'],'true')
        result=r.summary_row('s',m,[],q)
        self.assertEqual(result['candidate_status'],'TE_candidate')
        self.assertEqual(result['checkv_complete_call'],'true')

    def test_determined_noncomplete_is_not_proven_incomplete(self):
        for quality in ('High-quality','Medium-quality','Low-quality'):
            self.assertEqual(r.complete_call(quality,'exact_candidate'),'false')

    def test_unknown_and_absent_are_distinct(self):
        for quality in ('Not-determined','','NA'):
            self.assertEqual(r.complete_call(quality,'exact_candidate'),'NA')
        absent=list(r.quality_rows('s',metadata(),[]))[0]
        self.assertEqual(absent['assessment_status'],'no_matching_record')
        available=list(r.quality_rows('s',metadata(),[dict(tool='checkv',sequence_id='s__c1',checkv_quality='Not-determined')]))[0]
        self.assertEqual(available['assessment_status'],'available')
        self.assertEqual(available['checkv_quality'],'Not-determined')

    def test_parent_complete_is_never_child_complete(self):
        m=metadata(final_sequence_id='child',record_type='provirus',provirus_coordinates='10-40')
        q=list(r.quality_rows('s',m,[dict(tool='checkv',sequence_id='s__c1',checkv_quality='Complete')]))
        self.assertEqual(q[0]['assessment_scope'],'parent_context')
        self.assertEqual(q[0]['checkv_complete_call'],'NA')
        self.assertEqual(r.summary_row('s',m,[],q)['checkv_quality'],'NA')

    def test_exact_proviral_ID_does_not_hide_aggregate_quality(self):
        m=metadata(final_sequence_id='child',record_type='provirus')
        q=list(r.quality_rows('s',m,[dict(tool='checkv',sequence_id='child',parent_sequence_id='s__c1',
                                        record_type='provirus',coordinates='10-40',checkv_quality='Complete')]))
        self.assertEqual(q[0]['assessment_scope'],'proviral_aggregate_context')
        self.assertEqual(q[0]['checkv_complete_call'],'NA')
        self.assertEqual(q[0]['assessment_coordinates'],'NA')

    def test_disagreeing_assessments_not_arbitrarily_selected(self):
        q=[dict(assessment_scope='exact_candidate',checkv_quality='Complete',checkv_complete_call='true'),
           dict(assessment_scope='exact_candidate',checkv_quality='Low-quality',checkv_complete_call='false')]
        s=r.summary_row('s',metadata(),[],q)
        self.assertEqual(s['checkv_quality'],'NA')
        self.assertEqual(s['checkv_complete_call'],'NA')

    def test_parent_and_coordinates(self):
        self.assertEqual(r.summary_row('s',metadata(),[],[])['normalized_parent_id'],'NA')
        m=metadata(final_sequence_id='child',record_type='provirus',provirus_coordinates='10-40')
        row=r.summary_row('s',m,[],[])
        self.assertEqual(row['normalized_parent_id'],'s__c1')
        self.assertEqual(row['coordinates'],'10-40')
        self.assertEqual(row['original_sequence_id'],'original')

    def test_conflict_types_and_review_flags_preserved(self):
        m=metadata(cellular_conflict_tools='deep6,vicat',plasmid_conflict_tools='vicat')
        row=r.summary_row('s',m,['cellular_conflict','taxonomy_conflict'],[])
        self.assertEqual(row['conflict_tools'],'cellular:deep6;cellular:vicat;plasmid:vicat')
        self.assertEqual(row['review_flags'],'cellular_conflict,taxonomy_conflict')

    def test_taxonomy_is_projection_not_inference(self):
        m=metadata()
        row=r.taxonomy_row(m)
        self.assertEqual(row['domain'],'Viruses')
        self.assertEqual(row['family'],'NA')
        self.assertEqual(row['genus'],'NA')
        with self.assertRaises(ValueError):r.taxonomy_row(dict(m,strict_taxonomy='Viruses'))

    def test_database_domains_and_zero_support_kept_separate(self):
        sources=[dict(database_domain=d,database_version='232',candidate_assignment_status='clustered',
                      family_prediction='novel_family_0',family_evidence='novel',family_network_support='0.0')
                 for d in ('prokaryotes','eukaryotes')]
        rows=list(r.group_rows('s','s__c1',sources))
        self.assertEqual(len(rows),2)
        self.assertEqual({z['database_domain'] for z in rows},{'prokaryotes','eukaryotes'})
        self.assertTrue(all(z['prediction_type']=='novel_group' and z['network_support']=='0.0' for z in rows))

    def test_reference_named_prediction_not_replaced_by_group(self):
        rows=list(r.group_rows('s','s__c1',[dict(family_prediction='Baculoviridae',family_evidence='reference',
                  genus_prediction='novel_genus_1',genus_evidence='novel')]))
        self.assertEqual([z['prediction_type'] for z in rows],['named_prediction','novel_group'])

    def test_unplaced_and_unresolved_are_not_formal_taxa(self):
        rows=list(r.group_rows('s','s__c1',[dict(realm_prediction='Riboviria|Duplodnaviria',
                                              phylum_prediction='unplaced_phylum_1')]))
        self.assertEqual([z['prediction_type'] for z in rows],['unresolved','unplaced_group'])

    def test_singleton_not_an_assigned_realm(self):
        rows=list(r.group_rows('s','s__c1',[dict(candidate_assignment_status='singleton',realm_prediction='singleton')]))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['rank'],'NA')
        self.assertEqual(rows[0]['raw_prediction'],'NA')
        self.assertEqual(rows[0]['assignment_status'],'singleton')

    def test_no_group_record_not_singleton(self):
        self.assertEqual(list(r.group_rows('s','s__c1',[]))[0]['assignment_status'],'no_matching_record')

    def test_inputs_not_mutated(self):
        m=metadata();q=[dict(tool='checkv',sequence_id='s__c1',checkv_quality='Complete')]
        g=[dict(family_prediction='novel_family_0')]
        original=deepcopy((m,q,g))
        r.summary_row('s',m,[],list(r.quality_rows('s',m,q)));list(r.group_rows('s','s__c1',g))
        self.assertEqual((m,q,g),original)

    def test_zero_candidates_still_produces_headers(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'s';outputs={}
            for label in ('viral_candidates_summary','TE_candidates_summary'):
                outputs[label]=Path(d)/(label+'.tsv')
            result=r.finalize(p,'s',[],[],{}, {},outputs)
            for label,columns in [('taxonomy',r.TAXONOMY_COLUMNS),('review_report',r.REVIEW_COLUMNS),
                                   ('checkv_report',r.QUALITY_COLUMNS),('vcontact3_report',r.GROUP_COLUMNS)]:
                with result[label].open() as f:self.assertEqual(next(csv.reader(f,delimiter='\t')),columns)
            self.assertFalse(any('review_candidates.fasta' in z.name for z in Path(d).rglob('*')))

    def test_review_views_overlap_without_removing_candidate(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'s';m=metadata();original=deepcopy(m)
            outputs={k:Path(d)/(k+'.tsv') for k in ('viral_candidates_summary','TE_candidates_summary')}
            review=[dict(final_sequence_id='s__c1',review_reasons='cellular_conflict,taxonomy_conflict')]
            result=r.finalize(p,'s',[m],review,{}, {},outputs)
            self.assertEqual(m,original)
            for key in ('review_reason_cellular_conflict','review_reason_taxonomy_conflict'):
                with result[key].open() as f:self.assertEqual(len(list(csv.DictReader(f,delimiter='\t'))),1)
            with result['viral_candidates_summary'].open() as f:self.assertEqual(len(list(csv.DictReader(f,delimiter='\t'))),1)

    def test_module_and_workflow_paths_agree(self):
        root=Path(__file__).resolve().parents[1]
        module=(root/'modules/local/viharmony.nf').read_text()
        workflow=(root/'visum_nextflow.nf').read_text()
        self.assertNotIn('reviewCandidatesFasta',workflow)
        self.assertNotIn('.review_candidates.fasta',module)
        self.assertIn("'reporting_bundle.py'",module)
        self.assertIn("MessageDigest.getInstance('SHA-256')",module)
        self.assertIn('ext harmonizer_source_digest:',module)
        self.assertIn('# viHARMONY source fingerprint: ${task.ext.harmonizer_source_digest}',module)
        self.assertNotIn('def sourceDigest =',module)
        for path in ('audit/legacy/${prefix}.final.normalized.fasta','audit/${prefix}.harmonizer_manifest.json',
                     'quality/${prefix}.checkv.tsv','groups/${prefix}.vcontact3.tsv','review/by_reason/'):
            self.assertIn(path,module)


if __name__=='__main__':unittest.main(verbosity=2)
