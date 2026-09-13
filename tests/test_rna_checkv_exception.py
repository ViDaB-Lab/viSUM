import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from run_viharmony import rna_checkv_exception_loci

class RNAExceptionTests(unittest.TestCase):
    def setUp(self):
        self.locus=dict(sequence_id='p',locus_id='l1',locus_classification='viral_supported',coordinates='51-1950',
          protein_length='630',best_bitscore='1100',best_identity='95',best_query_coverage='95',best_subject_coverage='95')
        self.args=dict(input_type='rna',decision='ambiguous_review',record_type='input_contig',cellular=['checkv'],plasmid=[],
          supported_tools={'deep6','virbot','vicat'},interpretation='viral_candidate',loci=[self.locus],parent_id='p',length=2000)
    def test_supported(self):self.assertEqual(rna_checkv_exception_loci(**self.args),['l1'])
    def test_dna_regions_other_conflicts_blocked(self):
        for key,value in [('input_type','dna'),('record_type','provirus'),('cellular',['checkv','vicat']),
                          ('plasmid',['genomad']),('supported_tools',{'deep6','vicat'}),('interpretation','ambiguous_mobile_element')]:
            with self.subTest(key=key):self.assertEqual(rna_checkv_exception_loci(**dict(self.args,**{key:value})),[])
    def test_competing_locus_blocked(self):
        for state in ['cellular_supported','ambiguous']:
            loci=[self.locus,dict(self.locus,locus_classification=state)]
            self.assertEqual(rna_checkv_exception_loci(**dict(self.args,loci=loci)),[])
    def test_missing_or_weak_identity_does_not_pass(self):
        for identity in ['', 'NA','49']:
            self.assertEqual(rna_checkv_exception_loci(**dict(self.args,loci=[dict(self.locus,best_identity=identity)])),[])
    def test_out_of_bounds_fails(self):
        with self.assertRaises(ValueError):rna_checkv_exception_loci(**dict(self.args,loci=[dict(self.locus,coordinates='1-2001')]))
    def test_workflow_supplies_loci_without_pair_floor(self):
        workflow=(Path(__file__).resolve().parents[1]/'visum_nextflow.nf').read_text()
        self.assertIn("(runVicat ? ['vicat_loci'] : [])",workflow)
        self.assertNotIn("if( rnaPairHomologyFloor )",workflow)

if __name__=='__main__':unittest.main()
