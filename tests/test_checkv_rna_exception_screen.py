import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import test_checkv_rna_exception as screen

class CheckVExceptionScreenTests(unittest.TestCase):
    def setUp(self):
        self.row=dict(record_type='input_contig',viral_decision='ambiguous_review',
            cellular_conflict_tools='checkv',plasmid_conflict_tools='NA',
            qualified_tools='deep6,virbot,vicat',strong_tools='NA',
            sequence_interpretation='viral_candidate',original_length='2000')
        self.locus=dict(locus_classification='viral_supported',coordinates='51-1950',
            protein_length='632',best_bitscore='1100',best_identity='95',
            best_query_coverage='95',best_subject_coverage='95')

    def test_supported_candidate_passes_all_tiers_without_identifier(self):
        for tier in screen.TIERS.values():
            self.assertEqual(screen.eligible(self.row,[self.locus],tier),[self.locus])

    def test_conflicts_and_mobile_regions_cannot_bypass(self):
        for key,value in [('cellular_conflict_tools','checkv,vicat'),('plasmid_conflict_tools','genomad'),
                          ('record_type','provirus'),('sequence_interpretation','viral_retroelement_conflict')]:
            with self.subTest(key=key):
                self.assertEqual(screen.eligible(dict(self.row,**{key:value}),[self.locus],screen.TIERS['balanced']),[])

    def test_missing_additional_support_fails(self):
        self.assertEqual(screen.eligible(dict(self.row,qualified_tools='deep6,vicat'),[self.locus],screen.TIERS['balanced']),[])

    def test_cellular_or_ambiguous_locus_blocks_exception(self):
        for classification in ['cellular_supported','ambiguous']:
            self.assertEqual(screen.eligible(self.row,[self.locus,dict(self.locus,locus_classification=classification)],screen.TIERS['balanced']),[])

    def test_short_or_weak_matches_do_not_pass(self):
        for key,value in [('coordinates','1-100'),('protein_length','50'),('best_bitscore','50'),
                          ('best_identity','20'),('best_query_coverage','20'),('best_subject_coverage','20')]:
            with self.subTest(key=key):
                self.assertEqual(screen.eligible(self.row,[dict(self.locus,**{key:value})],screen.TIERS['balanced']),[])

if __name__=='__main__':unittest.main()
