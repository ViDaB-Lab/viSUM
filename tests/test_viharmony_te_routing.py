"""TE routing integration: frozen expectations and complete candidate evaluation."""
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import run_viharmony as harmony
import te_routing as cf
import te_routing as cr
corrected_route = cf.corrected_route


def audit(family='', order='Ortervirales'):
    values = [*cf.ANCESTRY.values(), family, '', '']
    values[5] = order
    prefixes = ('d__','r__','k__','p__','c__','o__','f__','g__','s__')
    return {'uncapped_strict_taxonomy': ';'.join(p+(v or 'unclassified') for p,v in zip(prefixes,values,strict=True))}


def te(superfamily='Gypsy', strength='strong', complete='yes', category='retroelement', **extra):
    return {'tesorter_order':'LTR','tesorter_superfamily':superfamily,
            'tesorter_clade':'unknown','evidence_strength':strength,
            'tesorter_complete':complete,'assignment_method':'direct_hmm',
            'tesorter_evidence_category':category,**extra}


def decision(value='retained_viral'):
    return {'viral_decision':value}


class RoutingTests(unittest.TestCase):
    def call(self, rows, taxonomy=None, ordinary='retained_viral', votes=()):
        return cf.route(decision('ambiguous_review'), decision(ordinary), rows,
                        taxonomy or audit(), votes)

    def test_Gypsy_with_broad_Ortervirales_is_TE(self):
        self.assertEqual(self.call([te()])['destination'],'TE_candidates')

    def test_Gypsy_with_Metaviridae_is_compatible(self):
        self.assertEqual(self.call([te()],audit('Metaviridae'))['destination'],'TE_candidates')

    def test_Copia_with_Pseudoviridae_is_compatible(self):
        self.assertEqual(self.call([te('Copia')],audit('Pseudoviridae'))['destination'],'TE_candidates')

    def test_BelPao_with_Belpaoviridae_is_compatible(self):
        self.assertEqual(self.call([te('Bel-Pao')],audit('Belpaoviridae'))['destination'],'TE_candidates')

    def test_strong_Gypsy_with_accepted_Retroviridae_is_flagged_viral(self):
        result=self.call([te()],audit('Retroviridae'))
        self.assertEqual((result['destination'],result['flag']),('viral_candidates','TE_viral_conflict'))

    def test_competing_family_does_not_rescue_independent_exclusion(self):
        self.assertEqual(self.call([te()],audit('Retroviridae'),ordinary='likely_nonviral')['destination'],'excluded')

    def test_partial_TE_does_not_veto_ordinary_eligibility(self):
        result=self.call([te(strength='qualified',complete='no')],ordinary='provisional_viral')
        self.assertEqual((result['destination'],result['flag']),('viral_candidates','partial_TE_support'))

    def test_partial_TE_does_not_create_eligibility(self):
        self.assertEqual(self.call([te(strength='qualified',complete='no')],ordinary='ambiguous_review')['destination'],'review_only')

    def test_mixed_domains_preserve_eligibility_with_warning(self):
        result=self.call([te('mixture',strength='qualified',complete='unknown',category='ambiguous_mobile_element')],audit('Retroviridae'))
        self.assertEqual((result['destination'],result['flag']),('viral_candidates','mixed_mobile_element_evidence'))

    def test_second_pass_is_not_strong_TE(self):
        result=self.call([te('unknown',strength='weak',complete='none',assignment_method='second_pass')])
        self.assertEqual((result['destination'],result['flag']),('viral_candidates','weak_TE_signal'))

    def test_explicit_Retrovirus_is_compatible_not_proven_exogenous(self):
        result=self.call([te('Retrovirus',strength='qualified',complete='unknown',category='viral_like_mobile_element')],audit('Retroviridae'))
        self.assertEqual(result['flag'],'retrovirus_compatible')
        self.assertIn('endogenous_status_unresolved',result['reason'])

    def test_no_TE_hit_does_not_create_eligibility(self):
        self.assertEqual(self.call([],ordinary='ambiguous_review')['destination'],'review_only')

    def test_nonselected_prediction_warns_but_does_not_block_TE(self):
        votes=[dict(tool='vicat',rank='family',taxon='Retroviridae',tier='moderate',valid=True,accepted=False)]
        result=self.call([te()],audit('Metaviridae'),votes=votes)
        self.assertEqual(result['destination'],'TE_candidates')
        self.assertIn('vicat:family:Retroviridae',result['nonselected_predictions'])

    def test_mixed_clades_within_Copia_still_support_TE_superfamily(self):
        self.assertEqual(self.call([te('Copia',tesorter_clade='mixture')],audit('Pseudoviridae'))['destination'],'TE_candidates')

    def test_windowed_hit_does_not_label_whole_contig_as_TE(self):
        self.assertNotEqual(self.call([te(tesorter_window_was_split='true')])['destination'],'TE_candidates')

    def test_mixed_TE_superfamilies_not_automatically_TE_only(self):
        self.assertNotEqual(self.call([te(),te('Copia')])['destination'],'TE_candidates')

    def test_router_does_not_change_any_input(self):
        args=(decision(),decision(),[te()],audit(),[])
        before=deepcopy(args)
        cf.route(*args)
        self.assertEqual(args,before)

    def test_reported_rank_cap_not_used_for_routing(self):
        source=audit('Retroviridae')
        source['reported_strict_taxonomy']=audit()['uncapped_strict_taxonomy']
        first=self.call([te()],source)
        source['reported_strict_taxonomy']=audit('Retroviridae')['uncapped_strict_taxonomy']
        self.assertEqual(self.call([te()],source),first)

    def test_unresolved_family_conflict_is_not_absence_of_family_evidence(self):
        unresolved=audit()
        unresolved.update(uncapped_taxonomy_conflict='true',
            uncapped_strict_taxonomy_conflict_rank='family',
            uncapped_strict_taxonomy_conflicting_taxa='Metaviridae,Retroviridae')
        result=corrected_route(decision('ambiguous_review'),decision(),[te()],unresolved,[])
        self.assertEqual((result['destination'],result['flag']),('viral_candidates','TE_viral_conflict'))

    def test_genus_conflict_does_not_invalidate_TE_superfamily(self):
        unresolved=audit('Metaviridae')
        unresolved.update(uncapped_taxonomy_conflict='true',
            uncapped_strict_taxonomy_conflict_rank='genus',
            uncapped_strict_taxonomy_conflicting_taxa='Metavirus,Errantivirus')
        result=corrected_route(decision('ambiguous_review'),decision(),[te()],unresolved,[])
        self.assertEqual(result['destination'],'TE_candidates')

    def test_consensus_conflict_does_not_bypass_ordinary_exclusion(self):
        unresolved=audit()
        unresolved.update(uncapped_taxonomy_conflict='true',
            uncapped_strict_taxonomy_conflict_rank='family',
            uncapped_strict_taxonomy_conflicting_taxa='Metaviridae,Retroviridae')
        result=corrected_route(decision('ambiguous_review'),decision('likely_nonviral'),[te()],unresolved,[])
        self.assertEqual(result['destination'],'excluded')


def evidence(tool='tesorter', **extra):
    return dict(tool=tool, evidence_application='used_for_final_decision',
                classification='retroelement' if tool == 'tesorter' else 'virus',
                evidence_strength='qualified', assignment_method='direct_hmm',
                tesorter_order='LTR', tesorter_superfamily='Gypsy',
                tesorter_evidence_category='retroelement', tesorter_complete='no', **extra)


class TEReportingTests(unittest.TestCase):
    def label(self, row, previous='partial_TE_support'):
        return cr.clarify_te(previous, [row])['proposed_flag']

    def test_incomplete_architecture_is_not_weak_TE(self):
        self.assertEqual(self.label(evidence()), 'incomplete_TE_architecture')

    def test_unknown_architecture_is_not_asserted_incomplete(self):
        row = evidence()
        row.update(tesorter_complete='unknown', tesorter_order='LINE', tesorter_superfamily='unknown')
        self.assertEqual(self.label(row), 'limited_TE_domain_support')

    def test_explicit_mixture_has_mixed_label(self):
        row = evidence()
        row.update(tesorter_superfamily='mixture', tesorter_evidence_category='ambiguous_mobile_element')
        self.assertEqual(self.label(row), 'mixed_domain_assignments')

    def test_Maverick_is_ambiguous_not_mixed(self):
        row = evidence()
        row.update(tesorter_order='Maverick', tesorter_superfamily='unknown', domains='ATPase|Maverick',
                   tesorter_evidence_category='ambiguous_mobile_element')
        self.assertEqual(self.label(row, 'mixed_mobile_element_evidence'), 'ambiguous_mobile_element')

    def test_Polinton_is_ambiguous_not_mixed(self):
        row = evidence()
        row.update(tesorter_order='Polinton', tesorter_superfamily='unknown',
                   tesorter_evidence_category='ambiguous_mobile_element')
        self.assertEqual(self.label(row), 'ambiguous_mobile_element')

    def test_different_clades_do_not_make_mixed_superfamily(self):
        row = evidence()
        row['tesorter_clade'] = 'mixture'
        self.assertEqual(self.label(row), 'incomplete_TE_architecture')

    def test_different_superfamilies_across_rows_are_mixed(self):
        a, b = evidence(), evidence()
        b['tesorter_superfamily'] = 'Copia'
        self.assertEqual(cr.clarify_te('partial_TE_support', [a, b])['proposed_flag'], 'mixed_domain_assignments')

    def test_split_window_reports_localization(self):
        row = evidence()
        row['tesorter_window_was_split'] = 'true'
        self.assertEqual(self.label(row), 'localized_TE_domain_support')

    def test_substantive_guarded_flags_preserved(self):
        for label in ('TE_supported', 'TE_viral_conflict', 'retrovirus_compatible', 'weak_TE_signal'):
            with self.subTest(label=label):
                self.assertEqual(self.label(evidence(), label), label)

    def test_no_TE_rows(self):
        self.assertEqual(cr.clarify_te('NA', [])['proposed_flag'], 'NA')

    def test_parent_context_does_not_supply_TE_flag(self):
        row = evidence()
        row['evidence_application'] = 'parent_context_only'
        self.assertEqual(cr.clarify_te('NA', [row])['proposed_flag'], 'NA')


class ScopeReportingTests(unittest.TestCase):
    def scope(self, rows, **record):
        return cr.viral_scope(dict(record_type='input_contig', **record), rows)

    def vicat(self, pattern='localized_viral_cluster', viral='3', **extra):
        row = evidence('vicat')
        row.update(origin_pattern=pattern, viral_supported_loci=viral, cellular_supported_loci='21',
                   ambiguous_loci='3', uninformative_loci='79', orf_loci='106',
                   viral_cluster_coordinates='97037-98328', applied_scope='region', **extra)
        return row

    def test_localized_cluster_does_not_describe_whole_scaffold(self):
        result = self.scope([self.vicat()])
        self.assertEqual(result['viral_evidence_scope'], 'localized_cluster_only')
        self.assertIn('entire_candidate_viral_origin_not_established', result['viral_scope_interpretation'])
        self.assertEqual(result['candidate_record_scope'], 'input_contig')
        self.assertEqual(result['vicat_viral_supported_loci'], 3)
        self.assertEqual(result['vicat_total_loci'], 106)

    def test_record_level_prediction_does_not_erase_localization(self):
        result = self.scope([self.vicat(), evidence('gianthunter')])
        self.assertEqual(result['viral_evidence_scope'], 'localized_cluster_with_record_level_predictions')
        self.assertEqual(result['record_level_viral_prediction_tools'], 'gianthunter')

    def test_single_locus_not_whole_genome_even_when_100_percent(self):
        row = self.vicat('predominantly_viral', '1')
        row['orf_loci'] = '1'
        result = self.scope([row])
        self.assertEqual(result['viral_evidence_scope'], 'single_locus_viral_support')
        self.assertEqual(result['vicat_total_loci'], 1)

    def test_multi_locus_is_not_complete_genome_assertion(self):
        result = self.scope([self.vicat('predominantly_viral')])
        self.assertEqual(result['viral_evidence_scope'], 'multi_locus_viral_support')
        self.assertIn('not_basewise_or_complete_genome_validation', result['viral_scope_interpretation'])

    def test_provirus_child_coordinates_remain_candidate_relative(self):
        row = self.vicat('resolved_provirus_with_vicat_support', coordinates='22172-34646')
        row['viral_cluster_coordinates'] = '7272-11725'
        result = cr.viral_scope({'record_type': 'viral_region', 'coordinates': '22172-34646'}, [row])
        self.assertEqual(result['candidate_record_scope'], 'viral_region')
        self.assertEqual(result['viral_cluster_coordinates'], '7272-11725')
        self.assertEqual(result['viral_cluster_coordinate_system'], 'candidate_record_1_based_inclusive')

    def test_zero_counts_and_missing_counts_differ(self):
        result = self.scope([evidence('deep6')])
        self.assertEqual(result['vicat_viral_supported_loci'], 'NA')
        row = self.vicat('predominantly_viral')
        row['cellular_supported_loci'] = '0'
        self.assertEqual(self.scope([row])['vicat_cellular_supported_loci'], 0)

    def test_no_qualified_rows_is_explicit(self):
        self.assertEqual(self.scope([])['viral_evidence_scope'], 'no_qualified_viral_support')

    def test_CheckV_or_ML_are_record_level_not_spatial(self):
        self.assertEqual(self.scope([evidence('checkv'), evidence('deepmicroclass2')])['viral_evidence_scope'],
                         'record_level_predictions_only')

    def test_superseded_parent_and_weak_annotations_not_counted(self):
        row, weak = self.vicat(), evidence('checkv')
        row['evidence_application'] = 'superseded_parent_vicat'
        weak['evidence_strength'] = 'weak'
        self.assertEqual(self.scope([row, weak])['viral_evidence_scope'], 'no_qualified_viral_support')

    def test_multiple_projected_summaries_not_summed(self):
        with self.assertRaises(ValueError):
            self.scope([self.vicat(), self.vicat()])

    def test_negative_count_fails(self):
        with self.assertRaises(ValueError):
            self.scope([self.vicat(viral='-1')])



def complete_example(te_row=None, cellular=False):
    catalog = harmony.make_msl_catalog([
        ('Riboviria','Pararnavirae','Artverviricota','Revtraviricetes','Ortervirales','Metaviridae','Metavirus','Metavirus one'),
        ('Riboviria','Pararnavirae','Artverviricota','Revtraviricetes','Ortervirales','Retroviridae','Gammaretrovirus','Gammaretrovirus one')])
    rows=[dict(sequence_id='c1', tool='deep6', classification='virus', evidence_strength='qualified')]
    if te_row is not None:
        rows.append(dict(te_row, sequence_id='c1', tool='tesorter', classification='retroelement'))
    if cellular:
        rows.append(dict(sequence_id='c1',tool='vicat',classification='cellular',evidence_strength='qualified'))
    task=harmony.CandidateTask('c1',2000,
        dict(record_type='input_contig',coordinates='',original_length='2000'),'c1',
        dict(original_id='original',original_header='original header'),'ambiguous',rows,[],())
    context=harmony.CandidateContext('sample','dna',1000,False,False,False,True)
    return task,context,catalog


class IntegratedCandidateTests(unittest.TestCase):
    def test_partial_TE_restores_only_provisional_single_tool_eligibility(self):
        task,context,catalog=complete_example(te(strength='qualified',complete='no'))
        result=harmony.evaluate_candidate(task,context,catalog)
        self.assertEqual(result.metadata['viral_decision'],'provisional_viral')
        self.assertEqual(result.routing_audit['pre_TE_routing_decision'],'ambiguous_review')
        self.assertEqual(result.metadata['TE_flag'],'incomplete_TE_architecture')
        self.assertEqual(result.metadata['viral_confidence'],'provisional')

    def test_partial_TE_cannot_overcome_cellular_origin_safeguard(self):
        task,context,catalog=complete_example(te(strength='qualified',complete='no'),cellular=True)
        result=harmony.evaluate_candidate(task,context,catalog)
        self.assertEqual(result.metadata['viral_decision'],'likely_nonviral')
        self.assertEqual(result.metadata['candidate_destination'],'excluded')

    def test_strong_compatible_TE_is_not_a_viral_FASTA_decision(self):
        task,context,catalog=complete_example(te())
        result=harmony.evaluate_candidate(task,context,catalog)
        self.assertEqual(result.metadata['viral_decision'],'TE_candidate')
        self.assertEqual(result.metadata['candidate_destination'],'TE_candidates')
        self.assertEqual(result.metadata['ordinary_eligibility_decision'],'provisional_viral')

    def test_Maverick_is_not_reported_as_mixed_domains(self):
        row=te('unknown',strength='qualified',complete='unknown',category='ambiguous_mobile_element')
        row['tesorter_order']='Maverick'
        task,context,catalog=complete_example(row)
        result=harmony.evaluate_candidate(task,context,catalog)
        self.assertEqual(result.metadata['TE_flag'],'ambiguous_mobile_element')
        self.assertEqual(result.metadata['candidate_destination'],'viral_candidates')

    def test_audit_mode_does_not_change_routing_or_scope(self):
        task,context,catalog=complete_example(te(strength='qualified',complete='no'))
        full=harmony.evaluate_candidate(task,context,catalog)
        compact=harmony.evaluate_candidate(task,replace(context,full_audit=False),catalog)
        self.assertEqual(full.metadata,compact.metadata)
        self.assertEqual(full.routing_audit,compact.routing_audit)
        self.assertEqual(compact.evidence_audit,[])

    def test_reporting_cap_does_not_change_destination(self):
        task,context,catalog=complete_example(te())
        for rank in harmony.RANKS:
            result=harmony.evaluate_candidate(task,replace(context,max_rank=rank),catalog)
            self.assertEqual(result.metadata['candidate_destination'],'TE_candidates')

    def test_input_records_not_mutated_by_integrated_routing(self):
        task,context,catalog=complete_example(te(strength='qualified',complete='no'))
        original=deepcopy(task)
        harmony.evaluate_candidate(task,context,catalog)
        self.assertEqual(task,original)

    def test_no_TE_still_has_original_ordinary_eligibility(self):
        task,context,catalog=complete_example()
        result=harmony.evaluate_candidate(task,context,catalog)
        self.assertEqual(result.metadata['viral_decision'],'provisional_viral')
        self.assertEqual(result.routing_audit['pre_TE_routing_decision'],'provisional_viral')
        self.assertEqual(result.metadata['TE_flag'],'NA')

    def test_review_annotation_is_preserved_for_excluded_TE_conflict(self):
        task,context,catalog=complete_example(te(strength='qualified',complete='no'),cellular=True)
        result=harmony.evaluate_candidate(task,context,catalog)
        self.assertIsNotNone(result.review)
        self.assertEqual(result.review['candidate_destination'],'excluded')

    def test_process_publishes_both_candidate_reports(self):
        text=(Path(__file__).resolve().parents[1]/'modules/local/viharmony.nf').read_text()
        for suffix in ('TE_candidates.fasta','TE_candidates.summary.tsv','viral_candidates.fasta','viral_candidates.summary.tsv'):
            self.assertIn('${prefix}.'+suffix,text)
        self.assertIn('emit: candidate_reports',text)
        self.assertIn('overwrite: true',text)


if __name__=='__main__':
    unittest.main(verbosity=2)
