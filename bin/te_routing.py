"""Validated TE candidate routing and descriptive evidence scope.

Pure, local functions: no searches, sequence modification or taxonomy inference.
The default reporting rank never changes the uncapped routing comparison.
"""
from __future__ import annotations

RANKS = ('domain', 'realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
ELIGIBLE = {'retained_viral', 'retained_provirus', 'provisional_viral', 'provisional_provirus'}
TE_VETO = {'likely_retroelement', 'viral_retroelement_conflict', 'ambiguous_mobile_element'}
FAMILIES = {'Gypsy': 'Metaviridae', 'Copia': 'Pseudoviridae', 'Bel-Pao': 'Belpaoviridae'}
ANCESTRY = dict(domain='Viruses', realm='Riboviria', kingdom='Pararnavirae',
                phylum='Artverviricota', **{'class': 'Revtraviricetes'}, order='Ortervirales')
MISSING = {'', 'NA', 'unknown', 'none', 'unclassified'}
LOCAL_PATTERNS = {'localized_viral_cluster', 'mixed_host_viral_signal'}
QUALIFIED = {'qualified', 'strong'}
POLICY_VERSION = 'TE-routing-v1'

METADATA_COLUMNS = [
    'pre_TE_routing_decision', 'ordinary_eligibility_decision',
    'candidate_destination', 'TE_flag', 'TE_architecture_status', 'TE_evidence_scope',
    'TE_routing_reason', 'viral_evidence_scope', 'viral_scope_interpretation',
    'record_level_viral_prediction_tools', 'viral_cluster_coordinates',
    'viral_cluster_coordinate_system',
]
AUDIT_COLUMNS = [
    'sample_id', 'final_sequence_id', 'original_sequence_id', 'normalized_parent_id',
    'record_type', 'coordinates', 'length_nt', 'pre_TE_routing_decision',
    'pre_TE_routing_confidence', 'ordinary_eligibility_decision', 'viral_confidence',
    'candidate_destination', 'previous_proposed_flag', 'TE_flag', 'TE_routing_reason',
    'TE_competing_consensus', 'TE_nonselected_predictions', 'TE_architecture_status',
    'TE_evidence_scope', 'TE_flag_explanation', 'candidate_record_scope',
    'viral_evidence_scope', 'viral_scope_interpretation', 'record_level_viral_prediction_tools',
    'vicat_origin_pattern', 'vicat_viral_supported_loci', 'vicat_cellular_supported_loci',
    'vicat_ambiguous_loci', 'vicat_uninformative_loci', 'vicat_total_loci',
    'viral_cluster_coordinates', 'viral_cluster_coordinate_system',
    'reported_taxonomy', 'reported_taxonomy_rank', 'uncapped_taxonomy',
    'strong_tools', 'qualified_tools', 'cellular_conflict_tools', 'plasmid_conflict_tools',
    'tesorter_strength', 'tesorter_orders', 'tesorter_superfamilies', 'tesorter_domains',
]
SUMMARY_COLUMNS = [
    'sample_id', 'final_sequence_id', 'original_sequence_id', 'normalized_parent_id',
    'record_type', 'coordinates', 'length_nt', 'ordinary_eligibility_decision',
    'viral_confidence', 'candidate_destination', 'TE_flag', 'TE_architecture_status',
    'TE_evidence_scope', 'tesorter_orders', 'tesorter_superfamilies',
    'viral_evidence_scope', 'viral_scope_interpretation',
    'vicat_viral_supported_loci', 'vicat_cellular_supported_loci', 'vicat_total_loci',
    'viral_cluster_coordinates', 'viral_cluster_coordinate_system',
    'strong_tools', 'qualified_tools', 'reported_taxonomy', 'reported_taxonomy_rank',
]


def lineage(value):
    result = {}
    for rank, item in zip(RANKS, value.split(';'), strict=True):
        name = item.split('__', 1)[-1]
        result[rank] = '' if name in ('unclassified', 'NA', '') else name
    return result


def te_profile(rows):
    direct = [r for r in rows if r.get('assignment_method') == 'direct_hmm'
              and r.get('evidence_strength') in ('qualified', 'strong')]
    full = [r for r in direct if r.get('evidence_strength') == 'strong'
            and r.get('tesorter_complete', '').lower() == 'yes']
    identities = {(r.get('tesorter_order'), r.get('tesorter_superfamily')) for r in rows}
    # Clade mixtures within one superfamily do not contradict superfamily identity.
    coherent = bool(full) and len(identities) == 1 and all(
        r.get('tesorter_order') == 'LTR' and r.get('tesorter_superfamily') in FAMILIES
        and r.get('tesorter_window_was_split', 'false') != 'true' for r in rows)
    viral = any(r.get('tesorter_evidence_category') == 'viral_like_mobile_element' for r in rows)
    ordinary = any(r.get('tesorter_evidence_category') == 'retroelement' for r in rows)
    ambiguous = any(r.get('tesorter_evidence_category') == 'ambiguous_mobile_element' for r in rows)
    return dict(coherent=coherent, superfamily=full[0]['tesorter_superfamily'] if coherent else '',
                viral=viral, ordinary=ordinary, ambiguous=ambiguous, direct=bool(direct))


def route(baseline, eligibility, te_rows, uncapped, votes):
    """Export decision only; input mappings and their taxonomies remain intact."""
    destination = 'viral_candidates' if eligibility['viral_decision'] in ELIGIBLE else (
        'excluded' if eligibility['viral_decision'] == 'likely_nonviral' else 'review_only')
    if not te_rows:
        return dict(destination=destination, flag='NA', reason='no_TE_evidence',
                    competing_consensus='NA', nonselected_predictions='NA')
    profile = te_profile(te_rows)
    rejected = sorted({f"{v['tool']}:{v['rank']}:{v['taxon']}:{v['tier']}" for v in votes
                       if str(v.get('valid')).lower() == 'true'
                       and str(v.get('accepted')).lower() == 'false'
                       and v.get('tier') in ('moderate', 'high')})
    conflicting = []
    if profile['coherent']:
        expected = {**ANCESTRY, 'family': FAMILIES[profile['superfamily']]}
        actual = lineage(uncapped['uncapped_strict_taxonomy'])
        conflicting = [f'{rank}:{actual[rank]}' for rank in expected
                       if actual.get(rank) and actual[rank] != expected[rank]]
        if not conflicting:
            destination, flag, reason = ('TE_candidates', 'TE_supported',
                'strong_coherent_TE;accepted_consensus_compatible_or_broader')
        else:
            flag, reason = 'TE_viral_conflict', 'strong_TE;accepted_consensus_competes'
    elif profile['ambiguous'] or (profile['viral'] and profile['ordinary']):
        flag, reason = 'mixed_mobile_element_evidence', 'mixed_or_ambiguous_domains;ordinary_eligibility_preserved'
    elif profile['viral']:
        flag, reason = 'retrovirus_compatible', 'explicit_viral_like_element;endogenous_status_unresolved'
    elif not profile['direct']:
        flag, reason = 'weak_TE_signal', 'similarity_only_TE_signal;ordinary_eligibility_preserved'
    else:
        flag, reason = 'partial_TE_support', 'partial_or_localized_TE_domains;ordinary_eligibility_preserved'
    if destination not in ('viral_candidates', 'TE_candidates'):
        reason += ';ordinary_origin_or_eligibility_safeguard_blocks_export'
    return dict(destination=destination, flag=flag, reason=reason,
                competing_consensus=','.join(conflicting) or 'NA',
                nonselected_predictions=','.join(rejected) or 'NA')


def corrected_route(baseline, eligibility, te_rows, uncapped, votes):
    result = route(baseline, eligibility, te_rows, uncapped, votes)
    profile = te_profile(te_rows)
    rank = uncapped.get('uncapped_strict_taxonomy_conflict_rank', 'NA')
    is_conflict = str(uncapped.get('uncapped_taxonomy_conflict', 'false')).lower() == 'true'
    if (result['destination'] == 'TE_candidates' and profile['coherent'] and is_conflict
            and rank in RANKS[:RANKS.index('family')+1]):
        expected = {**ANCESTRY, 'family': FAMILIES[profile['superfamily']]}
        alternatives = [name.strip() for name in uncapped.get('uncapped_strict_taxonomy_conflicting_taxa', '').split(',')
                        if name.strip() not in ('', 'NA', expected.get(rank))]
        if alternatives:
            result.update(destination=('viral_candidates' if eligibility['viral_decision'] in ELIGIBLE else
                'excluded' if eligibility['viral_decision'] == 'likely_nonviral' else 'review_only'),
                flag='TE_viral_conflict', reason='strong_TE;consensus_records_unresolved_competing_taxa',
                competing_consensus=rank+':'+','.join(alternatives))
    return result


def present(value):
    return str(value).strip().lower() not in {v.lower() for v in MISSING}


def joined(values):
    return ','.join(sorted({str(v) for v in values if present(v)})) or 'NA'


def integer(value):
    """Missing counts stay missing; malformed or negative counts must fail."""
    if value is None or not present(value):
        return None
    number = int(value)
    if number < 0:
        raise ValueError('Negative locus/domain count')
    return number


def active_rows(rows, tool=None):
    return [r for r in rows if r.get('evidence_application') == 'used_for_final_decision'
            and (tool is None or r.get('tool') == tool)]


def clarify_te(previous_flag, rows):
    """Describe architecture/assignment uncertainty, without changing routing."""
    rows = active_rows(rows, 'tesorter')
    if not rows:
        return dict(proposed_flag='NA', TE_architecture_status='NA',
                    TE_evidence_scope='NA', TE_flag_explanation='no_TE_evidence')
    orders = {r.get('tesorter_order', '') for r in rows if present(r.get('tesorter_order', ''))}
    superfamilies = {r.get('tesorter_superfamily', '') for r in rows
                     if present(r.get('tesorter_superfamily', ''))}
    # A clade mixture within one superfamily is not a conflicting TE identity.
    explicit_mix = (any(v.lower() == 'mixture' for v in orders | superfamilies)
                    or len(orders) > 1 or len(superfamilies) > 1)
    categories = {r.get('tesorter_evidence_category', '') for r in rows}
    ambiguous = 'ambiguous_mobile_element' in categories
    localized = any(r.get('tesorter_window_was_split', '').lower() == 'true' for r in rows)
    direct = any(r.get('assignment_method') == 'direct_hmm'
                 and r.get('evidence_strength') in QUALIFIED for r in rows)
    complete = {r.get('tesorter_complete', '').lower() for r in rows}
    architecture = ('mixed_assignments' if explicit_mix else 'incomplete' if 'no' in complete
                    else 'complete' if 'yes' in complete else 'not_established')
    scope = 'split_window_support' if localized else 'candidate_record_assessed'
    # Preserve the guarded router's substantive TE/consensus-conflict decisions.
    if previous_flag in {'TE_supported', 'TE_viral_conflict', 'retrovirus_compatible', 'weak_TE_signal'}:
        flag = previous_flag
        explanation = {
            'TE_supported': 'strong_complete_coherent_TE;guarded_destination_unchanged',
            'TE_viral_conflict': 'strong_TE_vs_accepted_or_explicitly_conflicting_consensus',
            'retrovirus_compatible': 'viral_like_domains;endogenous_or_exogenous_status_unresolved',
            'weak_TE_signal': 'no_qualified_direct_domain_support',
        }[flag]
    elif explicit_mix or {'retroelement', 'viral_like_mobile_element'} <= categories:
        flag, explanation = 'mixed_domain_assignments', 'different_TE_or_retroviral_assignments;not_proof_of_a_hybrid'
    elif ambiguous:
        flag, explanation = 'ambiguous_mobile_element', 'ambiguous_element_identity;multiple_conflicting_domains_not_required'
    elif localized:
        flag, explanation = 'localized_TE_domain_support', 'split_window_TE_evidence;whole_record_TE_identity_not_established'
    elif direct and 'no' in complete:
        flag, explanation = 'incomplete_TE_architecture', 'direct_TE_domains;TEsorter_expected_architecture_incomplete'
    else:
        flag, explanation = 'limited_TE_domain_support', 'TE_domains_present;complete_architecture_not_established'
    return dict(proposed_flag=flag, TE_architecture_status=architecture,
                TE_evidence_scope=scope, TE_flag_explanation=explanation)


def viral_scope(record, evidence):
    """Use actual locus patterns, never audit 'region' alone, as a spatial claim.

    Projected viCAT rows may describe an intact input contig despite their
    applied_scope='region'. Coordinates from viCAT clusters are 1-based inclusive
    on the current candidate, not necessarily the original parent scaffold.
    """
    rows = active_rows(evidence)
    viral = [r for r in rows if r.get('tool') != 'tesorter'
             and r.get('classification') == 'virus' and r.get('evidence_strength') in QUALIFIED]
    vicat = [r for r in viral if r.get('tool') == 'vicat']
    # Actual fixtures have at most one applied viCAT summary per candidate.
    # Refuse to add repeated/projected summaries as if they were unique loci.
    if len(vicat) > 1:
        raise ValueError('Multiple active viCAT summaries: scope aggregation requires review')
    row = vicat[0] if vicat else {}
    pattern = row.get('origin_pattern') or 'NA'
    basis = row.get('strength_basis', '')
    local = pattern in LOCAL_PATTERNS or basis in {
        'vicat_projected_localized_viral_cluster', 'vicat_projected_mixed_host_viral_signal'}
    record_tools = joined(r['tool'] for r in viral if r['tool'] != 'vicat')
    count = integer(row.get('viral_supported_loci'))
    if local:
        scope = ('localized_cluster_only' if record_tools == 'NA'
                 else 'localized_cluster_with_record_level_predictions')
        note = 'localized_viral_content;entire_candidate_viral_origin_not_established'
    elif count is not None and count >= 2:
        scope, note = 'multi_locus_viral_support', 'multiple_supported_loci;not_basewise_or_complete_genome_validation'
    elif count == 1:
        scope, note = 'single_locus_viral_support', 'one_supported_locus;not_independent_whole_sequence_confirmation'
    elif viral:
        scope, note = 'record_level_predictions_only', 'record_level_predictions;spatial_extent_of_viral_content_unresolved'
    else:
        scope, note = 'no_qualified_viral_support', 'no_qualified_viral_origin_row_in_applied_audit'
    coordinates = row.get('viral_cluster_coordinates') or 'NA'
    record_type = record.get('record_type', 'NA')
    return dict(
        candidate_record_scope='input_contig' if record_type == 'input_contig' else record_type,
        viral_evidence_scope=scope,
        viral_scope_interpretation=note,
        record_level_viral_prediction_tools=record_tools,
        vicat_origin_pattern=pattern,
        vicat_viral_supported_loci=count if count is not None else 'NA',
        vicat_cellular_supported_loci=(integer(row.get('cellular_supported_loci'))
                                      if integer(row.get('cellular_supported_loci')) is not None else 'NA'),
        vicat_ambiguous_loci=integer(row.get('ambiguous_loci')) if integer(row.get('ambiguous_loci')) is not None else 'NA',
        vicat_uninformative_loci=(integer(row.get('uninformative_loci'))
                                 if integer(row.get('uninformative_loci')) is not None else 'NA'),
        vicat_total_loci=integer(row.get('orf_loci')) if integer(row.get('orf_loci')) is not None else 'NA',
        viral_cluster_coordinates=coordinates if present(coordinates) else 'NA',
        viral_cluster_coordinate_system=('candidate_record_1_based_inclusive' if present(coordinates) else 'NA'),
    )


def describe_candidate(metadata, evidence, reporting, votes, sample_id,
                       previous_decision, previous_confidence):
    te_rows = active_rows(evidence, 'tesorter')
    selected = corrected_route({}, metadata, te_rows, reporting, votes)
    clarified = clarify_te(selected['flag'], evidence)
    scope = viral_scope(metadata, evidence)
    audit = dict(sample_id=sample_id,
        final_sequence_id=metadata['final_sequence_id'],
        original_sequence_id=metadata['original_contig_name'],
        normalized_parent_id=metadata['normalized_name'],
        record_type=metadata['record_type'], coordinates=metadata['provirus_coordinates'],
        length_nt=metadata['refined_length'], pre_TE_routing_decision=previous_decision,
        pre_TE_routing_confidence=previous_confidence,
        ordinary_eligibility_decision=metadata['viral_decision'],
        viral_confidence=metadata['viral_confidence'], candidate_destination=selected['destination'],
        previous_proposed_flag=selected['flag'], TE_flag=clarified['proposed_flag'],
        TE_routing_reason=selected['reason'], TE_competing_consensus=selected['competing_consensus'],
        TE_nonselected_predictions=selected['nonselected_predictions'],
        TE_architecture_status=clarified['TE_architecture_status'],
        TE_evidence_scope=clarified['TE_evidence_scope'], TE_flag_explanation=clarified['TE_flag_explanation'],
        **scope,
        reported_taxonomy=metadata['strict_taxonomy'], reported_taxonomy_rank=metadata['strict_taxonomy_rank'],
        uncapped_taxonomy=reporting['uncapped_strict_taxonomy'],
        strong_tools=metadata['strong_tools'], qualified_tools=metadata['qualified_tools'],
        cellular_conflict_tools=metadata['cellular_conflict_tools'], plasmid_conflict_tools=metadata['plasmid_conflict_tools'],
        tesorter_strength=metadata['tesorter_evidence_strength'], tesorter_orders=metadata['tesorter_orders'],
        tesorter_superfamilies=metadata['tesorter_superfamilies'],
        tesorter_domains=' / '.join(r.get('domains', 'NA') for r in te_rows) or 'NA')
    if set(audit) != set(AUDIT_COLUMNS):
        raise ValueError('TE routing audit schema mismatch')
    return audit
