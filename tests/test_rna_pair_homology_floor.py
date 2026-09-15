import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from run_viharmony import apply_rna_pair_homology_floor


def call(loci, **overrides):
    args = dict(decision='retained_viral', enabled=True, input_type='rna',
                strong_tools=set(), qualified_tools={'deep6', 'vicat'},
                loci=loci, parent_id='p')
    args.update(overrides)
    return apply_rna_pair_homology_floor(**args)


def locus(**overrides):
    row = dict(sequence_id='p', locus_classification='viral_supported',
               coordinates='1-303', protein_length='100', best_bitscore='100',
               best_query_coverage='50', best_subject_coverage='50')
    row.update(overrides)
    return row


def test_inclusive_boundaries_and_any_single_locus():
    assert call([locus()]) == ('retained_viral', 'passed')
    assert call([locus(best_bitscore='99'), locus()])[0] == 'retained_viral'


@pytest.mark.parametrize('field,value', [('protein_length','99'), ('best_bitscore','99.99'),
    ('best_query_coverage','49.99'), ('best_subject_coverage','49.99'),
    ('sequence_id','other'), ('locus_classification','cellular_supported')])
def test_each_required_condition(field, value):
    assert call([locus(**{field:value})]) == ('provisional_viral', 'demoted')


@pytest.mark.parametrize('overrides', [dict(enabled=False), dict(input_type='dna'),
    dict(strong_tools={'genomad'}), dict(qualified_tools={'deep6','vicat','genomad'})])
def test_no_change_outside_scope(overrides):
    assert call([], **overrides) == ('retained_viral', 'not_applicable')


def test_existing_conflict_not_promoted():
    assert call([locus()], decision='ambiguous_review') == ('ambiguous_review','not_applicable')


def test_region_cannot_borrow_parent_locus():
    assert call([locus()], decision='retained_provirus', coordinates='400-800') == ('provisional_provirus','demoted')
    assert call([locus()], decision='retained_provirus', coordinates='1-303')[0] == 'retained_provirus'


def test_missing_loci_demote():
    assert call([]) == ('provisional_viral','demoted')
