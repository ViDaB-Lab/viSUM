import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from standardize_virsorter2 import hallmark_only_evidence
from run_viharmony import rna_floor_qualified_tools, apply_rna_pair_homology_floor, evidence_applies_to_final

def evidence():
    return hallmark_only_evidence(dict(seqname='p||lt2gene', hallmark='1', length='300'), 's', 'p', 400)

def test_native_evidence_is_review_only_unscored_and_unlocalized():
    row=evidence()
    assert row['classification']=='virus' and row['evidence_strength']=='review'
    assert row['score']=='NA' and row['coordinates']=='' and row['max_score_group']=='NA'

def test_hallmark_required():
    with pytest.raises(ValueError):
        hallmark_only_evidence(dict(seqname='p||lt2gene',hallmark='0',length='300'),'s','p',400)

def test_no_rna_floor_bypass():
    original={'deep6','vicat','virsorter2'}
    qualified=rna_floor_qualified_tools([evidence()],original)
    assert original=={'deep6','vicat','virsorter2'}
    assert qualified=={'deep6','vicat'}
    assert apply_rna_pair_homology_floor('retained_viral',True,'rna',set(),qualified,[],'p')==('provisional_viral','demoted')

def test_scored_vs2_support_still_exempts_pair():
    row=dict(evidence(),strength_basis='virsorter2_default_cutoff',evidence_strength='qualified')
    assert rna_floor_qualified_tools([evidence(),row],{'deep6','vicat','virsorter2'})=={'deep6','vicat','virsorter2'}

def test_unlocalized_hallmark_does_not_leak_to_child_region():
    args=('p||lt2gene','p','short_hallmark_region',None)
    assert evidence_applies_to_final(*args,'p','p','input_contig',None,400,{'p'})
    assert not evidence_applies_to_final(*args,'child','p','provirus',(1,300),300,{'child'})

def test_hallmark_cannot_supply_a_primary_vote():
    from run_viharmony import adjudicate_viral_decision
    from evidence_schema import is_review_only_evidence
    row=evidence()
    assert is_review_only_evidence(row)
    qualified={'deep6'}
    if not is_review_only_evidence(row):qualified.add(row['tool'])
    assert adjudicate_viral_decision('likely_viral','input_contig',set(),qualified,[],[],'')=='provisional_viral'

def test_cached_modules_include_policy_change():
    root=Path(__file__).resolve().parents[1]
    assert 'Adapter policy v1.2' in (root/'modules/local/standardize_virsorter2.nf').read_text()
    assert 'viHARMONY decision policy v1.2' in (root/'modules/local/viharmony.nf').read_text()
    assert 'Routing policy v1.2' in (root/'modules/local/discovery_gate.nf').read_text()

@pytest.mark.parametrize('strength', ['review', 'qualified', 'strong'])
def test_gate_ignores_hallmark_review_rows_even_with_legacy_strength(tmp_path, strength):
    import csv
    from discovery_gate import load_evidence
    row=dict(evidence(),evidence_strength=strength)
    path=tmp_path/'evidence.tsv'
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
    assert not load_evidence([path],'s',{'p':{'length':'400'}})
