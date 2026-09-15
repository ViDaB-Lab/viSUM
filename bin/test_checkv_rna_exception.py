#!/usr/bin/env python3
"""Experimental terminal-decision replay; does not modify production policy.

Uses verified restored-context outputs and saved per-locus evidence. Preserves
CheckV conflict annotations in a separate exception audit. No sequence-ID rules.
"""
import csv
import json
from pathlib import Path
import argparse
import run_viharmony as harmony

TIERS = {
    'permissive': (30, 50, 50),
    'balanced': (50, 80, 70),
    'strict': (80, 90, 80),
}

def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def names(v):
    return {x for x in v.split(',') if x and x != 'NA'}

def eligible(row, loci, tier):
    if row['record_type'] != 'input_contig' or row['viral_decision'] != 'ambiguous_review':
        return []
    if names(row['cellular_conflict_tools']) != {'checkv'} or names(row['plasmid_conflict_tools']):
        return []
    qualified = names(row['qualified_tools']) | names(row['strong_tools'])
    if not {'deep6', 'virbot', 'vicat'} <= qualified:
        return []
    if row['sequence_interpretation'] in {'likely_retroelement', 'viral_retroelement_conflict', 'ambiguous_mobile_element'}:
        return []
    if any(l['locus_classification'] in {'cellular_supported', 'ambiguous'} for l in loci):
        return []
    identity, coverage, input_span = tier
    accepted=[]
    for l in loci:
        if l['locus_classification'] != 'viral_supported':
            continue
        a,b=map(int,l['coordinates'].split('-'))
        def n(k): return float(l[k]) if l[k] not in ('','NA') else 0
        if (n('protein_length') >= 100 and n('best_bitscore') >= 100
            and n('best_identity') >= identity
            and min(n('best_query_coverage'), n('best_subject_coverage')) >= coverage
            and 100*(b-a+1)/int(row['original_length']) >= input_span):
            accepted.append(l)
    return accepted

def evaluate(saved, evidence_root):
    summary=[];audit=[]
    for folder in sorted(saved.iterdir()):
        if not folder.is_dir():continue
        s=folder.name; output=folder/'restored_context'
        metadata=read(output/f'{s}.final_metadata.tsv')
        primary=read(output/f'{s}.database_candidates.tsv')
        baseline={r['normalized_name'] for r in primary}
        n=len(read(output/f'{s}.sequence_disposition.tsv'))
        loci_by_parent={}
        if 'RNA' in s:
            for l in read(evidence_root/f'{s}_results'/'vicat'/f'{s}.vicat_orf_evidence.tsv'):
                loci_by_parent.setdefault(l['sequence_id'],[]).append(l)
        result={'sample':s,'n':n,'baseline_primary':len(baseline)}
        for name,tier in TIERS.items():
            rescued=set()
            for row in metadata:
                matches=eligible(row,loci_by_parent.get(row['normalized_name'],[]),tier) if 'RNA' in s else []
                if not matches:continue
                common=(row['discovery_status'],row['record_type'],names(row['strong_tools']),names(row['qualified_tools']))
                before=harmony.adjudicate_viral_decision(*common,['checkv'],[],row['sequence_interpretation'])
                after=harmony.adjudicate_viral_decision(*common,[],[],row['sequence_interpretation'])
                if before!='ambiguous_review' or after!='retained_viral':
                    raise ValueError(f'Unexpected adjudication {s} {row["normalized_name"]}')
                # All eligible records are intact, non-mobile and meet a stricter
                # protein floor than the existing RNA pair floor. No region or
                # taxonomy result is reassigned. This is a terminal what-if only.
                rescued.add(row['normalized_name'])
                audit.append(dict(sample=s,tier=name,parent=row['normalized_name'],original_id=row['original_contig_name'],
                    previous_decision=before,experimental_decision=after,
                    preserved_cellular_conflict_tools=row['cellular_conflict_tools'],
                    supporting_loci=matches))
            if rescued & baseline:raise ValueError('Rescued record was already primary')
            result[name]={'primary':len(baseline|rescued),'rescued':len(rescued)}
        summary.append(result)
    return summary,audit

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--context-root',type=Path,required=True)
    p.add_argument('--evidence-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();args.output_dir.mkdir(parents=True,exist_ok=False)
    summary,audit=evaluate(args.context_root,args.evidence_root)
    (args.output_dir/'summary.json').write_text(json.dumps(summary,indent=2))
    (args.output_dir/'exception_audit.json').write_text(json.dumps(audit,indent=2))
    (args.output_dir/'thresholds.json').write_text(json.dumps({'columns':['minimum_identity_percent','minimum_query_and_subject_coverage_percent','minimum_ORF_span_percent_of_input'],'tiers':TIERS,'min_protein_aa':100,'min_bitscore':100},indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
