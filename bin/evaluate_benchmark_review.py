#!/usr/bin/env python3
"""Evaluate saved viSUM development benchmarks without relabelling suspected proviruses.

Counts standardized pre-gate discovery calls at original-input level. Not a
claim that every tool was run with its author's default configuration.
"""
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path

TOOLS=['genomad','virsorter2','cenotetaker3','deep6','deepmicroclass2','virbot','gianthunter','vicat']

def native_virsorter2_calls(rows, sample):
    """Include native hallmark-only lt2gene calls even without an ML score."""
    result=set()
    for row in rows:
        match=re.fullmatch(re.escape(sample)+r'(__c\d+)(?:\|\|(?:full|lt2gene|\d+_partial))?',row['seqname'])
        if match is None:
            raise ValueError(f"Unrecognized VirSorter2 parent ID: {row['seqname']}")
        result.add(sample+match.group(1))
    return result
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))
def write(path,rr):
    if not rr:return
    keys=list(dict.fromkeys(k for r in rr for k in r))
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(rr)
def evaluate(root,out):
    out.mkdir(parents=True,exist_ok=True)
    on=next(p for p in root.glob('rna_pair_floor_*') if (p/'BMARK_RNA_VIRAL_results').exists())
    off=next(root.glob('rna_pair_OFF_*'))
    metrics=[];calls=[];negative=[];deltas=[];integrity=[];details=[];regression=[]
    for folder in sorted(on.glob('*_results')):
        s=folder.name.removesuffix('_results');h=folder/'viharmony'
        typ='rna' if s.startswith('RNANEG') or s=='BMARK_RNA_VIRAL' else 'dna'
        positive=s in {'BMARK_DNA_VIRAL','BMARK_RNA_VIRAL'}
        disp=read(h/f'{s}.sequence_disposition.tsv')
        mapping={r['normalized_name']:r['original_contig_name'] for r in disp}
        assert len(mapping)==len(disp)==len(set(mapping.values()))
        ids=set(mapping);n=len(ids)
        primary=read(h/f'{s}.database_candidates.tsv')
        meta=read(h/f'{s}.final_metadata.tsv')
        manifest=json.loads((h/f'{s}.harmonizer_manifest.json').read_text())
        for entry in manifest['outputs'].values():
            p=h/entry['file']
            if p.exists():
                assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256'],p
                integrity.append(dict(sample=s,file=entry['file'],status='hash_match'))
            else:integrity.append(dict(sample=s,file=entry['file'],status='not_extracted_or_not_packaged'))
        called={r['normalized_name'] for r in primary}; assert called<=ids
        sets={'viSUM_ON':called}
        if (off/folder.name).exists():
            old=off/folder.name/'viharmony'
            sets['viSUM_OFF']={r['normalized_name'] for r in read(old/f'{s}.database_candidates.tsv')}
            if typ=='dna':
                for suffix in ['final_metadata.tsv','sequence_disposition.tsv','database_candidates.tsv','provisional_metadata.tsv']:
                    assert (old/f'{s}.{suffix}').read_bytes()==(h/f'{s}.{suffix}').read_bytes(),(s,suffix)
                regression.append(dict(sample=s,status='DNA_ON_OFF_identical_for_four_tables'))
        else:
            demoted=set(manifest['policy']['rna_pair_homology_floor_demoted_ids'])
            sets['viSUM_OFF']=called|{r['normalized_name'] for r in meta if r['final_sequence_id'] in demoted}
        for tool in TOOLS:
            p=folder/tool/f'{s}.{tool}_evidence.tsv'
            if not p.exists():
                unsupported=(typ=='dna' and tool in {'deep6','virbot'}) or (typ=='rna' and tool in {'deepmicroclass2','gianthunter'})
                if not unsupported:raise FileNotFoundError(f'Missing applicable discovery output: {p}')
                metrics.append(dict(sample=s,input_type=typ,positive=positive,method=tool,n=n,called='',rate='',status='not_applicable'));continue
            evidence=read(p)
            raw={r['parent_sequence_id'] if r['parent_sequence_id'] not in {'','NA'} else r['sequence_id'] for r in evidence if r['classification']=='virus'}
            assert raw<=ids,(s,tool,raw-ids)
            sets[tool]=raw
            if tool=='virsorter2':
                native=native_virsorter2_calls(read(folder/tool/f'{s}.virsorter2_final-viral-score.tsv'),s)
                assert native<=ids and raw<=native
                sets['virsorter2_native']=native
        for method,selected in sets.items():
            metrics.append(dict(sample=s,input_type=typ,positive=positive,method=method,n=n,called=len(selected),rate=len(selected)/n,status='evaluated'))
            for sid in sorted(ids):calls.append(dict(sample=s,original_id=mapping[sid],normalized_id=sid,positive=positive,method=method,called=int(sid in selected)))
            if method not in {'viSUM_ON','viSUM_OFF'}:
                deltas.append(dict(sample=s,method=method,n=n,both=len(called&selected),visum_only=len(called-selected),tool_only=len(selected-called),neither=len(ids-called-selected)))
        grouped=collections.defaultdict(list)
        for r in primary:grouped[r['normalized_name']].append(r)
        for sid,rr in grouped.items():
            extracted=[r for r in rr if r['record_type'] in {'provirus','viral_region'}]
            if not positive:
                evidence_class='extracted_region_candidate' if extracted else 'intact_retention_no_extracted_region'
                negative.append(dict(sample=s,input_type=typ,original_id=mapping[sid],normalized_id=sid,evidence_class=evidence_class,retained_records=len(rr),extracted_regions=len(extracted),independently_confirmed_provirus='not_established'))
                details.extend(dict(sample=s,evidence_class=evidence_class,**r) for r in rr)
        print(s,n,'viSUM',len(called),'regions',sum(r['record_type']!='input_contig' for r in primary),'TOOLS',[(k,len(v)) for k,v in sets.items() if k not in {'viSUM_ON','viSUM_OFF'}])
    for name,rr in [('metrics.tsv',metrics),('input_calls.tsv',calls),('negative_parent_audit.tsv',negative),('negative_record_audit.tsv',details),('paired_comparison.tsv',deltas),('integrity.tsv',integrity),('regression.tsv',regression)]:write(out/name,rr)
    summary=[]
    for s in sorted({r['sample'] for r in metrics if not r['positive']}):
        rr=[r for r in negative if r['sample']==s]
        n=next(r['n'] for r in metrics if r['sample']==s)
        extracted=sum(r['evidence_class']=='extracted_region_candidate' for r in rr)
        summary.append(dict(sample=s,n=n,all_retained=len(rr),extracted_candidate_parents=extracted,intact_only_parents=len(rr)-extracted,label_based_fpr=len(rr)/n,intact_only_rate_not_corrected_fpr=(len(rr)-extracted)/n))
    write(out/'negative_evidence_summary.tsv',summary)
    print('NEGATIVE_SUMMARY',summary)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--review-root',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();evaluate(a.review_root,a.outdir)
