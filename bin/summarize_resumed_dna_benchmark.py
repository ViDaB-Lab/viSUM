#!/usr/bin/env python3
"""Verify resumed DNA benchmark outputs and report original-label metrics."""
import argparse,csv,json,hashlib
from pathlib import Path

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def fasta(p):
    result={};key=None
    with p.open() as f:
        for line in f:
            if line.startswith('>'):
                key=line[1:].split()[0]
                if key in result:raise ValueError(f'Duplicate FASTA ID: {key}')
                result[key]=''
            elif line.strip():result[key]+=line.strip()
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);a=p.parse_args()
    panels=['CELLULAR','MITOCHONDRIAL','PLASMID','PLASTID','RETROELEMENT','VIRAL']
    totals={};summary=[]
    for panel in panels:
        s='BMARK_DNA_'+panel;root=a.results/f'{s}_results';h=root/'viharmony'
        manifest=json.loads((h/f'{s}.harmonizer_manifest.json').read_text())
        for item in manifest['outputs'].values():
            file=(h/item['file']).resolve()
            if not file.is_relative_to(h.resolve()):raise ValueError('Unsafe manifest output path')
            if hashlib.sha256(file.read_bytes()).hexdigest()!=item['sha256']:raise ValueError(f'Hash mismatch: {file}')
        refinement=read(root/'refinement'/f'{s}.provirus_refinement_summary.tsv')
        if len(refinement)!=1 or refinement[0]['allow_ct3_only_refinement']!='false':raise ValueError(f'Wrong trimming policy: {s}')
        disposition=read(h/f'{s}.sequence_disposition.tsv');ids={r['normalized_name'] for r in disposition}
        if len(ids)!=len(disposition):raise ValueError('Duplicate input IDs')
        primary=read(h/f'{s}.database_candidates.tsv');parents={r['normalized_name'] for r in primary}
        if not parents<=ids:raise ValueError('Unknown retained parent')
        normalized=fasta(root/'prep'/f'{s}.normalized.fasta');sequences=fasta(h/f'{s}.database_candidates.fasta')
        if set(normalized)!=ids or set(sequences)!={r['final_sequence_id'] for r in primary}:raise ValueError('FASTA/table ID mismatch')
        for row in primary:
            if row['rna_checkv_conflict_exception']!='false':raise ValueError('RNA exception applied to DNA')
            expected=normalized[row['normalized_name']]
            if row['record_type']!='input_contig':
                start,end=map(int,row['provirus_coordinates'].split('-'));expected=expected[start-1:end]
            actual=sequences[row['final_sequence_id']]
            if actual.upper()!=expected.upper() or len(actual)!=int(row['refined_length']):raise ValueError('Incorrect retained sequence')
        intact={r['normalized_name'] for r in primary if r['record_type']=='input_contig'}
        result=dict(panel=panel,inputs=len(ids),primary_parents=len(parents),intact_primary_parents=len(intact),checks='PASS')
        summary.append(result);totals[panel]=result
    positive=totals['VIRAL'];neg=[r for r in summary if r['panel']!='VIRAL']
    nn=sum(r['inputs'] for r in neg);fp=sum(r['primary_parents'] for r in neg)
    print(json.dumps(dict(panels=summary,dna_sensitivity_percent=100*positive['primary_parents']/positive['inputs'],
        dna_intact_positive_recovery_percent=100*positive['intact_primary_parents']/positive['inputs'],
        dna_nominal_negative_retained=fp,dna_nominal_negative_total=nn,dna_original_label_fpr_percent=100*fp/nn,
        note='No virus-like relabeling; regional-content review remains a separate endpoint.'),indent=2))

if __name__=='__main__':main()
