#!/usr/bin/env python3
"""Isolate restored CheckV context with saved boundaries and taxonomy held fixed.

Both baseline and variant use the same harmonizer. Baseline tables must reproduce
the saved run. This is not a combined validation of a changed trimming policy.
"""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor
import standardize_checkv as standardizer
from refine_proviral_regions import load_fasta

def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def execute(folder, args):
    s=folder.name.removesuffix('_results'); typ='rna' if 'RNA' in s else 'dna'
    dest=args.output_dir/s; dest.mkdir()
    fa=(args.fasta_root/folder.name) if args.fasta_root else folder
    cv=folder/'checkv'
    quality=cv/f'{s}.checkv_quality_summary.tsv'
    wanted={r['contig_id'] for r in read(quality)}
    original=load_fasta(fa/'prep'/f'{s}.normalized.fasta')
    assert wanted<=original.keys()
    candidate=dest/'candidates.fasta'
    with candidate.open('w',encoding='utf-8',newline='\n') as f:
        for sid in sorted(wanted):f.write(f'>{sid}\n{original[sid]}\n')
    restored=dest/'checkv_context.tsv'
    standardizer.run(SimpleNamespace(sample_id=s,input_type=typ,candidate_fasta=candidate,
        quality_summary=quality,completeness=cv/f'{s}.checkv_completeness.tsv',
        contamination=cv/f'{s}.checkv_contamination.tsv',complete_genomes=cv/f'{s}.checkv_complete_genomes.tsv',
        run_metadata=cv/f'{s}.checkv_run_metadata.tsv',output=restored))
    old_cv=read(cv/f'{s}.checkv_evidence.tsv'); new_cv=read(restored)
    new_by_id={r['sequence_id']:r for r in new_cv}
    assert all(new_by_id[r['sequence_id']]==r for r in old_cv),'Existing CheckV evidence changed'
    tools=['genomad','virsorter2','cenotetaker3','deep6','deepmicroclass2','virbot','gianthunter','vicat','checkv','tesorter','vitap','vcontact3']
    evidence=[folder/t/f'{s}.{t}_evidence.tsv' for t in tools if (folder/t/f'{s}.{t}_evidence.tsv').exists()]
    evidence.append(folder/'vicat'/f'{s}.vicat_refined_evidence.tsv')
    if typ=='rna':evidence.append(folder/'vicat'/f'{s}.vicat_orf_evidence.tsv')
    common=['--sample-id',s,'--input-type',typ,'--normalized-fasta',fa/'prep'/f'{s}.normalized.fasta',
        '--refined-fasta',fa/'refinement'/f'{s}.refined_candidates.fasta',
        '--header-map',folder/'prep'/f'{s}.header_map.tsv','--discovery-gate',folder/'discovery_gate'/f'{s}.discovery_gate.tsv',
        '--region-map',folder/'refinement'/f'{s}.provirus_region_map.tsv','--ictv-msl',args.ictv_msl,
        '--vcontact3-groups',folder/'vcontact3'/f'{s}.vcontact3_group_membership.tsv','--audit-mode','full']
    if typ=='rna':common.append('--rna-pair-homology-floor')
    for variant in ['baseline','restored_context']:
        out=dest/variant;out.mkdir()
        ev=[restored if variant=='restored_context' and p==cv/f'{s}.checkv_evidence.tsv' else p for p in evidence]
        cmd=[sys.executable,str(Path(__file__).with_name('run_viharmony.py')),*map(str,common),'--disable-rna-checkv-exception','--output-prefix',str(out/s),'--evidence',*map(str,ev)]
        run=subprocess.run(cmd,capture_output=True,text=True)
        (out/'run.log').write_text(run.stdout+run.stderr,encoding='utf-8')
        if run.returncode:raise RuntimeError(f'{s}: {run.stderr}')
    for suffix in ['database_candidates.tsv','final_metadata.tsv','provisional_metadata.tsv','sequence_disposition.tsv']:
        observed=read(dest/'baseline'/f'{s}.{suffix}')
        saved=read(folder/'viharmony'/f'{s}.{suffix}')
        assert len(observed)==len(saved) and all(all(new.get(k)==v for k,v in old.items()) for new,old in zip(observed,saved)),f'Baseline mismatch {s} {suffix}'
    before=read(dest/'baseline'/f'{s}.database_candidates.tsv');after=read(dest/'restored_context'/f'{s}.database_candidates.tsv')
    old_parents={r['normalized_name'] for r in before};new_parents={r['normalized_name'] for r in after}
    changes={'lost':[r for r in before if r['normalized_name'] not in new_parents], 'gained':[r for r in after if r['normalized_name'] not in old_parents]}
    (dest/'changes.json').write_text(json.dumps(changes,indent=2),encoding='utf-8')
    result=dict(sample=s,baseline_tables_match=4,old_checkv_rows=len(old_cv),new_checkv_rows=len(new_cv),
        primary_before=len(old_parents),primary_after=len(new_parents),lost_parents=len(old_parents-new_parents),gained_parents=len(new_parents-old_parents))
    print(json.dumps(result),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-root',type=Path,required=True)
    p.add_argument('--fasta-root',type=Path)
    p.add_argument('--ictv-msl',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    folders=sorted(a.results_root.glob('*_results'))
    if not folders:raise ValueError('No saved panels found')
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(lambda f:execute(f,a),folders))
    (a.output_dir/'summary.json').write_text(json.dumps(results,indent=2),encoding='utf-8')

if __name__=='__main__':main()
