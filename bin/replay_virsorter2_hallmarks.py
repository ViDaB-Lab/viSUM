#!/usr/bin/env python3
"""Replay the VS2 adapter and viHARMONY against saved, fixed downstream evidence.

Never runs discovery tools. Reports newly gated inputs separately: those require
downstream follow-up, so fixed-candidate results are not full-pipeline estimates.
"""
import argparse
import csv
import json
import subprocess
import sys
import tarfile
from pathlib import Path

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write(path,rows):
    if not rows:return
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def command(args):
    result=subprocess.run([sys.executable,*map(str,args)],capture_output=True,text=True)
    if result.returncode:raise RuntimeError(result.stdout+result.stderr)
    print(result.stdout.strip(),flush=True)

def run(a):
    root=a.review_root.resolve();out=a.outdir.resolve()
    out.mkdir(parents=True,exist_ok=False)
    on=next(root.glob('rna_pair_floor_*')); fasta=out/'fasta';fasta.mkdir()
    if a.archive:
        wanted=set()
        for p in on.glob('*_results'):
            s=p.name.removesuffix('_results');prefix=p.relative_to(root).as_posix()
            wanted.update([f'{prefix}/prep/{s}.normalized.fasta',f'{prefix}/refinement/{s}.refined_candidates.fasta'])
        with tarfile.open(a.archive,'r|gz') as tf:
            for member in tf:
                if member.name in wanted:
                    tf.extract(member,fasta,filter='data');wanted.remove(member.name)
                    print('Extracted',member.name,flush=True)
                if not wanted:break
        if wanted:raise FileNotFoundError(str(wanted))
    else:fasta=a.fasta_root.resolve() if a.fasta_root else root
    local=Path(__file__).resolve().parent
    summary=[];changes=[];newgate=[]
    for folder in sorted(on.glob('*_results')):
        s=folder.name.removesuffix('_results');typ='rna' if 'RNA' in s else 'dna'
        dest=out/s;dest.mkdir();vs=folder/'virsorter2';h=folder/'viharmony'
        command([local/'standardize_virsorter2.py','--sample-id',s,'--input-type',typ,
            '--header-map',folder/'prep'/f'{s}.header_map.tsv',
            '--score-table',vs/f'{s}.virsorter2_final-viral-score.tsv',
            '--boundary-table',vs/f'{s}.virsorter2_final-viral-boundary.tsv',
            '--run-metadata',vs/f'{s}.virsorter2_run_metadata.tsv',
            '--output',dest/f'{s}.virsorter2_evidence.tsv'])
        before=read(vs/f'{s}.virsorter2_evidence.tsv');after=read(dest/f'{s}.virsorter2_evidence.tsv')
        oldids={r['sequence_id'] for r in before}
        added=[r for r in after if r['sequence_id'] not in oldids]
        regions=read(folder/'refinement'/f'{s}.provirus_region_map.tsv')
        existing={r['parent_sequence_id'] or r['sequence_id'] for r in regions}
        for r in added:
            if r['evidence_strength'] in {'qualified','strong'} and r['parent_sequence_id'] not in existing:
                newgate.append(dict(sample=s,parent_id=r['parent_sequence_id'],reason='new_discovery_candidate_requires_downstream_followup'))
        evidence=[]
        for tool in ['genomad','virsorter2','cenotetaker3','deep6','deepmicroclass2','virbot','gianthunter','vicat','checkv','tesorter','vitap','vcontact3']:
            p=folder/tool/f'{s}.{tool}_evidence.tsv'
            if p.exists():evidence.append(p)
        evidence.append(folder/'vicat'/f'{s}.vicat_refined_evidence.tsv')
        if typ=='rna':evidence.append(folder/'vicat'/f'{s}.vicat_orf_evidence.tsv')
        gate=dest/f'{s}.discovery_gate.tsv'
        discovery=[dest/f'{s}.virsorter2_evidence.tsv' if p.name==f'{s}.virsorter2_evidence.tsv' else p for p in evidence if p.parent.name in {'genomad','virsorter2','cenotetaker3','deep6','deepmicroclass2','virbot','gianthunter','vicat'} and not p.name.endswith(('.vicat_refined_evidence.tsv','.vicat_orf_evidence.tsv'))]
        command([local/'discovery_gate.py','--sample-id',s,'--input-type',typ,
            '--fasta',fasta/folder.relative_to(root)/'prep'/f'{s}.normalized.fasta',
            '--header-map',folder/'prep'/f'{s}.header_map.tsv','--evidence',*discovery,
            '--output-candidates',dest/'gate_candidates.fasta','--output-noncandidates',dest/'gate_noncandidates.fasta',
            '--output-audit',gate,'--output-summary',dest/'gate_summary.tsv'])
        if a.require_unchanged:
            assert gate.read_bytes()==(folder/'discovery_gate'/f'{s}.discovery_gate.tsv').read_bytes(),f'Gate changed: {s}'
        common=['--sample-id',s,'--input-type',typ,
            '--normalized-fasta',fasta/folder.relative_to(root)/'prep'/f'{s}.normalized.fasta',
            '--refined-fasta',fasta/folder.relative_to(root)/'refinement'/f'{s}.refined_candidates.fasta',
            '--header-map',folder/'prep'/f'{s}.header_map.tsv',
            '--discovery-gate',folder/'discovery_gate'/f'{s}.discovery_gate.tsv',
            '--region-map',folder/'refinement'/f'{s}.provirus_region_map.tsv',
            '--ictv-msl',root/'assets/ICTV_VMR_MSL41.csv','--audit-mode','full',
            '--vcontact3-groups',folder/'vcontact3'/f'{s}.vcontact3_group_membership.tsv']
        if typ=='rna':common+=['--rna-pair-homology-floor']
        for label in ['baseline','updated']:
            runout=dest/label;runout.mkdir()
            ev=[dest/f'{s}.virsorter2_evidence.tsv' if label=='updated' and p.name==f'{s}.virsorter2_evidence.tsv' else p for p in evidence]
            command([local/'run_viharmony.py',*common,'--evidence',*ev,'--output-prefix',runout/s])
        for suffix in ['database_candidates.tsv','final_metadata.tsv','provisional_metadata.tsv','sequence_disposition.tsv']:
            assert (dest/'baseline'/f'{s}.{suffix}').read_bytes()==(h/f'{s}.{suffix}').read_bytes(),f'Baseline does not reproduce archive: {s}.{suffix}'
            if a.require_unchanged:
                assert (dest/'updated'/f'{s}.{suffix}').read_bytes()==(h/f'{s}.{suffix}').read_bytes(),f'Final policy changed decisions: {s}.{suffix}'
        old={r['normalized_name'] for r in read(dest/'baseline'/f'{s}.database_candidates.tsv')}
        new={r['normalized_name'] for r in read(dest/'updated'/f'{s}.database_candidates.tsv')}
        mapping={r['normalized_name']:r['original_contig_name'] for r in read(h/f'{s}.sequence_disposition.tsv')}
        for sid in sorted(old^new):changes.append(dict(sample=s,normalized_id=sid,original_id=mapping[sid],change='gained' if sid in new else 'lost'))
        summary.append(dict(sample=s,n=len(mapping),restored_vs2_records=len(added),baseline_primary=len(old),updated_primary=len(new),gained=len(new-old),lost=len(old-new),newly_gated_parents=len({r['parent_id'] for r in newgate if r['sample']==s})))
        write(out/'summary.tsv',summary);write(out/'changes.tsv',changes);write(out/'newly_gated.tsv',newgate)
        print(summary[-1],flush=True)
    (out/'scope.json').write_text(json.dumps(dict(scope='saved downstream evidence with regenerated discovery-gate verification',newly_gated=len(newgate),baseline_four_tables='byte-identical for every sample',require_unchanged=a.require_unchanged,updated_gate_and_four_tables='byte-identical for every sample' if a.require_unchanged else 'not required'),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--review-root',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);p.add_argument('--archive',type=Path);p.add_argument('--fasta-root',type=Path);p.add_argument('--require-unchanged',action='store_true')
    run(p.parse_args())
