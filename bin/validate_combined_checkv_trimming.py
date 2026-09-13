#!/usr/bin/env python3
"""Validate combined policies; refuse to reuse taxonomy for changed sequences."""
import argparse,csv,json,subprocess,sys
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor
import refine_proviral_regions as refiner

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def execute(folder,a):
    s=folder.name.removesuffix('_results');typ='rna' if 'RNA' in s else 'dna'
    d=a.output_dir/s;d.mkdir();fa=a.fasta_root/folder.name if a.fasta_root else folder
    cv=a.context_root/s/'checkv_context.tsv'
    ev=[folder/t/f'{s}.{t}_evidence.tsv' for t in ['genomad','virsorter2','cenotetaker3','vicat']]
    ev += [cv,folder/'vicat'/f'{s}.vicat_provirus_evidence.tsv']
    refiner.run(SimpleNamespace(sample_id=s,input_type=typ,candidate_fasta=a.context_root/s/'candidates.fasta',evidence=ev,
        allow_ct3_only_refinement=False,vicat_support_min_overlap_fraction=0.5,output_fasta=d/'refined.fasta',
        output_map=d/'map.tsv',output_audit=d/'audit.tsv',output_summary=d/'summary.tsv'))
    old=read(folder/'refinement'/f'{s}.provirus_region_map.tsv');new=read(d/'map.tsv')
    def signature(rows):
        return {(r['sequence_id'],r['parent_sequence_id'],r['coordinates'],r['record_type'],r['refined_length']) for r in rows}
    result={'sample':s,'refinement_validated':True,'sequence_mapping_unchanged':signature(old)==signature(new)}
    if not result['sequence_mapping_unchanged']:
        result['status']='pending_regenerated_projection_and_taxonomy_on_server'
        return result
    assert refiner.load_fasta(d/'refined.fasta')==refiner.load_fasta(fa/'refinement'/f'{s}.refined_candidates.fasta')
    tools=['genomad','virsorter2','cenotetaker3','deep6','deepmicroclass2','virbot','gianthunter','vicat','tesorter','vitap','vcontact3']
    evidence=[folder/t/f'{s}.{t}_evidence.tsv' for t in tools if (folder/t/f'{s}.{t}_evidence.tsv').exists()]
    evidence += [cv,folder/'vicat'/f'{s}.vicat_refined_evidence.tsv']
    if typ=='rna':evidence.append(folder/'vicat'/f'{s}.vicat_orf_evidence.tsv')
    common=['--sample-id',s,'--input-type',typ,'--normalized-fasta',fa/'prep'/f'{s}.normalized.fasta',
      '--refined-fasta',d/'refined.fasta','--header-map',folder/'prep'/f'{s}.header_map.tsv',
      '--discovery-gate',folder/'discovery_gate'/f'{s}.discovery_gate.tsv','--region-map',d/'map.tsv',
      '--ictv-msl',a.ictv_msl,'--vcontact3-groups',folder/'vcontact3'/f'{s}.vcontact3_group_membership.tsv',
      '--audit-mode','full','--output-prefix',d/s,'--evidence',*evidence]
    if typ=='rna':common.append('--rna-pair-homology-floor')
    r=subprocess.run([sys.executable,str(Path(__file__).with_name('run_viharmony.py')),*map(str,common)],capture_output=True,text=True)
    (d/'run.log').write_text(r.stdout+r.stderr)
    if r.returncode:raise RuntimeError(r.stderr)
    primary=read(d/f'{s}.database_candidates.tsv')
    result.update(status='combined_harmonizer_pass',primary=len({r['normalized_name'] for r in primary}),
      exceptions=[r['original_contig_name'] for r in primary if r['rna_checkv_conflict_exception']=='true'])
    for row in primary:
        if row['rna_checkv_conflict_exception']=='true':
            assert typ=='rna' and row['cellular_conflict_tools']=='checkv' and row['rna_checkv_exception_supporting_loci']!='NA'
    print(json.dumps(result),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['results-root','context-root','ictv-msl','output-dir']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--fasta-root',type=Path)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    with ThreadPoolExecutor(max_workers=2) as pool:r=list(pool.map(lambda f:execute(f,a),sorted(a.results_root.glob('*_results'))))
    (a.output_dir/'summary.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))

if __name__=='__main__':main()
