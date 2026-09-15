"""Verify expected RNA counts and unchanged DNA scientific outputs after a cached run."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--negative-results',type=Path,required=True)
    p.add_argument('--combined-results',type=Path,required=True)
    p.add_argument('--baseline-combined-results',type=Path,required=True)
    args=p.parse_args()
    for sample,n,expected,root in [
        ('RNANEG_CLEAN_CONFIRMED',1200,20,args.negative_results),
        ('RNANEG_CLEAN_BROAD',400,7,args.negative_results),
        ('RNANEG_MOBILE_CHALLENGE',200,0,args.negative_results),
        ('BMARK_RNA_VIRAL',874,788,args.combined_results)]:
        h=root/f'{sample}_results'/'viharmony'
        manifest=json.loads((h/f'{sample}.harmonizer_manifest.json').read_text())
        assert manifest['policy']['rna_pair_homology_floor'] is True,sample
        for entry in manifest['outputs'].values():
            assert hashlib.sha256((h/entry['file']).read_bytes()).hexdigest()==entry['sha256'],entry['file']
        disp=rows(h/f'{sample}.sequence_disposition.tsv')
        primary=rows(h/f'{sample}.database_candidates.tsv')
        retained={r['original_contig_name'] for r in primary}
        assert len(disp)==n and len({r['original_contig_name'] for r in disp})==n,sample
        assert len(retained)==expected,(sample,len(retained),expected)
        print(f'PASS {sample}: {len(retained)}/{n} primary inputs')
    dna=list(args.baseline_combined_results.glob('BMARK_DNA_*_results'))
    assert dna,'No baseline DNA samples found'
    for old in dna:
        sample=old.name.removesuffix('_results')
        new=args.combined_results/old.name/'viharmony'
        for suffix in ['final_metadata.tsv','sequence_disposition.tsv','final.normalized.fasta',
                       'database_candidates.tsv','provisional_metadata.tsv']:
            name=f'{sample}.{suffix}'
            assert (old/'viharmony'/name).read_bytes()==(new/name).read_bytes(),f'DNA changed: {name}'
        print(f'PASS {sample}: DNA outputs unchanged')

if __name__=='__main__':
    main()
