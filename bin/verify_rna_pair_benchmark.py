"""Verify expected RNA counts and unchanged DNA scientific outputs after a cached run."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from benchmark_bundle import Bundle, SUFFIX_LABELS

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
        bundle=Bundle(h,sample)
        manifest=bundle.manifest
        assert manifest['policy']['rna_pair_homology_floor'] is True,sample
        disp=bundle.table('disposition')
        primary=bundle.table('database_tsv')
        retained={r['original_contig_name'] for r in primary}
        assert len(disp)==n and len({r['original_contig_name'] for r in disp})==n,sample
        assert len(retained)==expected,(sample,len(retained),expected)
        print(f'PASS {sample}: {len(retained)}/{n} primary inputs')
    dna=list(args.baseline_combined_results.glob('BMARK_DNA_*_results'))
    assert dna,'No baseline DNA samples found'
    for old in dna:
        sample=old.name.removesuffix('_results')
        new=args.combined_results/old.name/'viharmony'
        previous=Bundle(old/'viharmony',sample)
        current=Bundle(new,sample)
        for suffix in ['final_metadata.tsv','sequence_disposition.tsv','final.normalized.fasta',
                       'database_candidates.tsv','provisional_metadata.tsv']:
            name=f'{sample}.{suffix}'
            label=SUFFIX_LABELS[suffix]
            assert previous.path(label).read_bytes()==current.path(label).read_bytes(),f'DNA changed: {name}'
        print(f'PASS {sample}: DNA outputs unchanged')

if __name__=='__main__':
    main()
