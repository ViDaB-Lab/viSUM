#!/usr/bin/env python3
"""Freeze reviewed RNA routes; defer ambiguous DNA record annotations to dual routing.

This conservative policy uses unopposed RNA source evidence. NCBI DNA annotations
alone do not establish a DNA viral genome and remain in the unresolved cohort.
No classifier predictions or benchmark outcomes are used to select a route.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from audit_genomad_routing_evidence import load_parents
from prepare_genomad_benchmark import sha256, table_writer


def finalize(package, evidence_path, ncbi_path, out):
    parents, source = load_parents(package)
    ncbi = json.loads(Path(ncbi_path).read_text(encoding='utf-8'))
    observed = {}
    with Path(evidence_path).open(encoding='utf-8') as handle:
        for line in handle:
            row = json.loads(line)
            key = row['parent']
            if key not in parents or key in observed:
                raise ValueError('Unknown/duplicate evidence parent')
            items = row['evidence']
            if not items or any(v.get('route') not in ('dna','rna') or not v.get('source') or not v.get('record') for v in items):
                raise ValueError('Invalid source evidence')
            for item in items:
                if item.get('basis') == 'exact_identity_to_published_RNA_reference':
                    if f"exact_full_length_canonical_sha256={parents[key]['canonical_sequence_sha256']}" not in item['record']:
                        raise ValueError('RNA sequence hash mismatch')
                elif item.get('basis') != 'exact_accession_version_and_length_with_NCBI_record_molecule_type':
                    raise ValueError('Unreviewed evidence basis')
            observed[key] = items
    out = Path(out)
    out.mkdir(parents=True,exist_ok=False)
    counts,fragments = Counter(),Counter()
    with (out/'routing.reviewed.tsv').open('w',encoding='utf-8',newline='') as f, (out/'routing_decisions.tsv').open('w',encoding='utf-8',newline='') as d:
        writer=table_writer(f,['parent_sequence_id','input_type','evidence_source','evidence_record'])
        decisions=table_writer(d,['parent_sequence_id','route','reason','fragment_count'])
        for key,parent in sorted(parents.items()):
            items=observed.get(key,[])
            types={v['route'] for v in items}
            exact_rna = any(v['basis']=='exact_identity_to_published_RNA_reference' for v in items)
            record = ncbi.get(key,{}).get('record',{})
            genomic_rna = (record.get('accessionversion') == key and record.get('moltype') == 'rna'
                           and record.get('biomol') in ('genomic','cRNA')
                           and int(record.get('slen',0)) == int(parent['length_nt']))
            if types == {'rna'} and (exact_rna or genomic_rna):
                route,reason='rna','unopposed_RNA_source_evidence'
                chosen=items[0]
                writer.writerow([key,route,chosen['source'],chosen['record']])
            else:
                route='unresolved'
                reason=('conflicting_source_annotations' if len(types)>1 else
                        'DNA_record_not_independent_viral_genome_type' if types == {'dna'} else
                        'RNA_record_without_genome_type_support' if types else 'no_exact_source_match')
            decisions.writerow([key,route,reason,parent['fragment_count']])
            counts[route]+=1
            fragments[route]+=int(parent['fragment_count'])
    report=dict(status='ROUTING_REVIEWED_FOR_EXPLICIT_DUAL_POLICY',unresolved_policy='both',
        parent_counts=dict(counts),fragment_counts=dict(fragments),truth_sha256=source['truth_sha256'],
        fingerprints_package_sha256=sha256(package),evidence_sha256=sha256(evidence_path),
        ncbi_records_sha256=sha256(ncbi_path),
        routing_sha256=sha256(out/'routing.reviewed.tsv'),script_sha256=sha256(__file__),
        note='Unresolved is not a DNA assignment. Its two workflow results must remain separate.')
    (out/'routing_policy.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fingerprints',type=Path,required=True)
    p.add_argument('--evidence',type=Path,required=True)
    p.add_argument('--ncbi-records',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    finalize(a.fingerprints,a.evidence,a.ncbi_records,a.out)
