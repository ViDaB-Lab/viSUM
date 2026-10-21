#!/usr/bin/env python3
"""Audit exact parent identity against saved RNA matches and NCBI summaries.

Offline: never guesses routes for missing records, resolves conflicting sources,
launches tools, or treats a nucleotide DNA annotation as proof of DNA-virus biology.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile


def load_parents(package):
    with tarfile.open(package) as archive:
        names = archive.getnames()
        for name in ('summary.json', 'viral_parent_fingerprints.tsv.gz', 'missing_parents.txt'):
            if names.count(name) != 1:
                raise ValueError(f'Missing or duplicate archive member: {name}')
        summary = json.load(archive.extractfile('summary.json'))
        packed = archive.extractfile('viral_parent_fingerprints.tsv.gz').read()
        if hashlib.sha256(packed).hexdigest() != summary['fingerprints_sha256']:
            raise ValueError('Fingerprint checksum mismatch')
        if summary['status'] != 'PASS' or archive.extractfile('missing_parents.txt').read().strip():
            raise ValueError('Incomplete fingerprint export')
    rows = list(csv.DictReader(io.StringIO(gzip.decompress(packed).decode()), delimiter='\t'))
    parents = {r['parent_sequence_id']: r for r in rows}
    if len(rows) != len(parents) or len(parents) != summary['required_parents']:
        raise ValueError('Duplicate or missing parent rows')
    if sum(int(r['fragment_count']) for r in rows) != summary['viral_fragments']:
        raise ValueError('Fragment-count mismatch')
    for row in rows:
        if int(row['fragment_count']) < 1 or int(row['length_nt']) < 1:
            raise ValueError('Nonpositive fragment count or sequence length')
    return parents, summary


def collect_evidence(parents, rna, cache_paths):
    evidence = defaultdict(list)
    for row in rna['matches']:
        parent = parents.get(row['parent'])
        if parent is None or row['sha256'] != parent['canonical_sequence_sha256'] or int(row['length']) != int(parent['length_nt']):
            raise ValueError('RNA identity match does not reconcile with fingerprints')
        evidence[row['parent']].append(dict(route='rna', source=rna['source_url'],
            record=f"{row['source_id']};exact_full_length_canonical_sha256={row['sha256']};source_sha256={rna['source_sha256']}",
            basis='exact_identity_to_published_RNA_reference'))
    rejected = []
    for path in sorted(cache_paths):
        raw = Path(path).read_bytes()
        cache_hash = hashlib.sha256(raw).hexdigest()
        data = json.loads(raw)['result']
        for uid in data.get('uids', []):
            record = data[uid]
            key = record.get('accessionversion')
            if key not in parents:
                continue
            if int(record.get('slen', 0)) != int(parents[key]['length_nt']):
                rejected.append(dict(parent=key, reason='NCBI_length_mismatch', cache=Path(path).name))
                continue
            route = record.get('moltype', '').lower()
            if route not in ('dna', 'rna'):
                rejected.append(dict(parent=key, reason='NCBI_molecule_type_not_DNA_or_RNA', cache=Path(path).name))
                continue
            item = dict(route=route, source='https://www.ncbi.nlm.nih.gov/nuccore/' + key,
                record=f"{key};moltype={route};length={record['slen']};cache_sha256={cache_hash}",
                basis='exact_accession_version_and_length_with_NCBI_record_molecule_type')
            if item not in evidence[key]:
                evidence[key].append(item)
    return evidence, rejected


def audit(package, rna_matches, cache, out):
    parents, original = load_parents(package)
    rna = json.loads(Path(rna_matches).read_text(encoding='utf-8'))
    evidence, rejected = collect_evidence(parents, rna, Path(cache).glob('*.json'))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    counts, fragments = Counter(), Counter()
    conflicts = []
    with (out/'routing.candidate.tsv').open('w',encoding='utf-8',newline='') as handle, (out/'unresolved.tsv').open('w',encoding='utf-8',newline='') as pending:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        unresolved = csv.writer(pending, delimiter='\t', lineterminator='\n')
        writer.writerow(['parent_sequence_id','input_type','evidence_source','evidence_record'])
        unresolved.writerow(['parent_sequence_id','fragment_count','reason','original_header'])
        for key, parent in sorted(parents.items()):
            items = evidence[key]
            types = {v['route'] for v in items}
            route = next(iter(types)) if len(types) == 1 else 'unresolved'
            if route == 'unresolved':
                reason = 'conflicting_source_annotations' if types else 'no_exact_source_match'
                unresolved.writerow([key,parent['fragment_count'],reason,parent['original_header']])
                if types: conflicts.append(dict(parent=key,evidence=items))
            else:
                writer.writerow([key,route,items[0]['source'],items[0]['record']])
            counts[route] += 1
            fragments[route] += int(parent['fragment_count'])
    report = dict(status='CANDIDATE_ROUTING_REQUIRES_REVIEW',launch_ready=False,
        parent_counts=dict(counts),fragment_counts=dict(fragments),conflict_count=len(conflicts),
        conflicts=conflicts,rejected_metadata=rejected,
        truth_sha256=original['truth_sha256'],
        note='NCBI DNA molecule annotations can describe cDNA or integrated copies of RNA viruses. Candidate assignments are not independently validated viral genome types. Unmatched parents are never defaulted to DNA.')
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    with (out/'evidence.jsonl').open('w',encoding='utf-8') as handle:
        for key in sorted(evidence):
            if evidence[key]: handle.write(json.dumps(dict(parent=key,evidence=evidence[key]))+'\n')
    with (out/'checksums.sha256').open('w',encoding='utf-8') as handle:
        for path in sorted(out.iterdir()):
            if path.name != 'checksums.sha256':
                handle.write(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('conflicts','rejected_metadata')},indent=2))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fingerprints',type=Path,required=True)
    p.add_argument('--rna-matches',type=Path,required=True)
    p.add_argument('--ncbi-cache',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a = p.parse_args()
    audit(a.fingerprints,a.rna_matches,a.ncbi_cache,a.out)


if __name__ == '__main__': main()
