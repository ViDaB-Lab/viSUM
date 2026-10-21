#!/usr/bin/env python3
"""Join exact viral parent IDs to saved NCBI GenBank XML molecule metadata.

No accession-prefix guessing, taxonomy inference, network access or cluster propagation.
Non-accession metagenomic identifiers remain explicitly unresolved.
"""
import argparse
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from prepare_genomad_benchmark import sha256, table_writer


def molecule_route(value):
    value = value.upper()
    dna, rna = 'DNA' in value, 'RNA' in value
    return 'dna' if dna and not rna else 'rna' if rna and not dna else None


def resolve(required, xml_paths, out, source):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    with Path(required).open(encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'parent_sequence_id' not in (reader.fieldnames or []):
            raise ValueError('Required table lacks parent_sequence_id')
        keys = [row['parent_sequence_id'] for row in reader]
    if len(set(keys)) != len(keys) or any(not key for key in keys):
        raise ValueError('Blank/duplicate required parent ID')
    wanted, found = set(keys), {}
    hashes = {}
    for path in xml_paths:
        path = Path(path)
        hashes[str(path.resolve())] = sha256(path)
        for event, node in ET.iterparse(path, events=('end',)):
            if node.tag not in {'GBSeq', 'INSDSeq'}:
                continue
            prefix = node.tag
            accession = node.findtext(f'{prefix}_accession-version', '').strip()
            primary = node.findtext(f'{prefix}_primary-accession', '').strip()
            molecule = node.findtext(f'{prefix}_moltype', '').strip()
            route = molecule_route(molecule)
            for key in wanted.intersection({accession, primary}):
                if not accession or route is None:
                    continue
                value = (route, source, f'{accession};moltype={molecule};xml_sha256={hashes[str(path.resolve())]}')
                if key in found and found[key] != value:
                    raise ValueError(f'Conflicting/duplicate metadata for {key}; use one frozen record version')
                found[key] = value
            node.clear()
    with (out / 'routing.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = table_writer(f, ['parent_sequence_id', 'input_type', 'evidence_source', 'evidence_record'])
        for key in sorted(found):
            writer.writerow([key, *found[key]])
    with (out / 'unresolved.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = table_writer(f, ['parent_sequence_id', 'reason'])
        writer.writerows((key, 'No exact accession match with unambiguous molecule metadata')
                         for key in sorted(wanted - found.keys()))
    report = {'required_parents': len(wanted), 'resolved_parents': len(found),
              'unresolved_parents': len(wanted) - len(found), 'xml_sha256': hashes,
              'required_sha256': sha256(required), 'source': source}
    (out / 'routing_summary.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--required', type=Path, required=True)
    p.add_argument('--genbank-xml', type=Path, nargs='+', required=True)
    p.add_argument('--source', required=True, help='Metadata source URL or versioned citation')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not a.source.strip():
        p.error('--source must not be blank')
    resolve(a.required, a.genbank_xml, a.out, a.source)


if __name__ == '__main__':
    main()
