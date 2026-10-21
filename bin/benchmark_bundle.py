"""Manifest-driven benchmark access; never modifies calls, labels or sequences."""
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath

PRIMARY = frozenset({'retained_viral', 'retained_provirus'})
PROVISIONAL = frozenset({'provisional_viral', 'provisional_provirus'})
ENDPOINTS = ('primary', 'viral_candidates', 'TE_candidates', 'TE_only')
SUFFIX_LABELS = {
    'final_metadata.tsv': 'metadata', 'sequence_disposition.tsv': 'disposition',
    'sequence_map.tsv': 'map', 'database_candidates.tsv': 'database_tsv',
    'database_candidates.fasta': 'database_fasta', 'final.normalized.fasta': 'normalized_fasta',
    'provisional_metadata.tsv': 'provisional_metadata',
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class Bundle:
    """Resolve BOTH historical flat and current nested layouts from their manifest.

    Missing archived artifacts are allowed only explicitly, reported in integrity,
    and still fail if a caller actually needs them. Never silently choose between
    a stale flat manifest and a new audit manifest in the same directory.
    """
    def __init__(self, directory, sample, allow_missing=False):
        self.root = Path(directory).resolve()
        self.sample = sample
        locations = [self.root / (sample + '.harmonizer_manifest.json'),
                     self.root / 'audit' / (sample + '.harmonizer_manifest.json')]
        found = [p for p in locations if p.is_file()]
        if len(found) != 1:
            raise ValueError('Require exactly one current manifest (missing or mixed layouts): ' + str(self.root))
        self.manifest_path = found[0]
        self.manifest = json.loads(found[0].read_text(encoding='utf-8'))
        if self.manifest.get('sample_id') != sample:
            raise ValueError('Manifest sample mismatch')
        if found[0].parent != self.root and self.manifest.get('reporting_bundle', {}).get('output_root') != '..':
            raise ValueError('Unsupported nested manifest output root')
        self.integrity = []
        names = set()
        for label, record in self.manifest['outputs'].items():
            name = record['file']
            if name in names:
                raise ValueError('Duplicate manifest output path')
            names.add(name)
            path = self.path(label, required=False)
            if path.is_file():
                if sha(path) != record['sha256']:
                    raise ValueError('Output checksum mismatch: ' + str(path))
                status = 'hash_match'
            elif allow_missing:
                status = 'not_extracted_or_not_packaged'
            else:
                raise FileNotFoundError(path)
            self.integrity.append(dict(sample=sample, file=name, status=status))

    def path(self, label, required=True):
        record = self.manifest['outputs'].get(label)
        if record is None:
            if required: raise ValueError('Missing manifest artifact label: ' + label)
            return None
        name = record['file']
        relative = PurePosixPath(name)
        if relative.is_absolute() or '..' in relative.parts or '\\' in name or ':' in name:
            raise ValueError('Unsafe manifest output path: ' + name)
        path = (self.root / name).resolve()
        if self.root not in path.parents or not path.name.startswith(self.sample + '.'):
            raise ValueError('Escaped/foreign manifest output: ' + name)
        if required and not path.is_file(): raise FileNotFoundError(path)
        return path

    def table(self, label):
        return read(self.path(label))

    def endpoints(self):
        disposition = self.table('disposition')
        mapping = {r['normalized_name']: r['original_contig_name'] for r in disposition}
        if len(mapping) != len(disposition) or len(set(mapping.values())) != len(mapping):
            raise ValueError('Duplicate normalized or original input IDs')
        metadata = self.table('metadata')
        primary_rows = self.table('database_tsv')
        seen = set()
        selected = {name: set() for name in ENDPOINTS}
        primary_ids = set()
        for row in metadata:
            final, parent = row['final_sequence_id'], row['normalized_name']
            if final in seen or parent not in mapping:
                raise ValueError('Duplicate candidate or unknown input parent')
            seen.add(final)
            if row.get('original_contig_name', mapping[parent]) != mapping[parent]:
                raise ValueError('Candidate original ID disagrees with input disposition')
            if 'candidate_destination' in row:
                destination, decision = row['candidate_destination'], row['viral_decision']
                if destination == 'viral_candidates':
                    if decision not in PRIMARY | PROVISIONAL:
                        raise ValueError('Viral destination has incompatible decision')
                    selected['viral_candidates'].add(parent)
                    if decision in PRIMARY:
                        selected['primary'].add(parent); primary_ids.add(final)
                elif destination == 'TE_candidates':
                    selected['TE_candidates'].add(parent)
                elif destination not in ('excluded', 'review_only'):
                    raise ValueError('Unknown candidate destination: ' + destination)
            else:
                # Historical policies predate the export/routing column. Their
                # primary table remains authoritative; no TE routing is invented.
                if row.get('viral_decision') in PRIMARY | PROVISIONAL:
                    selected['viral_candidates'].add(parent)
        table_ids = {r['final_sequence_id'] for r in primary_rows}
        if len(table_ids) != len(primary_rows) or not table_ids <= seen:
            raise ValueError('Primary table has duplicate/unknown candidate IDs')
        indexed = {r['final_sequence_id']: r for r in metadata}
        if any(row['normalized_name'] != indexed[row['final_sequence_id']]['normalized_name'] for row in primary_rows):
            raise ValueError('Primary artifact parent differs from candidate metadata')
        if metadata and 'candidate_destination' in metadata[0]:
            if table_ids != primary_ids:
                raise ValueError('Primary artifact disagrees with retained viral-destination metadata')
        else:
            selected['primary'] = {r['normalized_name'] for r in primary_rows}
            selected['viral_candidates'].update(selected['primary'])
        if not selected['primary'] <= selected['viral_candidates'] <= set(mapping):
            raise ValueError('Endpoint/input sets do not reconcile')
        # A parent with a viral child AND TE child is counted in both export
        # endpoints, but never in TE_only. Regions never increase denominators.
        selected['TE_only'] = selected['TE_candidates'] - selected['viral_candidates']
        return mapping, selected


def artifact(directory, sample, suffix):
    return Bundle(directory, sample, allow_missing=True).path(SUFFIX_LABELS[suffix])
