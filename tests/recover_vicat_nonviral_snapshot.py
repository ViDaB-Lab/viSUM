"""Offline recovery of the September reference inventory; not pipeline code.

Reads evidence without modifying it. Emits a JSON mapping of relative filenames
to UTF-8 contents for review/application with the normal file-editing tools.
No downloads, database builds, or analysis decisions occur here.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tarfile


SOURCE_SHA256 = 'dbc7cf8fffecad85e65d4670300691f2cc30f3a4e9ac9848c73e9640429dc9c8'
HISTORICAL_COMMIT = 'd8805bd8e83a38aff0c3d97f2b2306a4b74f97f1'
GROUPS = {'bacteria': 300, 'archaea': 80, 'fungi': 30, 'protozoa': 20,
          'plant': 20, 'invertebrate': 20, 'vertebrate_mammalian': 10}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def recover_catalogue(path):
    assemblies, classes, groups = {}, Counter(), Counter()
    notes = Counter()
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'source_accession', 'cellular_group', 'organism_name', 'taxid',
                    'reference_class', 'classification_note'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Source catalog is missing required columns')
        for row in reader:
            accession = row['source_accession']
            if not re.fullmatch(r'GCF_\d+\.\d+', accession):
                raise ValueError('Expected a versioned RefSeq assembly accession')
            annotation = tuple(row[key] for key in ('cellular_group', 'organism_name', 'taxid'))
            if not all(annotation) or not annotation[2].isdigit():
                raise ValueError('Missing or invalid assembly annotation')
            previous = assemblies.setdefault(accession, annotation)
            if previous != annotation:
                raise ValueError('Conflicting annotations for ' + accession)
            classes[row['reference_class']] += 1
            groups[row['cellular_group']] += 1
            notes[row['classification_note']] += 1
    manifest = io.StringIO(newline='')
    writer = csv.writer(manifest, delimiter='\t', lineterminator='\n')
    writer.writerow(('assembly_accession', 'cellular_group', 'organism_name', 'taxid'))
    writer.writerows((accession, *assemblies[accession]) for accession in sorted(assemblies))
    return manifest.getvalue(), {
        'assembly_count': len(assemblies), 'included_protein_count': sum(classes.values()),
        'assemblies_by_group': dict(Counter(a[0] for a in assemblies.values())),
        'proteins_by_class': dict(classes), 'proteins_by_group': dict(groups),
        'classification_notes': dict(notes),
    }


def summary_metrics(text):
    rows = list(csv.DictReader(io.StringIO(text), delimiter='\t'))
    result = {(row['metric'], row['category']): int(row['value']) for row in rows}
    if len(result) != len(rows):
        raise ValueError('Duplicate summary keys')
    return result


def recover(source, review, unmatched):
    if sha256(source) != SOURCE_SHA256:
        raise ValueError('Source catalog checksum does not match the frozen database')
    wanted = ('SHA256SUMS', 'vicat_nonviral_build_config.tsv',
              'vicat_nonviral_database_metadata.tsv', 'vicat_nonviral_cluster_summary.tsv',
              'vicat_nonviral_classification_summary.tsv', 'diamond_dbinfo.txt',
              'representative_audit.tsv', '.visum_db_complete')
    evidence = {}
    with tarfile.open(review, 'r:gz') as archive:
        for member in archive:
            name = Path(member.name).name
            if member.isfile() and name in wanted:
                if name in evidence:
                    raise ValueError('Duplicate archive member: ' + name)
                evidence[name] = archive.extractfile(member).read()
    if set(evidence) != set(wanted):
        raise ValueError('Incomplete build review')
    checksums = {line.split(maxsplit=1)[1]: line.split()[0]
                 for line in evidence['SHA256SUMS'].decode().splitlines() if line.strip()}
    for name, data in evidence.items():
        if name in checksums and hashlib.sha256(data).hexdigest() != checksums[name]:
            raise ValueError('Build evidence checksum mismatch: ' + name)
    if checksums['vicat_nonviral_source_metadata.tsv.gz'] != SOURCE_SHA256:
        raise ValueError('Build review belongs to a different source catalog')
    manifest, stats = recover_catalogue(source)
    matching = summary_metrics(evidence['vicat_nonviral_classification_summary.tsv'].decode())
    if (stats['assembly_count'] != 480 or stats['included_protein_count'] != 3731712
            or stats['assemblies_by_group'] != GROUPS
            or matching[('excluded_protein_count', 'all')] != 0):
        raise ValueError('Unexpected frozen inventory totals')
    for metric in ('assembly_count', 'included_protein_count'):
        if stats[metric] != matching[(metric, 'all')]:
            raise ValueError('Catalog and build summary disagree: ' + metric)
    for metric, values in (('reference_class_count', stats['proteins_by_class']),
                           ('cellular_group_included_count', stats['proteins_by_group'])):
        expected = {category: value for (key, category), value in matching.items() if key == metric}
        if values != expected:
            raise ValueError('Catalog and build summary disagree: ' + metric)
    build = dict(csv.reader(io.StringIO(evidence['vicat_nonviral_build_config.tsv'].decode()), delimiter='\t'))
    root = Path(__file__).resolve().parents[1]
    scripts = {'bin/build_vicat_nonviral_database.sh': build['builder_sha256'],
               'bin/prepare_vicat_nonviral_cluster_metadata.py': build['metadata_helper_sha256']}
    for name, expected in scripts.items():
        data = subprocess.check_output(['git', 'show', HISTORICAL_COMMIT + ':' + name], cwd=root)
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError('Historical script does not match recorded checksum: ' + name)
    unmatched_bytes = Path(unmatched).read_bytes()
    other = summary_metrics(unmatched_bytes.decode())
    if other[('excluded_protein_count', 'all')] != 404 or other == matching:
        raise ValueError('Unexpected unmatched summary')
    files = {'assembly_manifest.tsv': manifest}
    files.update({'build_evidence/' + name: data.decode('utf-8') for name, data in evidence.items()})
    files['unmatched_evidence/classification_summary.tsv'] = unmatched_bytes.decode('utf-8')
    provenance = {
        'snapshot_id': 'nonviral_classaware_v1_20260904',
        'status': 'historical_build_documented_not_a_new_database_release',
        'source_catalog': {'filename': Path(source).name, 'sha256': SOURCE_SHA256},
        'source_review': {'filename': Path(review).name, 'sha256': sha256(review)},
        'manifest': {'filename': 'assembly_manifest.tsv',
                     'sha256': hashlib.sha256(manifest.encode()).hexdigest(),
                     'recovery': 'Unique source_accession and consistent group/organism/taxid from every catalog row; accession-sorted.'},
        'inventory': stats,
        'historical_scripts': {'matching_git_commit': HISTORICAL_COMMIT, 'sha256': scripts,
                               'scope': 'These two files match; this does not identify the entire runtime checkout.'},
        'build': {'completed': '2026-09-04T17:25:40-05:00', 'mmseqs_version': '18.8cc5c',
                  'diamond_version': '2.2.5.185', 'representatives': 2708478,
                  'min_seq_id': 0.90, 'coverage': 0.80, 'cov_mode': 0,
                  'minimum_protein_length': None, 'recorded_threads': 52},
        'database_files_sha256': checksums,
        'verified_here': ['Complete source-catalog hash and row-level inventory/counts',
                          'Hashes of four build metadata/summary files present in the review',
                          'Historical builder and metadata-helper hashes against Git'],
        'not_verified_here': ['Current server database bytes', 'Absent DIAMOND/FASTA/parquet/membership files'],
        'unmatched_summary': {'filename': 'unmatched_evidence/classification_summary.tsv',
                              'source_filename': Path(unmatched).name,
                              'sha256': hashlib.sha256(unmatched_bytes).hexdigest(),
                              'status': 'unmatched_not_a_build_input_or_authoritative_summary',
                              'included_proteins': 3731308, 'excluded_proteins': 404,
                              'chronology': 'unknown; do not infer which classification state came first'},
        'unrecovered_provenance': ['Original assembly-selection criteria', 'Original download dates and commands',
                                   'Original protein FASTA and feature-table hashes/annotation revisions'],
        'redistribution_status': 'Not cleared by this audit; retain accession inventory and review applicable notices before hosting.',
    }
    files['provenance.json'] = json.dumps(provenance, indent=2, sort_keys=True) + '\n'
    return files


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--review', required=True, type=Path)
    parser.add_argument('--unmatched-summary', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(recover(args.source, args.review, args.unmatched_summary)))
