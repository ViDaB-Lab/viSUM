#!/usr/bin/env python3
"""Prepare fingerprinted NCBI inputs, then invoke the existing nonviral builder.

The calling Nextflow process owns the destination lock. Optional acquisition uses
the supplied versioned inventory without changing classification or clustering.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

from classify_vicat_nonviral_references import read_manifest
from download_vicat_nonviral_sources import download_sources


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def source_inventory(manifest, package_root, metadata_root):
    files = {'manifest': Path(manifest)}
    for row in read_manifest(Path(manifest)):
        accession = row['assembly_accession'].strip()
        group = row['cellular_group'].strip()
        for value in (accession, group):
            if not value or value in ('.', '..') or any(c in value for c in '/\\'):
                raise ValueError('Manifest accessions and groups must be single path components')
        files[f'{accession}/proteins'] = Path(package_root) / 'ncbi_dataset' / 'data' / accession / 'protein.faa'
        files[f'{accession}/features'] = Path(metadata_root) / group / accession / f'{accession}.feature_table.txt.gz'
    return {name: sha256(path) for name, path in sorted(files.items())}


def classified_inputs(args):
    if args.classified_dir:
        return Path(args.classified_dir).resolve()
    work = Path(args.work_dir).resolve()
    classified = work / 'classified'
    stamp = work / 'classification_provenance.json'
    script = Path(__file__).with_name('classify_vicat_nonviral_references.py')
    inventory = source_inventory(args.manifest, args.package_root, args.metadata_root)
    if getattr(args, 'download_provenance', None):
        inventory['download_sources'] = json.loads(Path(args.download_provenance).read_text())
    inventory['classifier_sha256'] = sha256(script)
    if stamp.exists():
        recorded = json.loads(stamp.read_text())
        if recorded['inputs'] != inventory:
            raise ValueError('Nonviral classification inputs changed; use a new work/destination directory')
        for name, digest in recorded['outputs'].items():
            if sha256(classified / name) != digest:
                raise ValueError(f'Classified output changed: {name}')
        return classified
    # The classifier safely replaces its own partial output files. No builder
    # is invoked until classification and its provenance record are complete.
    work.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, str(script), '--manifest', args.manifest,
                    '--package-root', args.package_root, '--metadata-root', args.metadata_root,
                    '--output-dir', str(classified)], check=True)
    shutil.copyfile(args.manifest, classified / 'vicat_nonviral_source_manifest.tsv')
    (classified / 'vicat_nonviral_source_provenance.json').write_text(
        json.dumps(inventory, indent=2) + '\n')
    outputs = {p.name: sha256(p) for p in sorted(classified.iterdir()) if p.is_file()}
    temporary = stamp.with_suffix('.tmp')
    temporary.write_text(json.dumps({'inputs': inventory, 'outputs': outputs}, indent=2) + '\n')
    temporary.replace(stamp)
    return classified


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--classified-dir')
    parser.add_argument('--manifest')
    parser.add_argument('--package-root')
    parser.add_argument('--metadata-root')
    parser.add_argument('--download-sources', action='store_true')
    parser.add_argument('--destination', required=True)
    parser.add_argument('--work-dir', required=True)
    parser.add_argument('--threads', required=True, type=int)
    args = parser.parse_args()
    raw = [args.manifest, args.package_root, args.metadata_root]
    valid_sources = (bool(args.classified_dir) and not any(raw) and not args.download_sources) or (
        not args.classified_dir and bool(args.manifest) and (
            (args.download_sources and not args.package_root and not args.metadata_root) or
            (not args.download_sources and bool(args.package_root) and bool(args.metadata_root))))
    if args.threads < 1 or not valid_sources:
        parser.error('Supply classified-dir OR all three NCBI source arguments OR manifest with download-sources, and positive threads')
    if Path(args.destination).exists():
        parser.error('Refusing to replace an existing database')
    if args.download_sources:
        args.package_root, args.metadata_root, args.download_provenance = map(str, download_sources(
            args.manifest, Path(args.work_dir).resolve() / 'sources', workers=min(args.threads, 4)))
    classified = classified_inputs(args)
    subprocess.run(['bash', str(Path(__file__).with_name('build_vicat_nonviral_database.sh')),
                    '--classified-dir', str(classified), '--destination', args.destination,
                    '--work-dir', str(Path(args.work_dir).resolve() / 'clustering'),
                    '--threads', str(args.threads)], check=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(f'ERROR: {error}') from error
