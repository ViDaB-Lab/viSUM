#!/usr/bin/env python3
"""Offline structural validation for geNomad 1.12 and vConTACT3 3.2.4 assets.

Contracts follow genomad/database.py (v1.12.0) and vcontact3 databases.py,
profiles.py, config.py (d57ae81e0be5ab57a8a14692a4cb8ab255ce2d52).
This checks completeness/readability, not cryptographic authenticity.
"""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess


def nonempty(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing or empty required file: {path}')
    with path.open('rb') as handle:
        handle.read(1)
    return path


def mmseqs(prefix, runtime=False, headers=False):
    prefix = Path(prefix)
    dbtype = nonempty(str(prefix) + '.dbtype')
    if dbtype.stat().st_size != 4:
        raise ValueError(f'Invalid MMseqs dbtype header: {dbtype}')
    nonempty(str(prefix) + '.index')
    if prefix.is_file() and prefix.stat().st_size:
        nonempty(prefix)
    else:
        parts = sorted((p for p in prefix.parent.glob(prefix.name + '.*')
                        if p.suffix[1:].isdigit()), key=lambda p: int(p.suffix[1:]))
        if not parts or [int(p.suffix[1:]) for p in parts] != list(range(len(parts))):
            raise ValueError(f'Missing MMseqs data or noncontiguous split files: {prefix}')
        for part in parts:
            nonempty(part)
    if runtime:
        subprocess.run(['mmseqs', 'dbtype', str(prefix)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if headers:
        mmseqs(str(prefix) + '_h', runtime)


def genomad(directory, runtime=False):
    root = Path(directory)
    version = float(nonempty(root / 'version.txt').read_text().strip())
    if not math.isfinite(version) or version <= 0:
        raise ValueError('Invalid geNomad database version')
    for name in ('genomad_db', 'genomad_mini_db', 'genomad_integrase_db'):
        mmseqs(root / name, runtime, headers=True)
    for name in ('nodes.dmp', 'names.dmp', 'genomad_marker_metadata.tsv'):
        nonempty(root / name)
    with (root / 'genomad_marker_metadata.tsv').open() as handle:
        header, first = handle.readline(), handle.readline()
        if '\t' not in header or '\t' not in first:
            raise ValueError('geNomad marker metadata is empty or malformed')


def local_asset(root, release, value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Reference asset must be a nonempty relative path')
    path = Path(value)
    if path.is_absolute() or '..' in path.parts or '\\' in value or ':' in value:
        raise ValueError(f'Unsafe/nonportable manifest path: {value}')
    resolved = (root / path).resolve()
    if not resolved.is_relative_to(release.resolve()):
        raise ValueError(f'Manifest asset is outside its release: {value}')
    return resolved


def vcontact3(manifest, domain='both', runtime=False):
    manifest = nonempty(manifest)
    if not re.fullmatch(r'[0-9]{3}', manifest.stem):
        raise ValueError('Expected a three-digit vConTACT3 release manifest')
    data = json.loads(manifest.read_text())
    release = manifest.parent / ('v' + manifest.stem)
    if not release.is_dir():
        raise ValueError(f'Missing release directory: {release}')
    info = data[manifest.stem]
    refseq = info['RefSeq']
    domains = ['prokaryotes', 'eukaryotes'] if domain == 'both' else [domain]
    for selected in domains:
        block = refseq['domains'][selected]
        for key in ('Genomes Report', 'gene2genome', 'proteins', 'nucleotides'):
            nonempty(local_asset(manifest.parent, release, block[key]))
        identities = block['identities']
        for identity in ('0.3', '0.4', '0.5', '0.6', '0.7'):
            cluster = identities[identity]
            for key in ('db', 'cluster'):
                mmseqs(local_asset(manifest.parent, release, cluster[key]), runtime,
                       headers=(key == 'db'))
            nonempty(local_asset(manifest.parent, release, cluster['cluster2genome']))
            # Also check other file references declared by this identity block.
            for key, value in cluster.items():
                if key not in ('db', 'cluster', 'cluster2genome') and isinstance(value, str) and '/' in value:
                    nonempty(local_asset(manifest.parent, release, value))
    mmseqs(local_asset(manifest.parent, release, info['VOGDB']['mmseqs_db']), runtime, headers=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tool', choices=['genomad', 'vcontact3'])
    parser.add_argument('path')
    parser.add_argument('--domain', choices=['both', 'prokaryotes', 'eukaryotes'], default='both')
    parser.add_argument('--check-mmseqs', action='store_true')
    args = parser.parse_args()
    try:
        if args.tool == 'genomad':
            genomad(args.path, args.check_mmseqs)
        else:
            vcontact3(args.path, args.domain, args.check_mmseqs)
    except (OSError, ValueError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        raise SystemExit(f'ERROR: {args.tool} database validation failed: {error}') from error


if __name__ == '__main__':
    main()
