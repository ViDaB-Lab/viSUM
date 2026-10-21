#!/usr/bin/env python3
"""Acquire the exact listed RefSeq assemblies; never select or substitute references.

Uses NCBI HTTPS directory indexes and their published MD5 checksums. The existing
classifier consumes the resulting protein-package/feature-table layout unchanged.
SHA256 receipts freeze acquired bytes, not an unavailable historical annotation.
The caller must hold the destination installation lock.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import shutil
import time
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen

from download_vicat_refseq_metadata import validate_feature_table

NCBI_ROOT = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'
ACCESSION = re.compile(r'GCF_([0-9]{9})\.([1-9][0-9]*)\Z')
COMPONENT = re.compile(r'[A-Za-z0-9_][A-Za-z0-9_.-]*\Z')


def digest(path, algorithm='sha256'):
    result = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.part')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    temporary.replace(path)


def read_inventory(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'assembly_accession', 'cellular_group'} <= set(reader.fieldnames or []):
            raise ValueError('Manifest requires assembly_accession and cellular_group')
        rows = list(reader)
    seen = set()
    for row in rows:
        accession, group = row['assembly_accession'], row['cellular_group']
        if not ACCESSION.fullmatch(accession) or not COMPONENT.fullmatch(group):
            raise ValueError(f'Unsafe/unversioned RefSeq assembly or group: {accession!r}, {group!r}')
        if accession in seen:
            raise ValueError(f'Duplicate assembly: {accession}')
        seen.add(accession)
    if not rows:
        raise ValueError('Assembly manifest is empty')
    return rows


def fetch(url, destination, retries=3, timeout=120):
    """Atomic, bounded download. Never treat a partial file as completed."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + '.part')
    for attempt in range(retries):
        try:
            request = Request(url, headers={'User-Agent': 'viSUM-nonviral-reference-setup/1'})
            with urlopen(request, timeout=timeout) as response, temporary.open('wb') as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
                expected = response.headers.get('Content-Length')
            if expected is not None and temporary.stat().st_size != int(expected):
                raise ValueError(f'Incomplete HTTP response: {url}')
            temporary.replace(destination)
            return
        except (OSError, ValueError):
            temporary.unlink(missing_ok=True)
            if attempt + 1 == retries:
                raise
            time.sleep(min(2 ** attempt, 4))


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls = []

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.urls.extend(value for key, value in attrs if key == 'href' and value)


def resolve_directory(accession, work, retries, timeout):
    digits = ACCESSION.fullmatch(accession).group(1)
    parent = NCBI_ROOT + 'GCF/' + '/'.join(digits[i:i + 3] for i in (0, 3, 6)) + '/'
    listing = work / 'directory.html'
    fetch(parent, listing, retries, timeout)
    links = Links()
    links.feed(listing.read_text(encoding='utf-8'))
    matches = set()
    for link in links.urls:
        url = urljoin(parent, link)
        name = urlparse(url).path.rstrip('/').rsplit('/', 1)[-1]
        if url.startswith(parent) and url.endswith('/') and name.startswith(accession + '_'):
            # Require a direct child, not a nested or external directory.
            if url[len(parent):].count('/') == 1:
                matches.add(url)
    if len(matches) != 1:
        raise ValueError(f'{accession}: expected exactly one matching NCBI directory; found {len(matches)}. No version substitution is allowed.')
    return matches.pop()


def md5_index(path):
    entries = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) == 2 and re.fullmatch(r'[a-fA-F0-9]{32}', parts[0]):
            name = parts[1].lstrip('*').removeprefix('./')
            if name in entries and entries[name] != parts[0].lower():
                raise ValueError('Conflicting NCBI checksum entries: ' + name)
            entries[name] = parts[0].lower()
    return entries


def validate_proteins(path):
    count = 0
    sequence = False
    with Path(path).open(encoding='utf-8') as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith('>'):
                if (count and not sequence) or not text[1:].strip():
                    raise ValueError('Empty FASTA record: ' + str(path))
                count += 1
                sequence = False
            else:
                if not count or not re.fullmatch(r'[A-Za-z*.-]+', text):
                    raise ValueError('Invalid protein FASTA: ' + str(path))
                sequence = True
    if not count or not sequence:
        raise ValueError('Empty/incomplete protein FASTA: ' + str(path))


def acquire(row, root, retries, timeout):
    accession, group = row['assembly_accession'], row['cellular_group']
    work = root / 'downloads' / accession
    receipt_path = work / 'receipt.json'
    protein = root / 'protein_package/ncbi_dataset/data' / accession / 'protein.faa'
    feature = root / 'feature_tables' / group / accession / f'{accession}.feature_table.txt.gz'
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        if receipt['assembly_accession'] != accession or receipt['cellular_group'] != group:
            raise ValueError('Downloaded assembly identity changed: ' + accession)
        # Completed sources are immutable; do not silently repair/replace them.
        for key, path in (('protein_sha256', protein), ('feature_sha256', feature)):
            if not path.is_file() or digest(path) != receipt[key]:
                raise ValueError('Completed downloaded source changed/missing; use a new source directory: ' + str(path))
        return receipt

    directory = resolve_directory(accession, work, retries, timeout)
    basename = urlparse(directory).path.rstrip('/').rsplit('/', 1)[-1]
    index = work / 'md5checksums.txt'
    fetch(directory + 'md5checksums.txt', index, retries, timeout)
    expected = md5_index(index)
    protein_name, feature_name = basename + '_protein.faa.gz', basename + '_feature_table.txt.gz'
    for name in (protein_name, feature_name):
        if name not in expected:
            raise ValueError(f'{accession}: required file/checksum unavailable: {name}')
    compressed = work / 'protein.faa.gz'
    feature_download = work / 'feature_table.txt.gz'
    for name, destination in ((protein_name, compressed), (feature_name, feature_download)):
        # A previous interrupted assembly may already have one complete file.
        if not destination.is_file() or digest(destination, 'md5') != expected[name]:
            fetch(directory + name, destination, retries, timeout)
        if digest(destination, 'md5') != expected[name]:
            raise ValueError(f'{accession}: NCBI MD5 mismatch: {name}')
    validate_feature_table(feature_download)
    # Read the entire gzip stream to detect truncation/CRC failures, not only its header.
    with gzip.open(feature_download, 'rb') as handle:
        for _ in iter(lambda: handle.read(1024 * 1024), b''):
            pass
    protein.parent.mkdir(parents=True, exist_ok=True)
    partial = protein.with_name('protein.faa.part')
    try:
        with gzip.open(compressed, 'rb') as source, partial.open('wb') as output:
            shutil.copyfileobj(source, output, length=1024 * 1024)
        validate_proteins(partial)
        partial.replace(protein)
    finally:
        partial.unlink(missing_ok=True)
    feature.parent.mkdir(parents=True, exist_ok=True)
    partial_feature = feature.with_name(feature.name + '.part')
    shutil.copyfile(feature_download, partial_feature)
    partial_feature.replace(feature)
    receipt = dict(assembly_accession=accession, cellular_group=group,
                   protein_url=directory + protein_name, feature_url=directory + feature_name,
                   protein_download_md5=expected[protein_name], feature_download_md5=expected[feature_name],
                   protein_download_sha256=digest(compressed),
                   protein_sha256=digest(protein), feature_sha256=digest(feature))
    atomic_json(receipt_path, receipt)
    return receipt


def download_sources(manifest, root, workers=4, retries=3, timeout=120):
    manifest, root = Path(manifest), Path(root)
    if workers < 1 or retries < 1 or timeout < 1:
        raise ValueError('workers, retries and timeout must be positive')
    rows = read_inventory(manifest)
    identity = dict(manifest_sha256=digest(manifest), assembly_count=len(rows))
    identity_path = root / 'manifest_identity.json'
    if identity_path.exists() and json.loads(identity_path.read_text()) != identity:
        raise ValueError('Download manifest changed; use a new work/destination directory')
    atomic_json(identity_path, identity)
    records, errors = [], []
    with ThreadPoolExecutor(max_workers=min(workers, 4)) as pool:
        jobs = {pool.submit(acquire, row, root, retries, timeout): row for row in rows}
        for job in as_completed(jobs):
            row = jobs[job]
            try:
                records.append(job.result())
                print('SOURCE_OK ' + row['assembly_accession'], flush=True)
            except Exception as error:
                errors.append(dict(assembly_accession=row['assembly_accession'], error=str(error)))
    atomic_json(root / 'download_summary.json', dict(status='FAIL' if errors else 'PASS',
                requested=len(rows), completed=len(records), errors=sorted(errors, key=lambda r: r['assembly_accession'])))
    if errors:
        raise ValueError(f'{len(errors)} required assemblies failed; no classification/build was started. See {root / "download_summary.json"}')
    provenance = dict(schema='visum-nonviral-download-v1', **identity,
        annotation_identity='acquired bytes recorded below; not a claim of historical annotation equivalence',
        assemblies=sorted(records, key=lambda r: r['assembly_accession']))
    atomic_json(root / 'download_provenance.json', provenance)
    return root / 'protein_package', root / 'feature_tables', root / 'download_provenance.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    download_sources(args.manifest, args.output_dir, args.workers)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f'ERROR: {error}') from error
