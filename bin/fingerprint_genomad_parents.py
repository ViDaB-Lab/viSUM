#!/usr/bin/env python3
"""Export exact viral-parent fingerprints from the downloaded author ZIP.

No classifiers are run, no routes are guessed, and no existing files are changed.
The small review archive replaces the need to upload the 8.6-GB source archive.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import zipfile

MEMBER = 'reference_sequences/reference_sequences.fna.gz'
COMPLEMENT = str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN')


def digest(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def fingerprints(sequence):
    normalized = sequence.upper().replace('U', 'T')
    if not normalized or set(normalized) - set('ACGTRYMKBDHVNSW'):
        raise ValueError('Empty sequence or invalid nucleotide alphabet')
    reverse = normalized.translate(COMPLEMENT)[::-1]
    return digest(sequence), min(digest(normalized), digest(reverse))


def required_parents(path):
    parents = collections.Counter()
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'sequence_id', 'source_class'} <= set(reader.fieldnames or []):
            raise ValueError('Missing required truth columns')
        for row in reader:
            if row['source_class'] != 'virus':
                continue
            parent, separator, number = row['sequence_id'].rpartition('|fragment_')
            if not parent or not separator or not number.isdigit():
                raise ValueError(f"Malformed fragment ID: {row['sequence_id']}")
            parents[parent] += 1
    if not parents:
        raise ValueError('No viral parents in truth metadata')
    return parents


def selected_records(handle, wanted):
    header = None
    key = None
    pieces = []
    for line in handle:
        line = line.rstrip('\r\n')
        if line.startswith('>'):
            if key in wanted:
                yield key, header, ''.join(pieces)
            header = line[1:]
            if not header.strip():
                raise ValueError('Empty FASTA header')
            key = header.split()[0]
            pieces = []
        elif line:
            if header is None:
                raise ValueError('Sequence before FASTA header')
            if key in wanted:
                pieces.append(line.strip())
    if key in wanted:
        yield key, header, ''.join(pieces)


def run(truth, archive_path, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    parents = required_parents(truth)
    print(f'Fingerprinting {len(parents):,} viral parents from the existing author archive.', flush=True)
    table = out / 'viral_parent_fingerprints.tsv.gz'
    seen = set()
    with zipfile.ZipFile(archive_path) as archive:
        if archive.namelist().count(MEMBER) != 1:
            raise ValueError('Missing or duplicated original reference FASTA member')
        info = archive.getinfo(MEMBER)
        with archive.open(info) as packed, gzip.GzipFile(fileobj=packed) as fasta:
            with io.TextIOWrapper(fasta, encoding='ascii') as handle, gzip.open(table, 'wt', encoding='utf-8', newline='') as output:
                writer = csv.writer(output, delimiter='\t', lineterminator='\n')
                writer.writerow(['parent_sequence_id', 'fragment_count', 'length_nt',
                                 'sequence_sha256', 'canonical_sequence_sha256', 'original_header'])
                for key, header, sequence in selected_records(handle, parents):
                    if key in seen:
                        raise ValueError(f'Duplicate required reference ID: {key}')
                    raw, canonical = fingerprints(sequence)
                    writer.writerow([key, parents[key], len(sequence), raw, canonical, header])
                    seen.add(key)
                    if len(seen) % 10000 == 0:
                        print(f'Fingerprinted {len(seen):,}/{len(parents):,} parents', flush=True)
    missing = sorted(parents.keys() - seen)
    summary = dict(status='PASS' if not missing else 'INCOMPLETE',
                   required_parents=len(parents), fingerprinted_parents=len(seen),
                   missing_parents=len(missing), viral_fragments=sum(parents.values()),
                   archive=str(Path(archive_path).resolve()), member=MEMBER,
                   member_crc32=f'{info.CRC:08x}', member_compressed_size=info.compress_size,
                   truth_sha256=file_hash(truth), script_sha256=file_hash(__file__),
                   fingerprints_sha256=file_hash(table),
                   canonicalization='uppercase; U to T; minimum of forward/reverse-complement SHA256',
                   note='Sequence identity evidence only; no molecule routes assigned; no pipeline launch.')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    with (out / 'missing_parents.txt').open('w', encoding='utf-8') as handle:
        handle.writelines(key + '\n' for key in missing)
    package = out / 'routing_evidence.tar.gz'
    with tarfile.open(package, 'w:gz') as tar:
        for name in ('viral_parent_fingerprints.tsv.gz', 'summary.json', 'missing_parents.txt'):
            tar.add(out / name, arcname=name)
    print(json.dumps(summary, indent=2), flush=True)
    print(f'UPLOAD: {package.resolve()}', flush=True)
    return 0 if not missing else 2


def file_hash(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while chunk := handle.read(1024 * 1024):
            value.update(chunk)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--truth', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    return run(args.truth, args.archive, args.out)


if __name__ == '__main__':
    raise SystemExit(main())
