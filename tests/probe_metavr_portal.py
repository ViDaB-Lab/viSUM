#!/usr/bin/env python3
"""Explicit, bounded live-format probe. NEVER download complete portal files.

Read HTTP headers, <=64 KiB of proteins, and <=1 MiB of Parquet footer metadata.
A response ignoring Range is rejected before reading its body. Schema/counts
and a FASTA sample are NOT a full frozen-source equivalence test.
"""
import argparse
import json
from pathlib import Path
import struct
import sys
import tempfile
from urllib.request import Request, urlopen
import zlib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from prepare_vicat_metavr_sources import atomic_json, load_manifest, REQUIRED_COLUMNS


def bounded_range(url, start, end, total):
    request = Request(url, headers={'Range': f'bytes={start}-{end}', 'Accept-Encoding': 'identity'})
    with urlopen(request, timeout=30) as response:
        if response.getcode() != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{total}':
            raise ValueError('Server ignored/misreported Range; refused to read a potentially complete database')
        payload = response.read(end - start + 2)
    if len(payload) != end - start + 1:
        raise ValueError('Unexpected bounded response length')
    return payload


def capture(manifest_path, directory):
    manifest = load_manifest(manifest_path, verified=False)
    headers = {}
    sizes = {}
    for name, entry in manifest['files'].items():
        with urlopen(Request(entry['url'], method='HEAD'), timeout=30) as response:
            sizes[name] = int(response.headers['Content-Length'])
            headers[name] = {'bytes': sizes[name], 'last_modified': response.headers.get('Last-Modified'),
                'accept_ranges': response.headers.get('Accept-Ranges')}
        if sizes[name] != entry['bytes']:
            raise ValueError('Portal file length differs from the recorded candidate: ' + name)
    protein = bounded_range(manifest['files']['proteins']['url'], 0, 65535, sizes['proteins'])
    url = manifest['files']['metadata']['url']
    total = sizes['metadata']
    tail = bounded_range(url, total - 8, total - 1, total)
    footer_length = struct.unpack('<I', tail[:4])[0]
    if tail[4:] != b'PAR1' or not 0 < footer_length <= 1024 * 1024 or footer_length > total - 12:
        raise ValueError('Parquet footer is invalid or exceeds the explicit 1-MiB probe budget')
    footer = bounded_range(url, total - 8 - footer_length, total - 9, total)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / 'protein-prefix.bgz.part').write_bytes(protein)
    (directory / 'parquet-footer-only.bin').write_bytes(b'PAR1' + footer + tail)
    result = {'file_headers': headers, 'binary_response_bytes': len(protein) + len(tail) + len(footer),
              'expected_uvigs': manifest['frozen']['uvigs']}
    atomic_json(directory / 'capture.json', result)
    return result


def inspect_capture(directory, output):
    directory = Path(directory)
    captured = json.loads((directory / 'capture.json').read_text())
    protein = (directory / 'protein-prefix.bgz.part').read_bytes()
    decoded = zlib.decompressobj(31).decompress(protein, 1024 * 1024).decode('ascii')
    protein_ids = [line[1:].split()[0] for line in decoded.splitlines() if line.startswith('>')]
    if not protein_ids or any('|' not in value for value in protein_ids):
        raise ValueError('Sampled protein IDs do not support the existing first-| UViG join')
    # A metadata-only local envelope: data offsets are deliberately unavailable.
    # Only parquet_schema/file_metadata are allowed; NO sequence rows are queried.
    import duckdb
    path = directory / 'parquet-footer-only.bin'
    with duckdb.connect() as connection:
        schema = connection.execute('SELECT name, type FROM parquet_schema(?)', [str(path)]).fetchall()
        rows = connection.execute('SELECT num_rows FROM parquet_file_metadata(?)', [str(path)]).fetchone()[0]
    columns = {name for name, _ in schema}
    if not set(REQUIRED_COLUMNS) <= columns or rows != captured['expected_uvigs']:
        raise ValueError('Portal metadata schema/row count is incompatible with the frozen builder expectation')
    result = {'status': 'PASS', 'scope': 'bounded schema/count/identifier-format probe ONLY',
        'full_source_equivalence': 'NOT VERIFIED', 'full_database_downloaded': False,
        'database_build_performed': False, 'student_runs_accessed': False,
        'file_headers': captured['file_headers'], 'protein_sample_ids': protein_ids[:10],
        'metadata_rows': rows, 'metadata_schema': schema,
        'binary_response_bytes': captured['binary_response_bytes'],
        'binary_response_budget': 65536 + 8 + 1024 * 1024}
    atomic_json(output, result)
    return result


def probe(manifest_path, output):
    with tempfile.TemporaryDirectory(prefix='metavr-probe-') as temporary:
        capture(manifest_path, temporary)
        return inspect_capture(temporary, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=ROOT / 'assets/vicat/metavr_portal_v1.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--capture-only', type=Path, help='Save only bounded byte ranges; do not load DuckDB')
    parser.add_argument('--inspect-capture', type=Path, help='Inspect a prior bounded capture entirely offline')
    args = parser.parse_args()
    if args.capture_only and args.inspect_capture:
        parser.error('Choose capture-only OR inspect-capture')
    if args.capture_only:
        result = capture(args.manifest, args.capture_only)
        atomic_json(args.output, result)
    elif args.inspect_capture:
        result = inspect_capture(args.inspect_capture, args.output)
    else:
        result = probe(args.manifest, args.output)
    print(json.dumps(result, indent=2))
