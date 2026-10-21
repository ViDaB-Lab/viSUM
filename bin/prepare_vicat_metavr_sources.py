#!/usr/bin/env python3
"""Acquire verified MetaVR portal sources; adapt format, never biological policy.

The packaged portal manifest binds verified native protein and patched-metadata
evidence to exact source/patch checksums. Pending manifests fail BEFORE network
access. `compare` is read-only on inputs and never builds a database.
The setup caller holds the existing destination installation lock.
"""
from __future__ import annotations

import argparse
import csv
import errno
import gzip
import hashlib
from http.client import IncompleteRead, HTTPException
import json
from pathlib import Path
import re
import ssl
import tempfile
import time
from urllib.parse import urlparse
from urllib.request import Request, urlopen
from urllib.error import HTTPError

from prepare_vicat_metadata import REQUIRED_COLUMNS, OPTIONAL_COLUMNS

HASH = re.compile(r'[a-f0-9]{64}\Z')
NAMES = {'proteins': 'METAVR.faa.bgz', 'metadata': 'METAVR_main_table.parquet'}


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + '.part')
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    partial.replace(path)


def load_manifest(path, frozen=None, verified=True):
    manifest = json.loads(Path(path).read_text(encoding='utf-8'))
    if manifest.get('schema') != 'visum-metavr-sources-v1':
        raise ValueError('Unsupported MetaVR source manifest schema')
    if frozen is not None and manifest.get('frozen') != frozen:
        raise ValueError('MetaVR manifest does not match the configured frozen source identities/counts')
    if verified and manifest.get('compatibility') != 'verified-full-source-equivalence':
        raise ValueError('MetaVR portal compatibility is NOT yet verified against the frozen sources. '
                         'No download or database build was started. Supply the original local source '
                         'pair/prepared database, or complete the read-only full-source comparison first.')
    for key, name in NAMES.items():
        entry = manifest['files'][key]
        parsed = urlparse(entry['url'])
        if entry['name'] != name or parsed.scheme != 'https' or not parsed.netloc or parsed.username:
            raise ValueError('Invalid MetaVR source name/HTTPS URL: ' + key)
        if not isinstance(entry['bytes'], int) or entry['bytes'] <= 0:
            raise ValueError('Invalid MetaVR source byte count: ' + key)
        if verified and not HASH.fullmatch(entry.get('sha256') or ''):
            raise ValueError('Verified MetaVR sources require full SHA256 checksums: ' + key)
    patch_bundle(path, manifest)
    if verified and manifest.get('validation_evidence'):
        evidence = manifest['validation_evidence']
        base = Path(path).resolve().parent
        evidence_path = (base / evidence['path']).resolve()
        if base not in evidence_path.parents or digest(evidence_path) != evidence['sha256']:
            raise ValueError('MetaVR source validation evidence missing/changed')
    return manifest


def patch_bundle(manifest_path, manifest):
    """Resolve and verify an optional packaged compatibility recipe; no network."""
    descriptor = manifest.get('metadata_compatibility_patch')
    if descriptor is None:
        return None
    base = Path(manifest_path).resolve().parent
    paths = {}
    for key in ('manifest', 'patch'):
        relative = Path(descriptor[key])
        path = (base / relative).resolve()
        if relative.is_absolute() or base not in path.parents or not path.is_file():
            raise ValueError('MetaVR patch asset missing or escaped manifest directory: ' + key)
        paths[key] = path
    if not HASH.fullmatch(descriptor.get('manifest_sha256') or '') or digest(paths['manifest']) != descriptor['manifest_sha256']:
        raise ValueError('MetaVR compatibility recipe checksum mismatch')
    recipe = json.loads(paths['manifest'].read_text(encoding='utf-8'))
    if (recipe.get('schema') != 'visum-metavr-votu-compatibility-v1' or
            recipe.get('original_metadata_sha256') != manifest['frozen']['metadata_sha256'] or
            recipe.get('source_identity', {}).get('original_metadata_sha256') != manifest['frozen']['metadata_sha256'] or
            recipe.get('source_identity', {}).get('portal_sha256') != manifest['files']['metadata']['sha256'] or
            recipe.get('expected_rows') != manifest['frozen']['uvigs']):
        raise ValueError('MetaVR compatibility recipe does not match pinned source identities/count')
    if not HASH.fullmatch(recipe.get('patch_sha256') or '') or digest(paths['patch']) != recipe['patch_sha256']:
        raise ValueError('MetaVR compatibility patch checksum mismatch')
    return dict(paths=paths, identity=dict(recipe_sha256=descriptor['manifest_sha256'],
        patch_sha256=recipe['patch_sha256'], expected_changes=recipe['expected_changes'],
        purpose=recipe['purpose']))


def prepare_metadata(manifest_path, manifest, portal, destination, work):
    """Shared metadata path for compare/setup; never accesses protein sources."""
    portal, destination, work = map(Path, (portal, destination, work))
    if destination.exists():
        raise FileExistsError('Refusing to overwrite prepared metadata')
    if destination.resolve() == portal.resolve():
        raise ValueError('Prepared metadata must be a new copy')
    bundle = patch_bundle(manifest_path, manifest)
    if bundle is None:
        convert_metadata(portal, destination, manifest['frozen']['uvigs'], work)
        return None
    entry = manifest['files']['metadata']
    if portal.stat().st_size != entry['bytes'] or digest(portal) != entry['sha256']:
        raise ValueError('Portal metadata does not match compatibility recipe source; no patch applied')
    from patch_vicat_metavr_metadata import apply_patch_copy, verify_non_votu_columns
    work.mkdir(parents=True, exist_ok=True)
    raw = work / 'METAVR.portal_metadata.tsv.gz'
    if raw.resolve() == destination.resolve():
        raise ValueError('Prepared metadata must not overwrite the raw conversion')
    if not raw.exists():
        convert_metadata(portal, raw, manifest['frozen']['uvigs'], work)
    # The validated recipe also pins the exact converted input. No automatic
    # rebinding of a patch to different bytes/schema/order is permitted.
    application = apply_patch_copy(raw, bundle['paths']['patch'], bundle['paths']['manifest'],
                                   destination, work, connect)
    untouched = verify_non_votu_columns(raw, destination, work, manifest['frozen']['uvigs'], connect)
    return dict(bundle['identity'], raw_metadata_sha256=digest(raw),
                patched_metadata_sha256=application['patched_sha256'], applied_rows=application['applied_rows'],
                non_votu_changed_rows=untouched['changed_rows'])


def prepare_metadata_cached(manifest_path, manifest, portal, destination, work, previous=None):
    """Shared metadata cache check; no protein I/O and no full-source approval."""
    identity = digest(manifest_path)
    bundle = patch_bundle(manifest_path, manifest)
    if previous is not None and previous.get('manifest_sha256') != identity:
        raise ValueError('MetaVR source manifest changed; use a new source directory')
    destination = Path(destination)
    if destination.exists():
        if previous is None or digest(destination) != previous.get('metadata_sha256'):
            raise ValueError('Completed converted MetaVR metadata changed/unreceipted; use a new source directory')
        metadata_patch = previous.get('metadata_compatibility_patch')
        if bundle is not None and (not isinstance(metadata_patch, dict) or any(
                metadata_patch.get(key) != value for key, value in bundle['identity'].items())):
            raise ValueError('Cached MetaVR patch provenance missing or changed')
        status = 'reused'
    else:
        metadata_patch = prepare_metadata(manifest_path, manifest, portal, destination, work)
        status = 'prepared'
    return dict(cache_status=status, metadata_sha256=digest(destination),
                metadata_compatibility_patch=metadata_patch)


def retryable_transfer_error(error):
    if isinstance(error, HTTPError):
        return error.code in (408, 425, 429, 500, 502, 503, 504)
    if isinstance(error, ssl.SSLCertVerificationError) or isinstance(getattr(error,'reason',None), ssl.SSLCertVerificationError):
        return False
    if isinstance(error, OSError) and error.errno in (errno.ENOSPC, errno.EDQUOT, errno.EACCES, errno.EPERM, errno.EIO):
        return False
    return True


def fetch(entry, destination, retries=10, timeout=120, chunk_bytes=256*1024**2):
    """Bounded Range transfers; preserve every completed prefix on disconnect.

    Retry budget counts CONSECUTIVE failed requests, resetting after a complete
    range. Identity/range/size violations are fatal, not transient. A server
    ignoring a resume range NEVER causes a partial file to be truncated.
    """
    if retries < 1 or chunk_bytes < 1:
        raise ValueError('Retries and chunk bytes must be positive')
    final = Path(destination)
    directory = final.parent
    expected_bytes = entry['bytes']
    url = entry['url']
    if final.exists():
        if final.stat().st_size != expected_bytes or digest(final) != entry['sha256']:
            raise ValueError('Completed MetaVR source changed; use a new source directory: ' + str(final))
        return
    directory.mkdir(parents=True, exist_ok=True)
    partial = final.with_name(final.name + '.part')
    identity_path = final.with_name(final.name + '.http_identity.json')
    if partial.exists() and partial.stat().st_size > 0 and not identity_path.is_file():
        raise ValueError('Existing MetaVR source partial lacks HTTP identity; refusing to append')
    with urlopen(Request(url, method='HEAD', headers={'Accept-Encoding':'identity'}), timeout=timeout) as response:
        if response.getcode() != 200 or int(response.headers.get('Content-Length', -1)) != expected_bytes:
            raise ValueError('Portal MetaVR source byte count changed; no transfer started')
        identity = dict(url=url, bytes=expected_bytes, etag=response.headers.get('ETag'),
                        last_modified=response.headers.get('Last-Modified'))
    if not identity['etag'] and not identity['last_modified']:
        raise ValueError('Portal lacks a stable HTTP validator for safe MetaVR source resume')
    if identity_path.exists():
        if json.loads(identity_path.read_text()) != identity:
            raise ValueError('Portal MetaVR source HTTP identity changed; existing partial was NOT modified')
    else:
        atomic_json(identity_path, identity)
    failures = 0
    while True:
        offset = partial.stat().st_size if partial.exists() else 0
        if offset > expected_bytes:
            raise ValueError('Oversized partial MetaVR source; refusing to use it')
        if offset == expected_bytes:
            break
        end = min(offset + chunk_bytes, expected_bytes) - 1
        headers = {'Accept-Encoding':'identity', 'User-Agent':'viSUM-MetaVR-setup/2',
                   'Range':f'bytes={offset}-{end}',
                   'If-Range':identity['etag'] or identity['last_modified']}
        try:
            with urlopen(Request(url, headers=headers), timeout=timeout) as response:
                status = response.getcode()
                if status == 206:
                    if response.headers.get('Content-Range') != f'bytes {offset}-{end}/{expected_bytes}':
                        raise ValueError('Unexpected MetaVR source HTTP range; existing partial was NOT truncated')
                elif status == 200 and offset == 0 and end == expected_bytes - 1:
                    pass  # Complete small file, e.g. offline regression fixtures.
                else:
                    raise ValueError('Server ignored/refused the bounded Range request; existing partial was NOT truncated')
                for name, key in (('ETag','etag'),('Last-Modified','last_modified')):
                    if response.headers.get(name) and identity[key] and response.headers[name] != identity[key]:
                        raise ValueError('MetaVR source HTTP identity changed during transfer; partial retained')
                length = response.headers.get('Content-Length')
                if length is not None and int(length) != end - offset + 1:
                    raise ValueError('MetaVR source response size differs from requested range')
                if (partial.stat().st_size if partial.exists() else 0) != offset:
                    raise ValueError('Partial file changed concurrently; stop duplicate download jobs')
                print(f'Downloading portal MetaVR source bytes {offset:,}-{end:,} / {expected_bytes:,}...', flush=True)
                last_progress = time.monotonic()
                with partial.open('ab') as output:
                    while True:
                        block = response.read(min(1024*1024, end + 1 - output.tell() + 1))
                        if not block:
                            break
                        if output.tell() + len(block) > end + 1:
                            raise ValueError('MetaVR source-only byte ceiling exceeded')
                        output.write(block)
                        if time.monotonic() - last_progress >= 30:
                            print(f'Portal download: {output.tell():,} / {expected_bytes:,} bytes ({100*output.tell()/expected_bytes:.1f}%)', flush=True)
                            last_progress = time.monotonic()
            if partial.stat().st_size != end + 1:
                raise OSError('Incomplete MetaVR source range; saved prefix retained')
        except (OSError, HTTPException) as error:
            if not retryable_transfer_error(error):
                raise
            failures += 1
            saved = partial.stat().st_size if partial.exists() else 0
            print(f'Transfer interrupted ({type(error).__name__}: {error}); retained {saved:,} bytes. Consecutive failure {failures}/{retries}.', flush=True)
            if failures >= retries:
                raise
            delay = min(2 ** (failures - 1), 30)
            print(f'Retrying from saved byte offset in {delay} seconds...', flush=True)
            time.sleep(delay)
            continue
        failures = 0
    if digest(partial) != entry['sha256']:
        raise ValueError('MetaVR SHA256 mismatch; no source was published: ' + str(final))
    partial.replace(final)

def connect(work):
    import duckdb
    connection = duckdb.connect()
    connection.execute("SET threads = 1")
    connection.execute("SET memory_limit = '1GB'")
    connection.execute('SET temp_directory = ' + sql(str(Path(work) / 'duckdb-tmp')))
    return connection


def sql(value):
    return "'" + str(value).replace("'", "''") + "'"


def convert_metadata(source, destination, expected_rows, work):
    """Export only existing compact-metadata columns with nulls as explicit \\N.

    The unchanged compact-metadata helper and taxonomy builder consume this TSV.
    No filtering, taxonomy translation, vOTU substitution or deduplication here.
    """
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.part')
    with connect(work) as connection:
        columns = [row[0] for row in connection.execute('DESCRIBE SELECT * FROM read_parquet(?)', [str(source)]).fetchall()]
        missing = set(REQUIRED_COLUMNS) - set(columns)
        if missing:
            raise ValueError('MetaVR Parquet is missing required columns: ' + ', '.join(sorted(missing)))
        rows, unique, empty = connection.execute('SELECT count(*), count(DISTINCT uvig), '
            "count(*) FILTER (WHERE uvig IS NULL OR trim(uvig) = '') FROM read_parquet(?)", [str(source)]).fetchone()
        if rows != expected_rows or rows != unique or empty:
            raise ValueError(f'MetaVR metadata row/unique-ID mismatch: {rows} rows, {unique} IDs, {empty} empty; expected {expected_rows}')
        retained = REQUIRED_COLUMNS + [name for name in OPTIONAL_COLUMNS if name in columns]
        selection = ', '.join('"' + name + '"' for name in retained)
        connection.execute(f'COPY (SELECT {selection} FROM read_parquet({sql(source)})) '
            f'TO {sql(partial)} (FORMAT CSV, DELIMITER \'\t\', HEADER true, COMPRESSION gzip, NULL \'\\N\')')
    # Complete the gzip read to catch truncation/CRC errors before publication.
    with gzip.open(partial, 'rb') as source_handle:
        for _ in iter(lambda: source_handle.read(1024 * 1024), b''):
            pass
    partial.replace(destination)
    return rows


def protein_fingerprint(path):
    """Order-sensitive exact ID/sequence digest; ignore only description/wrapping.

    Preserving input order is conservative: order can affect representative choice.
    Sequence case/content and first-token protein identifiers are not rewritten.
    """
    result = hashlib.sha256()
    count = 0
    has_sequence = False
    with gzip.open(path, 'rt', encoding='ascii') as source:
        for raw in source:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if count and not has_sequence:
                    raise ValueError('Empty protein record')
                identifier = line[1:].split(maxsplit=1)
                if not identifier:
                    raise ValueError('Empty protein identifier')
                result.update(b'\n>' + identifier[0].encode('ascii') + b'\n')
                count += 1
                has_sequence = False
            else:
                if not count or not re.fullmatch(r'[A-Za-z*.-]+', line):
                    raise ValueError('Invalid protein FASTA')
                result.update(line.encode('ascii'))
                has_sequence = True
    if not count or not has_sequence:
        raise ValueError('Empty/incomplete protein FASTA')
    return {'records': count, 'ordered_id_sequence_sha256': result.hexdigest()}


def taxonomy_fingerprint(path, work):
    # Reuse the ACTUAL builder's normalization instead of creating new rules.
    from build_vicat_taxonomy_lookup import build_metadata_sql
    result = hashlib.sha256()
    count = 0
    with connect(work) as connection:
        reader = connection.execute('SELECT * FROM (' + build_metadata_sql(Path(path)) + ') ORDER BY uvig')
        while rows := reader.fetchmany(10000):
            for row in rows:
                result.update((json.dumps(row, ensure_ascii=True, separators=(',', ':')) + '\n').encode('ascii'))
                count += 1
    return {'records': count, 'builder_normalized_taxonomy_sha256': result.hexdigest()}


def compare_sources(template, frozen_proteins, frozen_metadata, portal_proteins, portal_metadata, output, work):
    """Compare complete local inputs. NO download, build, changes to source files.

    Emit a usable manifest only if the frozen compressed-file pins and all exact
    protein IDs/sequences/order and builder-normalized metadata agree.
    """
    output, work = Path(output), Path(work)
    if output.exists():
        raise FileExistsError(output)
    manifest = load_manifest(template, verified=False)
    frozen = manifest['frozen']
    for key, path in (('protein_sha256', frozen_proteins), ('metadata_sha256', frozen_metadata)):
        if digest(path) != frozen[key]:
            raise ValueError('Frozen original source checksum mismatch: ' + key)
    work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='metavr-compare-', dir=work) as temporary:
        converted = Path(temporary) / 'portal.tsv.gz'
        metadata_patch = prepare_metadata(template, manifest, portal_metadata, converted, temporary)
        old_protein, new_protein = protein_fingerprint(frozen_proteins), protein_fingerprint(portal_proteins)
        old_metadata, new_metadata = taxonomy_fingerprint(frozen_metadata, temporary), taxonomy_fingerprint(converted, temporary)
        if old_protein != new_protein or old_protein['records'] != frozen['proteins']:
            raise ValueError('Portal protein IDs/sequences/order differ from the frozen source; no manifest approved')
        if old_metadata != new_metadata or old_metadata['records'] != frozen['uvigs']:
            raise ValueError('Portal taxonomy/vOTU/method metadata differ from the frozen source; no manifest approved')
    for key, path in (('proteins', portal_proteins), ('metadata', portal_metadata)):
        manifest['files'][key].update(bytes=Path(path).stat().st_size, sha256=digest(path))
    manifest['compatibility'] = 'verified-full-source-equivalence'
    manifest['comparison'] = {'proteins': old_protein, 'metadata': old_metadata,
        'scope': 'all ordered first-token protein IDs/exact sequences; all builder-normalized metadata rows',
        'database_build_performed': False}
    if metadata_patch is not None:
        manifest['comparison']['metadata_compatibility_patch'] = metadata_patch
        # Asset paths are relative to the manifest. A moved approved manifest
        # must be accompanied by its assets; refuse to emit an unusable bundle.
        patch_bundle(output, manifest)
    manifest.pop('note', None)
    atomic_json(output, manifest)
    return manifest


def prepare_sources(manifest_path, root, frozen):
    # Must precede directory creation, GET, conversion, and any builder invocation.
    manifest = load_manifest(manifest_path, frozen)
    root = Path(root)
    identity = digest(manifest_path)
    receipt = root / 'download_provenance.json'
    previous = None
    if receipt.exists():
        previous = json.loads(receipt.read_text(encoding='utf-8'))
        if previous.get('manifest_sha256') != identity:
            raise ValueError('MetaVR source manifest changed; use a new source directory')
    for key in NAMES:
        fetch(manifest['files'][key], root / NAMES[key])
    converted = root / 'METAVR.vicat_metadata.tsv.gz'
    cached = prepare_metadata_cached(manifest_path, manifest, root / NAMES['metadata'], converted, root, previous)
    metadata_patch = cached['metadata_compatibility_patch']
    if cached['cache_status'] == 'prepared':
        if metadata_patch is not None:
            approved = manifest.get('comparison', {}).get('metadata', {})
            if 'prepared_metadata_sha256' in approved:
                valid = (approved.get('records') == frozen['uvigs'] and
                         cached['metadata_sha256'] == approved['prepared_metadata_sha256'])
            else:
                valid = taxonomy_fingerprint(converted, root) == approved
            if not valid:
                raise ValueError('Patched metadata differs from approved full-source comparison')
    provenance = {'schema': 'visum-metavr-download-receipt-v1', 'manifest_sha256': identity,
        'compatibility': manifest['compatibility'], 'frozen': frozen, 'files': manifest['files'],
        'proteins': str((root / NAMES['proteins']).resolve()), 'metadata': str(converted.resolve()),
        'protein_sha256': manifest['files']['proteins']['sha256'], 'metadata_sha256': cached['metadata_sha256']}
    if metadata_patch is not None:
        provenance['metadata_compatibility_patch'] = metadata_patch
    if previous != provenance:
        atomic_json(receipt, provenance)
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    check = sub.add_parser('check-manifest', help='No download or build; reject an unverified manifest')
    acquire = sub.add_parser('prepare', help='Download verified sources and convert metadata; NEVER build')
    for command in (check, acquire):
        command.add_argument('--manifest', required=True, type=Path)
        command.add_argument('--frozen-protein-sha256', required=True)
        command.add_argument('--frozen-metadata-sha256', required=True)
        command.add_argument('--expected-uvigs', type=int, required=True)
        command.add_argument('--expected-proteins', type=int, required=True)
    acquire.add_argument('--source-dir', type=Path, required=True)
    compare = sub.add_parser('compare', help='Read-only complete-source equivalence check; no network or build')
    compare.add_argument('--template', type=Path, required=True)
    for name in ('frozen-proteins', 'frozen-metadata', 'portal-proteins', 'portal-metadata', 'output-manifest', 'work-dir'):
        compare.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'compare':
        result = compare_sources(args.template, args.frozen_proteins, args.frozen_metadata,
            args.portal_proteins, args.portal_metadata, args.output_manifest, args.work_dir)
    else:
        frozen = {'protein_sha256': args.frozen_protein_sha256, 'metadata_sha256': args.frozen_metadata_sha256,
            'uvigs': args.expected_uvigs, 'proteins': args.expected_proteins}
        result = load_manifest(args.manifest, frozen) if args.command == 'check-manifest' else prepare_sources(args.manifest, args.source_dir, frozen)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
