"""Exercise actual approved source preparation and reuse, never a database build."""
import argparse
from contextlib import ExitStack
import json
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import traceback
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'bin'))
import prepare_vicat_metavr_sources as sources


def identity(path):
    stat=Path(path).stat()
    return dict(bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,ctime_ns=stat.st_ctime_ns,
                inode=stat.st_ino,device=stat.st_dev)


def validate(proteins, metadata, output, project=ROOT, minimum_free_bytes=64*1024**3):
    proteins,metadata,output,project=map(lambda p:Path(p).resolve(),(proteins,metadata,output,project))
    private=project.parent
    if private not in output.parents or project == output or project in output.parents:
        raise ValueError('Output must be inside the private test directory, outside package')
    if output.exists(): raise FileExistsError('Use a new private results directory')
    for source in (proteins,metadata):
        if not source.is_file(): raise FileNotFoundError(source)
        if output == source or output in source.parents:
            raise ValueError('Output overlaps existing source files')
    output.mkdir(parents=True)
    report=dict(status='FAIL',network_allowed=False,protein_sequence_rescan_performed=False,
        database_build_performed=False,student_analyses_accessed=False,deployment_performed=False,
        resource_policy=dict(duckdb_threads=1,memory_limit='1GB',spill_limit='32GB'))
    before={str(path):identity(path) for path in (proteins,metadata)}
    boundaries={}
    try:
        import duckdb
        if duckdb.__version__!='1.5.4': raise ValueError('Existing DuckDB 1.5.4 environment required')
        manifest_path=project/'assets/vicat/metavr_portal_v1.json'
        manifest=sources.load_manifest(manifest_path)
        report.update(manifest_sha256=sources.digest(manifest_path),duckdb_version=duckdb.__version__,
            python_version=sys.version,compatibility=manifest['compatibility'])
        if shutil.disk_usage(output).free < minimum_free_bytes:
            raise ValueError('At least 64 GiB free private scratch disk required; no preparation started')
        cache=output/'source_cache'
        cache.mkdir()
        # Source aliases only: never copy or modify the 50-GB input.
        for key,path in (('proteins',proteins),('metadata',metadata)):
            entry=manifest['files'][key]
            if path.stat().st_size!=entry['bytes']: raise ValueError('Existing source size differs: '+key)
            (cache/sources.NAMES[key]).symlink_to(path)
        with ExitStack() as stack:
            for name in ('urlopen','protein_fingerprint','taxonomy_fingerprint'):
                boundaries[name]=stack.enter_context(patch.object(sources,name,
                    side_effect=AssertionError('Forbidden network/repeated full comparison: '+name)))
            for owner,name in ((socket,'create_connection'),(socket,'getaddrinfo'),(subprocess,'Popen')):
                boundaries[name]=stack.enter_context(patch.object(owner,name,
                    side_effect=AssertionError('No network or external build/search programs')))
            real_connect=sources.connect
            def bounded(work):
                connection=real_connect(work)
                connection.execute("SET max_temp_directory_size='32GB'")
                return connection
            stack.enter_context(patch.object(sources,'connect',side_effect=bounded))
            print('First actual source preparation: verify existing compressed checksums; prepare private metadata.',flush=True)
            first=sources.prepare_sources(manifest_path,cache,manifest['frozen'])
            expected=manifest['comparison']['metadata']
            if first['metadata_sha256']!=expected['prepared_metadata_sha256']:
                raise ValueError('Prepared metadata approval mismatch')
            applied=first['metadata_compatibility_patch']
            if (applied['applied_rows']!=manifest['comparison']['metadata_compatibility_patch']['applied_rows'] or
                    applied['non_votu_changed_rows']!=0):
                raise ValueError('Unexpected patch effects')
            paths=(cache/'download_provenance.json',cache/'METAVR.portal_metadata.tsv.gz',
                   cache/'METAVR.vicat_metadata.tsv.gz')
            snapshot={p.name:dict(identity=identity(p),sha256=sources.digest(p)) for p in paths}
            print('Second actual preparation: checksum verification and cached reuse; no conversion or patching allowed.',flush=True)
            with (patch.object(sources,'prepare_metadata',side_effect=AssertionError('No repatching')),
                  patch.object(sources,'convert_metadata',side_effect=AssertionError('No reconversion'))):
                second=sources.prepare_sources(manifest_path,cache,manifest['frozen'])
            after={p.name:dict(identity=identity(p),sha256=sources.digest(p)) for p in paths}
            if first!=second or snapshot!=after:
                raise ValueError('Unchanged preparation altered provenance/cache contents or timestamps')
            report.update(first_preparation=first,reuse=second,cached_hashes_and_timestamps_unchanged=True,
                existing_compressed_sources_checksum_verified_twice=True)
        report['forbidden_boundary_calls']={name:value.call_count for name,value in boundaries.items()}
        if any(report['forbidden_boundary_calls'].values()): raise ValueError('Forbidden boundary attempted')
        report['status']='PASS'
    except Exception as error:
        report.update(status='FAIL',error=str(error),traceback=traceback.format_exc())
        traceback.print_exc()
    finally:
        report['inputs_unchanged']={path:identity(Path(path))==value for path,value in before.items()}
        if not all(report['inputs_unchanged'].values()):
            report.update(status='FAIL',error='Existing source identity changed during test')
        sources.atomic_json(output/'report.json',report)
    print('Integrated source preparation/reuse: '+report['status'],flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proteins',type=Path,required=True)
    parser.add_argument('--metadata',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    result=validate(args.proteins,args.metadata,args.output_dir)
    raise SystemExit(0 if result['status']=='PASS' else 1)
