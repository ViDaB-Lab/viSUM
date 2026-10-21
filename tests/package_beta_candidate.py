"""Deterministic local beta review archive; never call Git or publish anything."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

RELEASE='viSUM-v0.9.0-beta.2'
VERSION='0.9.0-beta.2'
TOP=('README.md','LICENSE','CITATION.cff','CHANGELOG.md','visum','visum.config',
     'nextflow.config','visum_nextflow.nf','.gitignore')
TREES=('bin','envs','modules','workflows','assets','docs','tests','release')
EXCLUDE={'__pycache__','.pytest_cache','.git','outputs','node_modules'}

def digest(value):return hashlib.sha256(value).hexdigest()

def package(project,output,validation):
    project,output,validation=map(lambda p:Path(p).resolve(),(project,output,validation))
    if output.exists() or output.with_name(output.name+'.sha256').exists():raise FileExistsError(output)
    validated=json.loads(validation.read_text())
    if validated.get('status')!='PASS':raise ValueError('Passing local validation required')
    records={};data={}
    paths=[project/n for n in TOP]
    paths += [p for folder in TREES for p in (project/folder).rglob('*')
              if p.is_file() and not EXCLUDE.intersection(p.relative_to(project).parts)]
    for path in sorted(set(paths)):
        if path.is_symlink() or project not in path.resolve().parents:raise ValueError('Source escaped approved tree')
        name=path.relative_to(project).as_posix()
        original=path.read_bytes()
        payload=original.replace(b'\r\n',b'\n') if path.suffix=='.sh' or name=='visum' else original
        data[name]=payload
        mode=0o755 if path.suffix=='.sh' or name=='visum' or (path.suffix=='.py' and payload.startswith(b'#!')) else 0o644
        records[name]=dict(sha256=digest(payload),checkout_sha256=digest(original),mode=oct(mode),
            role='reference_asset' if name.startswith('assets/') else
                 'environment_specification' if name.startswith('envs/') else 'source')
    for name,value in validated['scientific_sources_unchanged'].items():
        if records[name]['sha256']!=value:raise ValueError('Validated scientific source changed: '+name)
    data['release/validation-summary.json']=(json.dumps(validated,indent=2,sort_keys=True)+'\n').encode()
    records['release/validation-summary.json']=dict(sha256=digest(data['release/validation-summary.json']),
        checkout_sha256=None,mode='0o644',role='validation_provenance')
    manifest=dict(schema='visum-beta-review-source-v1',release_id=RELEASE,version=VERSION,
        review_date='2026-10-21',git_base_commit=validated.get('public_main_commit'),
        source_control_state='local uncommitted review snapshot; no beta.2 tag/release exists',
        base_candidate=validated.get('base_candidate'),files=records,review_only=True,
        biological_payload_changed=False,external_database_bytes_included=False,
        resolved_environment_builds_included=False,publication_performed=False,
        deployment_performed=False,benchmark_launched=False,student_runs_accessed=False)
    data['SOURCE_MANIFEST.json']=(json.dumps(manifest,indent=2,sort_keys=True)+'\n').encode()
    data['checksums.sha256']=''.join(digest(v)+'  '+n+'\n' for n,v in sorted(data.items())).encode()
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('wb') as stream,gzip.GzipFile(fileobj=stream,filename='',mode='wb',mtime=0,compresslevel=9) as zipped:
        with tarfile.open(fileobj=zipped,mode='w',format=tarfile.PAX_FORMAT) as archive:
            for name,value in sorted(data.items()):
                item=tarfile.TarInfo(RELEASE+'/'+name);item.size=len(value);item.mtime=0
                item.uid=item.gid=0;item.uname=item.gname='';item.mode=int(records.get(name,{}).get('mode','0o644'),8)
                archive.addfile(item,io.BytesIO(value))
    checksum=digest(output.read_bytes())
    output.with_name(output.name+'.sha256').write_bytes((checksum+'  '+output.name+'\n').encode())
    return dict(archive=str(output),sha256=checksum,bytes=output.stat().st_size,files=len(data),release_id=RELEASE)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--validation-report',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(package(args.project,args.out,args.validation_report),indent=2))
