"""Create a deterministic private RC snapshot; never commit, deploy or launch.

Only curated runtime/source trees are packaged. Environment installations,
databases, analysis results and sibling/legacy worktrees are excluded.
"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile

TREES=('bin','envs','modules','workflows','assets','docs','tests','release')
TOP=('README.md','visum','visum.config','nextflow.config','visum_nextflow.nf','.gitignore')
EXCLUDE={'__pycache__','.pytest_cache','.git','outputs','node_modules'}
RELEASE='visum-alpha-rc3-20261020'


def digest(data):return hashlib.sha256(data).hexdigest()


def package(project,output,validation):
    project,output,validation=map(lambda p:Path(p).resolve(),(project,output,validation))
    if output.exists() or output.with_name(output.name+'.sha256').exists():raise FileExistsError(output)
    verified=json.loads(validation.read_text())
    if verified['status']!='PASS':raise ValueError('Freeze requires a passing validation report')
    data={};records={}
    files=[project/name for name in TOP]
    files += [p for folder in TREES for p in (project/folder).rglob('*')
              if p.is_file() and not EXCLUDE.intersection(p.relative_to(project).parts)]
    for path in sorted(set(files)):
        if path.is_symlink() or project not in path.resolve().parents:raise ValueError('Source escaped approved project')
        name=path.relative_to(project).as_posix()
        original=path.read_bytes()
        payload=original.replace(b'\r\n',b'\n') if path.suffix=='.sh' or name=='visum' else original
        mode=0o755 if path.suffix=='.sh' or name=='visum' or (path.suffix=='.py' and payload.startswith(b'#!')) else 0o644
        role='reference_asset' if name.startswith('assets/') else 'environment_specification' if name.startswith('envs/') else 'source'
        records[name]=dict(sha256=digest(payload),checkout_sha256=digest(original),mode=oct(mode),role=role)
        data[name]=payload
    scientific=verified['scientific_sources_unchanged']
    for name,value in scientific.items():
        if records[name]['sha256']!=value:raise ValueError('Validated scientific source changed before freeze: '+name)
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=project,text=True).strip()
    data['release/validation-summary.json']=(json.dumps(verified,indent=2,sort_keys=True)+'\n').encode()
    records['release/validation-summary.json']=dict(sha256=digest(data['release/validation-summary.json']),
        checkout_sha256=None,mode='0o644',role='validation_provenance')
    manifest=dict(schema='visum-private-alpha-source-v1',release_id=RELEASE,freeze_date='2026-10-20',git_base_commit=revision,
        source_control_state='approved isolated working-tree snapshot; not a committed/tagged public release',
        files=records,external_database_bytes_included=False,resolved_environment_builds_included=False,
        deployment_performed=False,benchmark_launched=False,student_runs_accessed=False)
    data['SOURCE_MANIFEST.json']=(json.dumps(manifest,indent=2,sort_keys=True)+'\n').encode()
    data['checksums.sha256']=''.join(digest(value)+'  '+name+'\n' for name,value in sorted(data.items())).encode()
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('wb') as stream:
        with gzip.GzipFile(filename='',fileobj=stream,mode='wb',mtime=0,compresslevel=9) as compressed:
            with tarfile.open(fileobj=compressed,mode='w',format=tarfile.PAX_FORMAT) as archive:
                for name,value in sorted(data.items()):
                    item=tarfile.TarInfo(RELEASE+'/'+name)
                    item.size=len(value);item.mtime=0;item.uid=0;item.gid=0;item.uname='';item.gname=''
                    item.mode=int(records.get(name,{}).get('mode','0o644'),8)
                    archive.addfile(item,io.BytesIO(value))
    value=digest(output.read_bytes())
    output.with_name(output.name+'.sha256').write_bytes((value+'  '+output.name+'\n').encode())
    return dict(archive=str(output),sha256=value,files=len(data),bytes=output.stat().st_size,release_id=RELEASE)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--validation-report',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(package(args.project,args.out,args.validation_report),indent=2))
