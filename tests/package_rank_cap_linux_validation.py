#!/usr/bin/env python3
"""Freeze a small saved-input/native-cache gate, not a viSUM deployment."""
import argparse
import ast
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT / 'tests'))
from run_rank_cap_linux_validation import scientific_snapshot

SAMPLES = ('BMARK_DNA_CELLULAR','BMARK_RNA_VIRAL')
BUNDLE = 'assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz'
FILES = ('bin/run_viharmony.py','bin/evidence_schema.py','bin/ictv_reconciliation.py',
    'modules/local/viharmony.nf','visum.config', 'assets/ICTV_VMR_MSL41.csv', BUNDLE, BUNDLE+'.manifest.json',
    'tests/run_rank_cap_linux_validation.py','tests/test_rank_cap_linux_validation.py',
    'tests/test_viharmony_reporting_cap.py','tests/test_viharmony_reconciliation.py')


def package(args):
    report_bytes = (args.replay / 'report.json').read_bytes()
    report = json.loads(report_bytes)
    if report['status'] != 'PASS' or report['candidate_records'] != 3068:
        raise ValueError('Expected the passed 3,068-candidate reporting-cap replay')
    if hashlib.sha256((ROOT / 'bin/run_viharmony.py').read_bytes()).hexdigest() != report['updated_script_sha256']:
        raise ValueError('Harmonizer source differs from the completed replay')
    data = {name:(ROOT/name).read_bytes() for name in FILES}
    data['native_validation.nf'] = (ROOT/'tests/integrated_reconciliation_native.nf').read_bytes()
    data['provenance/local_saved_replay_report.json'] = report_bytes
    specs=[]
    for sample in SAMPLES:
        source=args.saved/'inputs'/(sample+'_results')
        fasta=args.saved/'fasta'/(sample+'_results')
        published=args.saved/'published'/sample
        old=json.loads((published/'updated'/(sample+'.harmonizer_manifest.json')).read_text())
        paths=dict(normalized=fasta/'prep'/(sample+'.normalized.fasta'),
            refined=fasta/'refinement'/(sample+'.refined_candidates.fasta'),
            header=source/'prep'/(sample+'.header_map.tsv'), gate=source/'discovery_gate'/(sample+'.discovery_gate.tsv'),
            regions=source/'refinement'/(sample+'.provirus_region_map.tsv'),
            groups=source/'vcontact3'/(sample+'.vcontact3_group_membership.tsv'))
        evidence=[]
        for tool in ('genomad','virsorter2','cenotetaker3','deep6','deepmicroclass2','virbot','gianthunter','vicat','checkv','tesorter','vitap','vcontact3'):
            path=published/(sample+'.virsorter2_evidence.tsv') if tool=='virsorter2' else source/tool/(sample+'.'+tool+'_evidence.tsv')
            if path.is_file(): evidence.append(path)
        evidence.append(source/'vicat'/(sample+'.vicat_refined_evidence.tsv'))
        if old['input_type']=='rna': evidence.append(source/'vicat'/(sample+'.vicat_orf_evidence.tsv'))
        files=[*paths.values(),*evidence]
        if len({p.name for p in files}) != len(files): raise ValueError('Input basename collision')
        folder='saved/'+sample
        for path in files: data[folder+'/'+path.name]=path.read_bytes()
        specs.append(dict(sample=sample,type=old['input_type'],folder=folder,bundle=BUNDLE,
            **{key:path.name for key,path in paths.items()}, evidence=[p.name for p in evidence],
            files=[folder+'/'+p.name for p in files], minimum_length=old['policy']['vcontact3_min_taxonomy_length'],
            pair_floor=old['policy']['rna_pair_homology_floor']))
    expected={rank:{sample:scientific_snapshot(args.replay/sample/rank,sample) for sample in SAMPLES}
              for rank in ('family','genus','species')}
    data['saved_specs.json']=(json.dumps(specs,indent=2)+'\n').encode()
    data['expected_outputs.json']=(json.dumps(expected,indent=2)+'\n').encode()
    data['README.txt']=b'''Reporting-cap Linux validation. NOT a viSUM installation or deployment.
Two saved benchmark panels, 1,077 candidates. Imports only the real VIHARMONY
module and a saved-file copy process. No setup, discovery or search modules.
Fresh private launch/cache/work/output directories; installed Nextflow framework
is copied read-only into a private home. No Conda or package downloads.
Stages: actual family default; unchanged resume; genus; species; restore family.
First and rank-changed harmonizers execute; unchanged/restored family tasks must
be CACHED with the original task hashes. Saved-input copies must remain cached.
Actual task outputs AND published files are verified after each completed run.
Publication must match the CURRENT task byte-for-byte, including the manifest
and every output. The verifier never copies/repairs published files. Explicit
overwrite/failOnError belongs to the real harmonizer module, not a test override.
Audit row/column order and gzip timestamps are ignored, but values, multiplicity,
all sequence bytes, decisions, review reasons and uncapped consensus are checked.
The native fixture is NOT the full main viSUM workflow and contains no searches.
Local-only mode cannot produce a native PASS. Linux native execution is the gate.
'''
    test_count=sum(sum(isinstance(n,ast.FunctionDef) and n.name.startswith('test_') for n in ast.walk(ast.parse(data['tests/'+name].decode())))
                   for name in ('test_viharmony_reporting_cap.py','test_rank_cap_linux_validation.py'))
    manifest=dict(schema='visum-reporting-cap-linux-v2-publication',samples=list(SAMPLES),candidate_records=1077,test_count=test_count,
        scope='Actual harmonizer module + saved copy process only; no searches, deployment or student runs',
        source_replay_sha256=hashlib.sha256(report_bytes).hexdigest(),
        files={name:hashlib.sha256(value).hexdigest() for name,value in data.items()})
    data['RANK_CAP_PACKAGE.json']=(json.dumps(manifest,indent=2)+'\n').encode()
    if args.archive.exists(): raise ValueError('Choose a new archive path; never overwrite a package')
    args.archive.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(args.archive,'w:gz') as tar:
        for name,value in data.items():
            item=tarfile.TarInfo('viSUM-rank-cap-check/'+name)
            item.size,item.mode,item.mtime=len(value),0o644,0
            tar.addfile(item,io.BytesIO(value))
    print(json.dumps(dict(archive=str(args.archive),bytes=args.archive.stat().st_size,
        sha256=hashlib.sha256(args.archive.read_bytes()).hexdigest(),files=len(data),test_count=test_count),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--saved',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True)
    parser.add_argument('--archive',type=Path,required=True)
    package(parser.parse_args())
