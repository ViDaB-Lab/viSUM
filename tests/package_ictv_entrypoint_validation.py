#!/usr/bin/env python3
"""Package only the new test harness; reuse the checksum-verified saved package."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[1]
FILES = ('tests/entrypoint_reconciliation_fixture.py', 'tests/run_ictv_entrypoint_validation.py',
         'tests/test_ictv_entrypoint_fixture.py')


def package(review, original, archive):
    from entrypoint_reconciliation_fixture import BASE_MANIFEST, interface_contracts, verify_original
    verify_original(original)
    # Main/config/scientific files in the isolated checkout must still be exactly
    # the files that passed the earlier saved replay and native module gates.
    for name in ('visum_nextflow.nf','visum.config','nextflow.config','modules/local/viharmony.nf',
                 'bin/run_viharmony.py','bin/evidence_schema.py','bin/ictv_reconciliation.py'):
        if (ROOT/name).read_bytes() != (original/name).read_bytes():
            raise ValueError('Isolated source diverged from tested source: '+name)
    with tarfile.open(review) as tar:
        status_bytes=tar.extractfile('results/validation_status.json').read()
    status=json.loads(status_bytes)
    if status['status']!='PASS' or not status['native_nextflow']['status']=='PASS':
        raise ValueError('Earlier native module gate did not pass')
    data={name:(ROOT/name).read_bytes() for name in FILES}
    manifest=dict(schema='visum-main-entrypoint-saved-input-test-v1',
        scope='Unchanged main workflow and real harmonizer; inert adapters replace all upstream modules',
        original_package_manifest_sha256=BASE_MANIFEST,
        passed_native_status_sha256=hashlib.sha256(status_bytes).hexdigest(),
        adapter_interfaces_verified=interface_contracts(ROOT),
        files={name:hashlib.sha256(value).hexdigest() for name,value in data.items()})
    data['ENTRYPOINT_PACKAGE.json']=(json.dumps(manifest,indent=2)+'\n').encode()
    if archive.exists():raise ValueError('Use a new archive filename')
    archive.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(archive,'w:gz') as tar:
        for name,value in data.items():
            item=tarfile.TarInfo('ictv-entrypoint-check/'+name)
            item.size,item.mode=len(value),0o644
            tar.addfile(item,io.BytesIO(value))
    print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size,
        sha256=hashlib.sha256(archive.read_bytes()).hexdigest()),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review',type=Path,required=True)
    parser.add_argument('--original-package',type=Path,required=True)
    parser.add_argument('--archive',type=Path,required=True)
    args=parser.parse_args()
    package(args.review,args.original_package,args.archive)
