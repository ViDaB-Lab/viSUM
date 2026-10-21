#!/usr/bin/env python3
"""Independently reconcile prepared FASTAs, unique truth and paired execution map.

Does not launch tools or calculate pooled/union performance metrics.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
from prepare_genomad_benchmark import execution_routes, fasta_records, sha256


def verify(directory):
    root=Path(directory).resolve()
    report=json.loads((root/'preparation_summary.json').read_text())
    if report['status']!='INPUTS_PREPARED_NOT_LAUNCHED':
        raise ValueError('Preparation is incomplete or failed')
    checks={}
    for line in (root/'checksums.sha256').read_text().splitlines():
        digest,name=line.split('  ',1)
        path=(root/name).resolve()
        if path.parent!=root or name in checks:
            raise ValueError('Invalid checksum path or duplicate')
        if sha256(path)!=digest:
            raise ValueError(f'Checksum mismatch: {name}')
        checks[name]=digest
    for name in ('preparation_summary.json','truth.routed.tsv','execution_map.tsv','manifest.csv'):
        if name not in checks: raise ValueError(f'Missing checksum: {name}')
    counts=Counter()
    with tempfile.TemporaryDirectory(prefix='visum-verify-') as temp:
        db=sqlite3.connect(Path(temp)/'verify.sqlite')
        try:
            db.executescript('''CREATE TABLE expected(id TEXT,route TEXT,stratum TEXT,hash TEXT,length INTEGER,cls TEXT,parent TEXT,seen INTEGER DEFAULT 0,PRIMARY KEY(id,route));
                CREATE TABLE truth_ids(id TEXT PRIMARY KEY);
                CREATE TABLE executions(batch TEXT,id TEXT,route TEXT,hash TEXT,length INTEGER,seen INTEGER DEFAULT 0,PRIMARY KEY(batch,id),UNIQUE(id,route));''')
            with (root/'truth.routed.tsv').open(encoding='utf-8',newline='') as handle:
                for row in csv.DictReader(handle,delimiter='\t'):
                    if row['source_class'] not in ('virus','chromosome','plasmid'):
                        raise ValueError('Invalid truth class')
                    if row['source_class']!='virus' and row['input_type']!='dna':
                        raise ValueError('Negative genomic controls must use DNA only')
                    db.execute('INSERT INTO truth_ids VALUES (?)',(row['sequence_id'],))
                    counts['unique_inputs']+=1
                    for route,stratum in execution_routes(row['input_type'],report['unresolved_policy']):
                        db.execute('INSERT INTO expected(id,route,stratum,hash,length,cls,parent) VALUES (?,?,?,?,?,?,?)',
                                   (row['sequence_id'],route,stratum,row['sequence_sha256'],int(row['length_nt']),row['source_class'],row['parent_sequence_id']))
                    if counts['unique_inputs']%10000==0: db.commit()
            with (root/'execution_map.tsv').open(encoding='utf-8',newline='') as handle:
                for row in csv.DictReader(handle,delimiter='\t'):
                    expected=db.execute('SELECT stratum,hash,length,cls,parent,seen FROM expected WHERE id=? AND route=?',
                                        (row['sequence_id'],row['input_type'])).fetchone()
                    actual=(row['evaluation_stratum'],row['sequence_sha256'],int(row['length_nt']),row['source_class'],row['parent_sequence_id'])
                    if expected is None or expected[:5]!=actual or expected[5]:
                        raise ValueError('Unexpected, duplicate or mismatched execution-map row')
                    db.execute('UPDATE expected SET seen=1 WHERE id=? AND route=?',(row['sequence_id'],row['input_type']))
                    db.execute('INSERT INTO executions(batch,id,route,hash,length) VALUES (?,?,?,?,?)',
                               (row['benchmark_batch'],row['sequence_id'],row['input_type'],row['sequence_sha256'],int(row['length_nt'])))
                    counts[row['evaluation_stratum']]+=1
                    if sum(counts.values())%10000==0: db.commit()
            if db.execute('SELECT COUNT(*) FROM expected WHERE seen=0').fetchone()[0]:
                raise ValueError('Missing planned executions, including an incomplete paired arm')
            db.commit()
            batches=set()
            with (root/'manifest.csv').open(encoding='utf-8',newline='') as handle:
                reader=csv.DictReader(handle)
                if reader.fieldnames!=['prefix','type','fasta']: raise ValueError('Invalid viSUM manifest header')
                for row in reader:
                    batch=row['prefix']
                    path=(root/row['fasta']).resolve()
                    if batch in batches or path.parent!=root or row['fasta'] not in checks:
                        raise ValueError('Duplicate manifest prefix or unsafe/unhashed FASTA path')
                    batches.add(batch)
                    batch_count=0
                    for key,seq in fasta_records(path):
                        expected=db.execute('SELECT route,hash,length,seen FROM executions WHERE batch=? AND id=?',(batch,key)).fetchone()
                        if expected is None or expected[3] or expected[:3]!=(row['type'],hashlib.sha256(seq.encode('ascii')).hexdigest(),len(seq)):
                            raise ValueError('FASTA does not match execution map exactly')
                        db.execute('UPDATE executions SET seen=1 WHERE batch=? AND id=?',(batch,key))
                        batch_count+=1
                        if batch_count%10000==0: db.commit()
                    if batch_count!=report['batches'].get(batch): raise ValueError('Batch count mismatch')
            if batches!=set(report['batches']) or db.execute('SELECT COUNT(*) FROM executions WHERE seen=0').fetchone()[0]:
                raise ValueError('Incomplete manifest/FASTA execution coverage')
            emitted=db.execute('SELECT COUNT(*) FROM executions').fetchone()[0]
            if counts['unique_inputs']!=report['unique_input_sequences'] or emitted!=report['planned_execution_sequences']:
                raise ValueError('Unique versus execution counts do not reconcile')
        finally: db.close()
    result=dict(status='PREPARATION_VERIFIED',unique_inputs=counts['unique_inputs'],execution_sequences=emitted,
                strata={k:v for k,v in counts.items() if k!='unique_inputs'},batches=len(batches),
                note='Unresolved DNA and RNA are paired executions, not additional independent observations. No pipeline launched.')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    verify(p.parse_args().directory)
