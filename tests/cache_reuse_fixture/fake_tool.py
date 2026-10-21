"""Tiny, offline tool substitutes for cache tests; never classify sequences."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time

HARMONY_SUFFIXES = (
    'final.normalized.fasta', 'final.original_ids.fasta', 'final_metadata.tsv',
    'review_queue.tsv', 'sequence_disposition.tsv', 'sequence_map.tsv',
    'harmonizer_manifest.json', 'database_candidates.fasta', 'database_candidates.tsv',
    'all_candidates.fasta', 'provisional.normalized.fasta',
    'provisional.original_ids.fasta', 'provisional_metadata.tsv',
    'TE_candidates.fasta', 'TE_candidates.summary.tsv',
    'viral_candidates.fasta', 'viral_candidates.summary.tsv', 'taxonomy.tsv',
    'checkv.tsv', 'vcontact3.tsv', 'review.tsv',
)
HARMONY_PATHS = tuple(
    ('audit/legacy/' if suffix in {'final.normalized.fasta','final.original_ids.fasta',
        'database_candidates.fasta','database_candidates.tsv','all_candidates.fasta',
        'provisional.normalized.fasta','provisional.original_ids.fasta','provisional_metadata.tsv'} else
     'audit/' if suffix in {'final_metadata.tsv','review_queue.tsv','sequence_disposition.tsv',
        'sequence_map.tsv','harmonizer_manifest.json'} else
     'TE_candidates/' if suffix.startswith('TE_candidates.') else
     'quality/' if suffix=='checkv.tsv' else 'groups/' if suffix=='vcontact3.tsv' else
     'review/' if suffix=='review.tsv' else '') + '{sample}.' + suffix
    for suffix in HARMONY_SUFFIXES)

# Consumer-valid empty geNomad tables. Keep these explicit so the real adapter
# tests detect schema drift rather than silently borrowing its required fields.
GENOMAD_TABLE_COLUMNS = {
    'virus_summary.tsv': ('seq_name', 'length', 'topology', 'coordinates', 'n_genes',
                          'genetic_code', 'virus_score', 'fdr', 'n_hallmarks',
                          'marker_enrichment', 'taxonomy'),
    'plasmid_summary.tsv': ('seq_name', 'length', 'topology', 'n_genes', 'genetic_code',
                            'plasmid_score', 'fdr', 'n_hallmarks', 'marker_enrichment'),
    'virus_genes.tsv': ('gene', 'uscg'),
}


def option(args, name):
    return args[args.index(name) + 1]


def record(tool, sample, **extra):
    root = Path(os.environ['VISUM_CACHE_TEST_ROOT'])
    # Separate files prevent concurrent tasks from interleaving event writes.
    event = root / 'events' / f'{time.time_ns()}-{os.getpid()}.pending'
    event.write_text(json.dumps(dict(tool=tool, sample=sample, cwd=str(Path.cwd()), **extra)))
    event.rename(event.with_suffix('.json'))
    return root


def main(argv=None):
    tool, *args = sys.argv[1:] if argv is None else argv
    if tool == 'genomad':
        if args == ['--version']:
            print('genomad cache-test-fixture')
            return 0
        assert args[0] == 'end-to-end', args
        sample = Path(args[-3]).stem
        record('RUN_GENOMAD', sample)
        out = Path(args[-2]) / f'{sample}_summary'
        out.mkdir(parents=True)
        for suffix in ('virus_summary.tsv', 'virus.fna', 'virus_genes.tsv',
                       'virus_proteins.faa', 'plasmid_summary.tsv'):
            columns = GENOMAD_TABLE_COLUMNS.get(suffix)
            (out / f'{sample}_{suffix}').write_text('\t'.join(columns) + '\n' if columns else '')
    elif tool == 'virsorter':
        assert args[0] == 'run', args
        sample = Path(option(args, '-i')).stem
        record('RUN_VIRSORTER2', sample)
        out = Path(option(args, '-w'))
        out.mkdir()
        for name in ('final-viral-score.tsv', 'final-viral-boundary.tsv', 'final-viral-combined.fa'):
            (out / name).write_text('seqname\n' if name.endswith('.tsv') else '')
    elif tool == 'diamond':
        assert args[0] == 'blastp', args
        out = Path(option(args, '--out'))
        process = 'RUN_VICAT_NONVIRAL_DIAMOND' if 'nonviral' in out.name else 'RUN_VICAT_DIAMOND'
        record(process, Path(option(args, '--query')).name.split('.orfs')[0])
        out.write_text('')
    elif tool == 'mmseqs':
        assert args[0] == 'dbtype', args
        assert Path(args[1] + '.dbtype').is_file(), args
        print('Aminoacid')
    elif tool in ('python3', 'python'):
        if not args or Path(args[0]).name != 'run_viharmony.py':
            return subprocess.call([sys.executable, *args])
        sample = option(args, '--sample-id')
        root = record('VIHARMONY', sample, workers=int(option(args, '--workers')))
        if os.environ.get('VISUM_TEST_REAL_HARMONY') == '1':
            return subprocess.call([sys.executable, *args])
        # Only this fixture holds tasks. Production viHARMONY is not edited.
        deadline = time.monotonic() + 60
        while (root / 'hold-harmony').exists():
            if time.monotonic() >= deadline:
                return 75
            time.sleep(0.05)
        for pattern in HARMONY_PATHS:
            output = Path(pattern.format(sample=sample))
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(f'cache-test-only:{sample}\n')
    else:
        raise ValueError(f'Unexpected fixture executable: {tool}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
