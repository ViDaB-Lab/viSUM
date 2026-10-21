#!/usr/bin/env python3
"""Shared read-only file contracts for prepared viSUM assets.

Contracts follow genomad/database.py (v1.12.0) and vcontact3 databases.py,
profiles.py, config.py (d57ae81e0be5ab57a8a14692a4cb8ab255ce2d52).
This checks completeness/readability, not cryptographic authenticity.
Setup and analysis preflight use these same structural requirements. Setup
retains its additional native-tool, checksum and semantic validation.
"""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess


TOOLS = ('genomad', 'virsorter2', 'cenotetaker3', 'deep6', 'deepmicroclass2',
         'virbot', 'gianthunter', 'checkv', 'vicat', 'vicat_nonviral', 'vitap', 'vcontact3')


def prepared_directory(path, allow_installing=False):
    root = Path(path).resolve(strict=True)
    if not root.is_dir():
        raise ValueError('Expected a prepared directory: {}'.format(root))
    marker = root / '.visum_installing'
    if not allow_installing and (marker.exists() or marker.is_symlink()):
        raise ValueError('Installation is incomplete or in progress: {}. '
                         'Finish setup before analysis; the marker was not removed.'.format(marker))
    return root


def required(root, names):
    for name in names:
        nonempty(root / name)


def vicat_files(root):
    """Select a complete supported layout, not a mixture of two layouts."""
    for diamond_dir, taxonomy_dir in (('', ''), ('diamond', 'taxonomy')):
        diamond = root / diamond_dir / 'IMGVR5_UViG_representatives.dmnd'
        taxonomy = root / taxonomy_dir / 'IMGVR5_UViG.vicat_taxonomy_lookup.parquet'
        if all(p.is_file() and p.stat().st_size for p in (diamond, taxonomy)):
            return nonempty(diamond), nonempty(taxonomy)
    raise ValueError('Missing complete viCAT runtime or development layout: {}'.format(root))


def deep6_installation(root):
    required(root, ['Master/deep6.py', 'Master/deep6_functions.py'])


def deep6_models(directory):
    models = prepared_directory(directory)
    for length in (250, 500, 1000, 1500):
        files = [p for p in models.glob('model_{}_*.h5'.format(length)) if p.is_file() and p.stat().st_size]
        if len(files) != 1:
            raise ValueError('Expected exactly one model_{}_*.h5 in {}'.format(length, models))
        with files[0].open('rb') as handle:
            if handle.read(8) != b'\x89HDF\r\n\x1a\n':
                raise ValueError('Invalid HDF5 model: {}'.format(files[0]))


def validate_asset(tool, spec, allow_installing=False, runtime=False):
    """No downloads, repair, writes, imports of tool runtimes or full-file hashes.

    allow_installing is internal to installers validating their own build while
    holding the destination lock. Analysis never enables it.
    """
    root = prepared_directory(spec['path'], allow_installing)
    if tool == 'genomad':
        genomad(root, runtime)
    elif tool == 'virsorter2':
        required(root, ['hmm/viral/combined.hmm'] + ['group/{}/model'.format(g) for g in
                 ('dsDNAphage', 'ssDNA', 'NCLDV', 'RNA', 'lavidaviridae')])
        histories = list(root.glob('conda_envs/**/conda-meta/history'))
        if not histories:
            raise ValueError('Missing prepared VirSorter2 Conda environments: {}'.format(root / 'conda_envs'))
        for history in histories:
            nonempty(history)
    elif tool == 'cenotetaker3':
        version = spec['hmm_version']
        required(root, ['hmmscan_DBs/{}/{}.h3m'.format(version, h) for h in
                 ('Virion_HMMs', 'DNA_rep_HMMs', 'RDRP_HMMs', 'Useful_Annotation_HMMs', 'phrogs_for_ct')])
        for name in ('ct3_hallmark.taxDB', 'refseq_virus_prot_taxDB', 'CDD'):
            required(root, ['mmseqs_DBs/{}'.format(name), 'mmseqs_DBs/{}.dbtype'.format(name)])
            if name != 'CDD':
                required(root, ['mmseqs_DBs/{}_mapping'.format(name), 'mmseqs_DBs/{}_taxonomy'.format(name)])
        required(root, ['viral_cdds_and_pfams_191028.txt'])
    elif tool == 'deep6':
        deep6_installation(root)
        deep6_models(spec.get('models') or root / 'Models')
    elif tool == 'deepmicroclass2':
        required(root, ['predict.py', 'models.py', 'utils.py', 'selftest.py', 'model/300/DeepMicroClass-best.ckpt'] +
                 ['model/8-class/DeepMicroClass-best-{}-8class.ckpt'.format(n) for n in (500, 1000, 3000)])
    elif tool == 'virbot':
        names = ['virbot/VirBot.py', 'virbot/__init__.py'] + ['virbot/data/ref/' + name for name in
                 ('VirBot.hmm', 'VirBot.dmnd', 'VirBot_hmm_threshold.txt', 'VirBot_hmm_taxa_full.txt', 'VirBot_RNAvirus_acc.txt')]
        required(root, names)
        for name in names:
            with (root / name).open('rb') as handle:
                if handle.read(100).startswith(b'version https://git-lfs.github.com/spec/'):
                    raise ValueError('Unresolved Git-LFS pointer: {}'.format(root / name))
    elif tool == 'gianthunter':
        required(root, ['database.dmnd', 'names.csv', 'nodes.csv', 'ProkaryoticGroup.pkl', 'proteins.csv',
                        'RefVirus.dmnd', 'RefVirus.faa', 'RefVirus_anno.pkl', 'taxid.csv', 'transformer.pth'])
    elif tool == 'checkv':
        required(root, ['genome_db/checkv_reps.' + ext for ext in ('faa', 'fna', 'tsv', 'dmnd')] +
                 ['hmm_db/checkv_hmms.tsv', 'hmm_db/genome_lengths.tsv'])
        hmm = next((p for p in root.glob('hmm_db/checkv_hmms/*.hmm') if p.is_file() and p.stat().st_size), None)
        if hmm is None:
            raise ValueError('Missing CheckV HMM files: {}'.format(root / 'hmm_db/checkv_hmms'))
        nonempty(hmm)
    elif tool == 'vicat':
        vicat_files(root)
    elif tool == 'vicat_nonviral':
        required(root, ['vicat_nonviral.dmnd', 'vicat_nonviral_representatives.faa.gz',
                        'vicat_nonviral_representative_metadata.parquet', 'vicat_nonviral_database_metadata.tsv',
                        'vicat_nonviral_cluster_summary.tsv', 'SHA256SUMS', '.visum_db_complete'])
    elif tool == 'vitap':
        gffs = [p for p in root.glob('VMR_genome_*.gff') if p.is_file() and p.stat().st_size]
        if len(gffs) != 1:
            raise ValueError('Expected exactly one VMR_genome_*.gff in {}'.format(root))
        required(root, [gffs[0].with_suffix(ext).name for ext in ('.gff', '.fasta', '.faa', '.dmnd')] +
                 ['uniref90.dmnd', 'uniref90.accession2taxid', 'taxdmp/nodes.dmp', 'taxdmp/names.dmp'] +
                 [rank + '_genome.threshold' for rank in ('Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')])
        csvs = [p for p in root.glob('*.csv') if p.is_file() and p.stat().st_size]
        if not csvs:
            raise ValueError('Missing VITAP taxonomy CSV in {}'.format(root))
        for csv in csvs:
            nonempty(csv)
    elif tool == 'vcontact3':
        release = spec.get('release')
        if not release:
            manifests = sorted(root.glob('[0-9][0-9][0-9].json'))
            if not manifests:
                raise ValueError('Missing vConTACT3 release manifest: {}'.format(root))
            release = manifests[-1].stem
        vcontact3(root / (release + '.json'), spec.get('domain', 'both'), runtime)
    else:
        raise ValueError('Unknown asset tool: {}'.format(tool))
    return str(root)


def nonempty(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing or empty required file: {path}')
    with path.open('rb') as handle:
        handle.read(1)
    return path


def mmseqs(prefix, runtime=False, headers=False):
    prefix = Path(prefix)
    dbtype = nonempty(str(prefix) + '.dbtype')
    if dbtype.stat().st_size != 4:
        raise ValueError(f'Invalid MMseqs dbtype header: {dbtype}')
    nonempty(str(prefix) + '.index')
    if prefix.is_file() and prefix.stat().st_size:
        nonempty(prefix)
    else:
        parts = sorted((p for p in prefix.parent.glob(prefix.name + '.*')
                        if p.suffix[1:].isdigit()), key=lambda p: int(p.suffix[1:]))
        if not parts or [int(p.suffix[1:]) for p in parts] != list(range(len(parts))):
            raise ValueError(f'Missing MMseqs data or noncontiguous split files: {prefix}')
        for part in parts:
            nonempty(part)
    if runtime:
        subprocess.run(['mmseqs', 'dbtype', str(prefix)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if headers:
        mmseqs(str(prefix) + '_h', runtime)


def genomad(directory, runtime=False):
    root = Path(directory)
    version = float(nonempty(root / 'version.txt').read_text().strip())
    if not math.isfinite(version) or version <= 0:
        raise ValueError('Invalid geNomad database version')
    for name in ('genomad_db', 'genomad_mini_db', 'genomad_integrase_db'):
        mmseqs(root / name, runtime, headers=True)
    for name in ('nodes.dmp', 'names.dmp', 'genomad_marker_metadata.tsv'):
        nonempty(root / name)
    with (root / 'genomad_marker_metadata.tsv').open() as handle:
        header, first = handle.readline(), handle.readline()
        if '\t' not in header or '\t' not in first:
            raise ValueError('geNomad marker metadata is empty or malformed')


def local_asset(root, release, value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Reference asset must be a nonempty relative path')
    path = Path(value)
    if path.is_absolute() or '..' in path.parts or '\\' in value or ':' in value:
        raise ValueError(f'Unsafe/nonportable manifest path: {value}')
    resolved = (root / path).resolve()
    if not resolved.is_relative_to(release.resolve()):
        raise ValueError(f'Manifest asset is outside its release: {value}')
    return resolved


def vcontact3(manifest, domain='both', runtime=False):
    manifest = nonempty(manifest)
    if not re.fullmatch(r'[0-9]{3}', manifest.stem):
        raise ValueError('Expected a three-digit vConTACT3 release manifest')
    data = json.loads(manifest.read_text())
    release = manifest.parent / ('v' + manifest.stem)
    if not release.is_dir():
        raise ValueError(f'Missing release directory: {release}')
    info = data[manifest.stem]
    refseq = info['RefSeq']
    domains = ['prokaryotes', 'eukaryotes'] if domain == 'both' else [domain]
    for selected in domains:
        block = refseq['domains'][selected]
        for key in ('Genomes Report', 'gene2genome', 'proteins', 'nucleotides'):
            nonempty(local_asset(manifest.parent, release, block[key]))
        identities = block['identities']
        for identity in ('0.3', '0.4', '0.5', '0.6', '0.7'):
            cluster = identities[identity]
            for key in ('db', 'cluster'):
                mmseqs(local_asset(manifest.parent, release, cluster[key]), runtime,
                       headers=(key == 'db'))
            nonempty(local_asset(manifest.parent, release, cluster['cluster2genome']))
            # Also check other file references declared by this identity block.
            for key, value in cluster.items():
                if key not in ('db', 'cluster', 'cluster2genome') and isinstance(value, str) and '/' in value:
                    nonempty(local_asset(manifest.parent, release, value))
    mmseqs(local_asset(manifest.parent, release, info['VOGDB']['mmseqs_db']), runtime, headers=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tool', choices=TOOLS + ('deep6_installation', 'deep6_models'))
    parser.add_argument('path')
    parser.add_argument('--domain', choices=['both', 'prokaryotes', 'eukaryotes'], default='both')
    parser.add_argument('--check-mmseqs', action='store_true')
    parser.add_argument('--models')
    parser.add_argument('--hmm-version', default='v3.1.1')
    parser.add_argument('--allow-installing', action='store_true', help='Installer-only: validate under the destination lock')
    args = parser.parse_args()
    try:
        if args.tool == 'deep6_models':
            deep6_models(args.path)
        elif args.tool == 'deep6_installation':
            deep6_installation(prepared_directory(args.path))
        else:
            path = Path(args.path)
            spec = dict(path=str(path), domain=args.domain, models=args.models, hmm_version=args.hmm_version)
            if args.tool == 'vcontact3' and path.suffix == '.json':
                spec.update(path=str(path.parent), release=path.stem)
            validate_asset(args.tool, spec, args.allow_installing, args.check_mmseqs)
    except (OSError, ValueError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        raise SystemExit(f'ERROR: {args.tool} database validation failed: {error}') from error


if __name__ == '__main__':
    main()
