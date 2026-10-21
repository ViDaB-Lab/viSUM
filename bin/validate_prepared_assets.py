#!/usr/bin/env python3
"""Read-only structural preflight. Runtime checks belong to setup/smoke testing."""

import argparse
import base64
import json
from pathlib import Path
import re

from validate_reference_database import nonempty, validate_asset

SETUP_HINT = 'Run Nextflow with --setup in a separate run, or supply paths to complete prepared assets.'


def validate(spec):
    samples = spec['samples']
    if not samples:
        raise ValueError('No input samples were supplied')
    prefixes = set()
    for sample in samples:
        prefix = sample['prefix']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', prefix):
            raise ValueError(f'Invalid sample prefix: {prefix}')
        if prefix in prefixes:
            raise ValueError(f'Duplicate sample prefix: {prefix}')
        prefixes.add(prefix)
        if sample['type'] not in ('dna', 'rna'):
            raise ValueError(f'Invalid input type: {sample["type"]}')
        nonempty(sample['fasta'])
    types = {sample['type'] for sample in samples}
    results, errors = {}, []
    for tool, asset in spec['assets'].items():
        if tool in ('deep6', 'virbot') and 'rna' not in types:
            continue
        if tool in ('deepmicroclass2', 'gianthunter') and 'dna' not in types:
            continue
        try:
            results[tool] = validate_asset(tool, asset)
        except (OSError, ValueError, KeyError, TypeError) as error:
            errors.append(f'{tool}: {error}')
    if errors:
        raise ValueError('Prepared assets are missing or incomplete:\n' + '\n'.join(errors) + '\n' + SETUP_HINT)
    return {'validation': 'passed', 'level': 'read-only-structural',
            'runtime_validation': 'not-performed', 'sample_count': len(samples), 'assets': results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec-base64', required=True)
    parser.add_argument('--report', required=True, type=Path)
    args = parser.parse_args()
    try:
        spec = json.loads(base64.b64decode(args.spec_base64, validate=True))
        report = validate(spec)
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise SystemExit(f'ERROR: viSUM preflight: {error}') from error
    args.report.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
