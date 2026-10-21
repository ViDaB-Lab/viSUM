#!/usr/bin/env python3
"""Check vConTACT3's Parquet runtime with a tiny in-memory table.

No reference files, sample data, environments or scientific settings are changed.
The default pandas reader is exercised just as in vConTACT3's reference loading.
"""
import importlib
import io
import sys


def check_parquet_runtime():
    try:
        arrow = importlib.import_module('pyarrow')
        pandas = importlib.import_module('pandas')
    except ImportError as exc:
        raise RuntimeError('vConTACT3 requires a working pandas/pyarrow Parquet runtime; '
                           'repair or recreate only the selected vConTACT3 environment.') from exc
    expected = pandas.DataFrame({
        'gene_id': ['test_gene_1', 'test_gene_2'],
        'genome_id': ['test_genome_1', 'test_genome_1'],
        'length': [42, 84],
        'optional_label': ['test', None],
    })
    buffer = io.BytesIO()
    expected.to_parquet(buffer, engine='pyarrow', index=False)
    buffer.seek(0)
    observed = pandas.read_parquet(buffer)
    pandas.testing.assert_frame_equal(observed, expected)
    return {'pandas': pandas.__version__, 'pyarrow': arrow.__version__}


def main():
    try:
        versions = check_parquet_runtime()
    except Exception as exc:
        print('ERROR: vConTACT3 Parquet runtime check failed: ' + str(exc), file=sys.stderr)
        return 1
    print('PASS: vConTACT3 Parquet round-trip; pandas={pandas}, pyarrow={pyarrow}'.format(**versions))
    return 0


if __name__ == '__main__':
    sys.exit(main())
