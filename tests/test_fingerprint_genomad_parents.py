import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import fingerprint_genomad_parents as target


class FingerprintTests(unittest.TestCase):
    def test_reverse_complement_and_uracil(self):
        self.assertEqual(target.fingerprints('aacu')[1], target.fingerprints('AGTT')[1])
        self.assertNotEqual(target.fingerprints('AACT')[1], target.fingerprints('AAAT')[1])

    def test_invalid_sequence(self):
        for value in ('', 'AC-G', 'ACTZ'):
            with self.assertRaises(ValueError):
                target.fingerprints(value)

    def test_skip_unrequested_large_records(self):
        self.assertEqual(list(target.selected_records(io.StringIO('>skip\nNNNN\n>keep note\nAC\nGT\n'), {'keep'})),
                         [('keep', 'keep note', 'ACGT')])

    def test_archive_reconciliation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            truth = root / 'truth.tsv'
            truth.write_text('sequence_id\tsource_class\nv|x|fragment_1\tvirus\nc|fragment_1\tchromosome\n')
            archive = root / 'author.zip'
            with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as output:
                output.writestr(target.MEMBER, gzip.compress(b'>c\nAAA\n>v|x original\nACGT\n'))
            self.assertEqual(target.run(truth, archive, root / 'out'), 0)
            summary = json.loads((root / 'out/summary.json').read_text())
            self.assertEqual(summary['fingerprinted_parents'], 1)
            with gzip.open(root / 'out/viral_parent_fingerprints.tsv.gz', 'rt') as handle:
                row = next(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(row['parent_sequence_id'], 'v|x')
            self.assertEqual(row['original_header'], 'v|x original')
            self.assertTrue((root / 'out/routing_evidence.tar.gz').is_file())
            with self.assertRaises(FileExistsError):
                target.run(truth, archive, root / 'out')

    def test_missing_and_duplicate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            truth = root / 'truth.tsv'
            truth.write_text('sequence_id\tsource_class\nv|fragment_1\tvirus\n')
            for name, fasta in [('missing', b'>c\nACGT\n'), ('duplicate', b'>v\nACGT\n>v\nACGT\n')]:
                archive = root / (name + '.zip')
                with zipfile.ZipFile(archive, 'w') as output:
                    output.writestr(target.MEMBER, gzip.compress(fasta))
                if name == 'missing':
                    self.assertEqual(target.run(truth, archive, root / name), 2)
                else:
                    with self.assertRaisesRegex(ValueError, 'Duplicate'):
                        target.run(truth, archive, root / name)


if __name__ == '__main__':
    unittest.main()
