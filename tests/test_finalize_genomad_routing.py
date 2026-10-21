import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
import finalize_genomad_routing as mod


class FinalizeRoutingTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.parents={'A.1':dict(length_nt='10',fragment_count='2',canonical_sequence_sha256='abc')}
        self.package=self.root/'package'
        self.package.write_text('fixture')
        self.evidence=self.root/'evidence.jsonl'
        self.ncbi=self.root/'ncbi.json'
        self.ncbi.write_text('{}')

    def run_case(self,items,record=None):
        self.evidence.write_text(json.dumps(dict(parent='A.1',evidence=items))+'\n' if items else '')
        if record:self.ncbi.write_text(json.dumps({'A.1':{'record':record}}))
        with patch.object(mod,'load_parents',return_value=(self.parents,{'truth_sha256':'truth'})),contextlib.redirect_stdout(io.StringIO()):
            return mod.finalize(self.package,self.evidence,self.ncbi,self.root/'out')

    def source(self,route='rna',exact=False):
        return dict(route=route,source='https://example.org',
            record='exact_full_length_canonical_sha256=abc' if exact else 'A.1',
            basis='exact_identity_to_published_RNA_reference' if exact else 'exact_accession_version_and_length_with_NCBI_record_molecule_type')

    def test_dna_record_is_deferred(self):
        self.assertEqual(self.run_case([self.source('dna')])['fragment_counts'],{'unresolved':2})

    def test_exact_rna_is_accepted(self):
        self.assertEqual(self.run_case([self.source(exact=True)])['fragment_counts'],{'rna':2})

    def test_conflict_remains_dual(self):
        self.assertEqual(self.run_case([self.source(exact=True),self.source('dna')])['fragment_counts'],{'unresolved':2})

    def test_genomic_rna_is_accepted(self):
        r=dict(accessionversion='A.1',moltype='rna',biomol='genomic',slen=10)
        self.assertEqual(self.run_case([self.source()],r)['fragment_counts'],{'rna':2})

    def test_mrna_not_assumed_rna_viral_genome(self):
        r=dict(accessionversion='A.1',moltype='rna',biomol='mRNA',slen=10)
        self.assertEqual(self.run_case([self.source()],r)['fragment_counts'],{'unresolved':2})

    def test_missing_source_remains_dual(self):
        self.assertEqual(self.run_case([])['fragment_counts'],{'unresolved':2})


if __name__=='__main__':unittest.main()
