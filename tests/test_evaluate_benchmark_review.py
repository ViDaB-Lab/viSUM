import importlib.util
from pathlib import Path

import pytest

spec=importlib.util.spec_from_file_location('benchmark_review',Path(__file__).parents[1]/'bin/evaluate_benchmark_review.py')
review=importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)

def test_native_virsorter_includes_hallmark_only_and_deduplicates_regions():
    rows=[{'seqname':name} for name in ['S__c000001||full','S__c000002||lt2gene','S__c000003||0_partial','S__c000003||1_partial']]
    assert review.native_virsorter2_calls(rows,'S')=={'S__c000001','S__c000002','S__c000003'}

def test_native_virsorter_rejects_foreign_or_malformed_identifiers():
    for name in ['OTHER__c000001||full','S__c000001garbage','S__c000001||unknown']:
        with pytest.raises(ValueError):review.native_virsorter2_calls([{'seqname':name}],'S')
