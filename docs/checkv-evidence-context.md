# CheckV evidence when completeness is undetermined

CheckV completeness estimation and gene-content evidence are separate outputs.
The adapter now preserves non-provirus records even when `checkv_quality` is
`Not-determined`. Measured viral/host gene counts and warnings remain available to
viHARMONY. An undetermined quality category is not itself negative evidence.

Existing interpretation rules remain unchanged:

- No informative viral or host genes: unclassified, weak evidence.
- Host genes but no viral genes: qualified cellular evidence.
- Viral genes with undetermined quality: weak viral evidence, with host context
  retained when present.
- Previously emitted records, including provirus records, are unchanged.

This change restores evidence delivery; it does not introduce a new whole-parent
rejection threshold or change CT3 discovery votes. Existing harmonization rules
decide how the added context affects retention or review status.

## Audited RNA-only conflict exception

The harmonizer now permits an intact RNA candidate to remain primary despite a
CheckV-only cellular conflict when Deep6, VirBot and viCAT all supply qualified or
strong support, there is no plasmid/mobile conflict, and no viCAT locus is cellular
or ambiguous. At least one viral-supported protein must be >=100 aa, bitscore
>=100, identity >=50%, query and subject coverage >=80%, and its coding interval
must span >=70% of the input. These are the balanced exploratory thresholds, not
an independently optimized or calibrated classifier.

The original CheckV conflict remains in `cellular_conflict_tools`. Added metadata
fields `rna_checkv_conflict_exception` and `rna_checkv_exception_supporting_loci`
record application of the exception. The harmonizer manifest records thresholds
and affected IDs. Other downstream safeguards can still demote a candidate.
No DNA sequence or extracted region is eligible. Missing supporting loci cannot
trigger rescue. The direct harmonizer comparison switch
`--disable-rna-checkv-exception` disables this exception; the pipeline enables it
by default and supplies RNA locus evidence independently of the pair-floor option.

`bin/replay_checkv_context.py` tests this change with saved boundaries, sequences,
and taxonomy held fixed. It first reproduces four saved baseline tables, checks
that previously emitted CheckV rows are unchanged, and then substitutes only the
restored CheckV evidence. Use a new output directory:

```bash
python bin/replay_checkv_context.py \
  --results-root /path/to/saved/results \
  --ictv-msl assets/ICTV_VMR_MSL41.csv \
  --output-dir /path/to/new/checkv_context_review
```

If FASTAs were packaged separately, supply `--fasta-root` pointing to the directory
containing the corresponding sample result folders. This isolated replay does not
replace an end-to-end test combining the new trimming default with restored
CheckV context. Sequence-changing downstream stages must be regenerated there.
