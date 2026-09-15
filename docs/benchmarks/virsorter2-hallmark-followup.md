# VirSorter2 hallmark-only policy: implementation and saved-evidence replay

## Final release policy: review-only, not a consensus vote

The qualified-vote experiment below is historical and is **not the chosen release policy**. Following that experiment, the adapter now emits `evidence_strength=review` while retaining `classification=virus`, the native hallmark evidence basis, missing score and unlocalized coordinates. This preserves the native prediction for audit and standalone comparisons without making it a discovery-routing or primary-retention vote.

The discovery gate and viHARMONY explicitly recognize the hallmark-only evidence basis, even if an older intermediate file labels it qualified/strong. viHARMONY's full evidence audit labels applicable rows `review_only_no_consensus_vote`. An input without any other qualifying evidence does not enter downstream analysis just because of this review-only prediction; its prediction is still retained in the published VirSorter2 evidence table.

Replay verification uses `--require-unchanged`: each regenerated discovery-gate table and four decision tables must exactly match the original archived run. The original promotion-policy result and its one newly gated input below must not be confused with this finalized policy.

### Finalized-policy verification: PASS

All ten panels passed the replay and an additional completed-output verification:

- 10 discovery-routing tables byte-identical to the archive.
- 40 decision tables byte-identical to the archive (primary candidates, final metadata including taxonomy, provisional metadata, and per-input disposition).
- 30 primary/provisional FASTA comparisons byte-identical between baseline and finalized-policy replay.
- 147 native hallmark-only records preserved with review status and missing scores/boundaries; 146 appear explicitly as review-only in the harmonizer audit. The remaining record belongs to the input that does not advance to harmonization and remains preserved in its standardized VirSorter2 evidence table.
- Zero newly gated inputs and zero changes in primary parent retention. The earlier downstream-data gap is eliminated under this finalized policy.
- Focused tests: 77 passed plus two subtests. No discovery programs or downstream biological tools were rerun; this was a saved-evidence regression, not a Nextflow installation smoke test.

| Endpoint | Original baseline | Final review-only policy |
| --- | --- | --- |
| DNA positives | 477/539 | 477/539 |
| RNA positives | 788/874 | 788/874 |
| DNA labelled negatives | 101/2,390 | 101/2,390 |
| Clean-RNA labelled negatives | 27/1,600 | 27/1,600 |
| RNA mobile challenge | 0/200 | 0/200 |

The two additional clean-RNA retentions from the rejected qualified-vote trial are no longer primary retentions. The next step is the clean-install smoke test, not another discovery benchmark. Source changes remain local and have not been pushed or released.

```bash
python bin/replay_virsorter2_hallmarks.py \
  --review-root /path/to/extracted/consolidated_review \
  --fasta-root /path/to/previous/replay/fasta \
  --outdir /path/to/new/review_only_replay \
  --require-unchanged
```

`--fasta-root` reuses the previously extracted FASTAs. The script can alternatively use `--archive` as described below. No discovery tool or database rebuild is involved.

---

## Historical experiment: qualified primary-retention votes

Date: 9 September 2026. This is a follow-up to the [original consolidated review](preliminary-2026-09.md), not a replacement for its pre-change measurements.

## Implemented policy

- Preserve native unscored `lt2gene` predictions with at least one hallmark as `qualified` viral evidence.
- Preserve the distinction through `strength_basis=virsorter2_hallmark_only`, `score_type=virsorter2_hallmark_only`, `score=NA`, and empty coordinates. No classifier group is promoted to taxonomy.
- A hallmark-only call by itself remains provisional under the existing origin policy. It can participate in qualified multi-tool support.
- Unlocalized short-hallmark evidence applies to the intact parent, not arbitrary extracted child regions. It does not invent or initiate a boundary.
- For RNA-floor exemption eligibility only, disregard VirSorter2 when all its applicable viral support is hallmark-only. A Deep6+viCAT pair cannot escape the optional floor by acquiring this extra tool label. Ordinary scored VirSorter2 support retains its existing behavior.
- Both affected Nextflow module command bodies carry changed policy markers, so `-resume` can invalidate the adapter/harmonizer tasks even though Python implementation files are referenced by project path. Discovery executable/model commands were not changed.
- The existing server benchmark helper now uses the `visum` launcher, including its `NEXTFLOW_BIN`/sibling-executable resolution.

## Quantified effect

The replay re-ran the adapter and viHARMONY, holding saved discovery-gate, region selection and downstream evidence fixed. It did **not** run discovery tools, CheckV, TEsorter, VITAP, vConTACT3 or model/database construction. Before accepting each comparison, it reproduced the archived baseline byte-for-byte for `final_metadata.tsv`, `database_candidates.tsv`, `provisional_metadata.tsv` and `sequence_disposition.tsv`: all four tables for all ten panels.

| Panel | Restored native VS2 records | Primary parents before | Primary parents after, fixed downstream evidence | Newly gated parents outside replay |
| --- | --- | --- | --- | --- |
| DNA viral | 1 | 477/539 | 477/539 | 0 |
| RNA viral | 141 | 788/874 | 788/874 | 0 |
| DNA cellular | 0 | 15/500 | 15/500 | 0 |
| DNA mitochondrial | 0 | 11/500 | 11/500 | 0 |
| DNA plasmid | 0 | 67/500 | 67/500 | 0 |
| DNA plastid | 0 | 1/500 | 1/500 | 0 |
| DNA retroelement | 0 | 7/390 | 7/390 | 0 |
| RNA clean-confirmed | 4 | 20/1,200 | 22/1,200 | 1 |
| RNA clean-broad | 0 | 7/400 | 7/400 | 0 |
| RNA mobile challenge | 1 | 0/200 | 0/200 | 0 |

Thus, **no primary positive sensitivity gain was observed**. Fixed-candidate clean-RNA FPR rises from 27/1,600 (1.6875%) to 29/1,600 (1.8125%). DNA primary detection/FPR and mobile-challenge retention remain unchanged. Restoring legitimate native evidence does not necessarily increase final sensitivity: many affected positives already survive through other tools; others remain below the integrated support/conflict requirements.

The two promoted labelled negatives are:

| Input | Source annotation | New qualified support |
| --- | --- | --- |
| `RNANEG_CCONF_000876` | *Hydra vulgaris* collagen alpha-1(IV) chain-like mRNA, `NM_001309661.1` | Deep6 + VirSorter2 hallmark-only |
| `RNANEG_CCONF_001067` | Ig-like domain-containing protein CDS, `WP_051658902.1` | Deep6 + VirSorter2 hallmark-only |

Neither has strong support or the Deep6+viCAT pair, so these promotions are **not a bypass of the RNA homology floor**. They follow the existing two-qualified-tools rule. Both remain false positives under the benchmark labels; this replay does not independently establish their biological origin.

### One remaining end-to-end uncertainty

`RNANEG_CLEAN_CONFIRMED__c000653` becomes eligible for discovery routing but was absent from the saved refined-candidate set. Its downstream results do not exist in the supplied archive. It is reported in `newly_gated.tsv`, not silently scored as a completed new negative. The 29/1,600 result is therefore explicitly a **fixed-candidate replay**, not a fully updated end-to-end benchmark.

The saved positive panels have no newly gated parents. The adapter was checked against native outputs; 147 records are restored in total. The focused policy/gate/refinement/workflow suite passed 74 tests plus two subtests.

## Reproduce without rerunning discovery tools

Use a new output directory. Either provide an extracted review tree with its original `prep/*.normalized.fasta` and `refinement/*.refined_candidates.fasta`, or let the script extract just those FASTAs from the supplied archive:

```bash
python bin/replay_virsorter2_hallmarks.py \
  --review-root /path/to/extracted/consolidated_review \
  --archive /path/to/consolidated_review_20260909_115020.tar.gz \
  --outdir /path/to/new/hallmark_replay
```

Outputs include `summary.tsv`, `changes.tsv`, `newly_gated.tsv`, `scope.json`, regenerated VirSorter2 evidence and baseline/updated viHARMONY results. The script refuses to overwrite an existing destination. This review-layout helper requires the archived all-tools configuration and fixed downstream files; it is not a general Nextflow replacement.

## Release decision

The adapter correction is semantically justified, but **this experiment does not support marketing it as improved sensitivity or specificity**. Keep its effect visible. Do not silently tune an additional rule against the two named negative examples and then describe the same benchmark as independent validation.

Before releasing this policy, perform a cached Nextflow resume on the server using the same parameters, database paths and work directory, writing to a new output directory. The source update must first be installed there. The revised adapter and downstream tasks may run; unchanged discovery tasks should be cached. Confirm that in the trace. This resolves the newly admitted negative and verifies real workflow integration without rebuilding the benchmark or rerunning every discovery program.

The existing `bin/benchmark_rna_pair.sh` is a **lab-specific** helper with the server's prior database/cache paths and resource settings; it can be used for this verification only on that configured server. It is not the public installation recipe. This implementation has not been pushed or released as part of the review.

viSUM's intended value remains broader detection plus standardized, auditable classification for reference curation and cross-study comparison. Confidence tiers are evidence summaries, not calibrated probabilities; novel gene-sharing groups are candidate groups, not automatically validated taxa. Accuracy on genuinely novel/non-model viruses needs reference-aware held-out evaluation. These qualifications distinguish what is already demonstrated from what the pipeline is designed to achieve; they do not demand perfect sensitivity or zero false positives.
