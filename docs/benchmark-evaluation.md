# Frozen alpha benchmark endpoints

Endpoint contract: `visum-alpha-endpoints-v1`. These definitions affect evaluation
only. They do not alter viHARMONY decisions, taxonomic votes, retention, or TE routing.

The counting unit is an original INPUT FRAGMENT in its execution batch, not an
extracted region and not the longer source genome from which a fragment arose.
The unique execution key is `(benchmark_batch, sequence_id)`.

| Endpoint | Candidate records selected | Input-level interpretation |
| --- | --- | --- |
| `primary` | `candidate_destination=viral_candidates` AND `viral_decision` is `retained_viral` or `retained_provirus` | At least one primary viral record from the input |
| `viral_candidates` | Every record with viral destination, including provisional decisions | Broader provisional-inclusive viral export |
| `TE_candidates` | Every record with TE destination | At least one TE export; NOT a viral detection |
| `TE_only` | TE-export input IDs minus broader viral-export input IDs | TE export with no viral export from that same input |

`viral_confidence=supported` alone does NOT imply primary retention: a provisional
boundary decision can retain that confidence. Neither a review flag nor taxonomy
on an excluded record creates a positive call. Multiple exported regions from one
input count once for each applicable endpoint. A mixed parent can contribute to
both viral and TE export counts, but never to `TE_only`. Intact-fragment recovery
and boundary accuracy remain separate outcomes from input-level detection.

Keep the original `virus`, `chromosome`, and `plasmid` truth labels unchanged.
For primary and broader viral endpoints, report original-label sensitivity and
FPR separately. TE endpoints produce export rates, not viral sensitivity/FPR.
Contextual review of possible proviruses/phage-plasmids is a secondary analysis;
it must not silently change denominators or remove inconvenient labelled negatives.
These endpoint choices are fixed before the full benchmark is examined.

## Prepared geNomad-derived collection

The approved preparation plan contains 1,445,955 original fragments and 1,960,294
planned route executions. The three strata are `source_routed`, `unresolved_dna`,
and `unresolved_rna`. Never pool unresolved paired arms, take their union, or
choose each sequence's better-performing route. Keep source classes and DNA/RNA
routes separate. Missing results are errors, not negative predictions.
The genomic chromosome/plasmid controls do not provide an RNA-negative transcript
panel; clean-RNA FPR requires the separate RNA controls, not a pooled DNA denominator.
The `author_weight_split_*` fields remain training-weight metadata; they are not
automatically evaluation weights or proof of training holdout. Current pretrained
tools may overlap this collection; disclose reference/training overlap rather
than claiming independent unseen validation or exact paper cross-validation.

The preparation utilities preserve IDs, labels, sequence hashes and routing.
Existing prepared inputs should be verified, not regenerated unnecessarily:

```bash
python3 -B /path/to/candidate/bin/verify_genomad_preparation.py /path/to/prepared_inputs
```

After a complete analysis, the new dependency-free evaluator reads both flat and
nested manifests, verifies all result artifacts, and writes per-execution endpoint
flags, stratified original-label metrics, and input/result provenance:

```bash
python3 -B /path/to/candidate/bin/evaluate_prepared_benchmark.py \
  --prepared /path/to/prepared_inputs --results /path/to/run/results \
  --outdir /path/to/NEW_evaluation_directory
```

This does not launch tools or reread the multi-GB input FASTAs. It verifies the
prepared metadata checksums, execution/route coverage, and result checksums;
the separate preparation verifier checks the FASTAs. It uses disk-backed SQLite
for truth/execution joins and reads results one batch at a time. Native-tool
comparative metrics, rank-specific taxonomy accuracy, boundary accuracy and
resource-use analysis are separate publication analyses, NOT established by
these binary export metrics. Preserve native per-tool outputs and their exact
configurations for those comparisons. Do not assume integrated adapters reproduce
every native call or that specialist tools apply to every DNA/RNA subgroup.

## Historical development readers

`evaluate_benchmark_review.py`, `summarize_resumed_dna_benchmark.py`, and
`verify_rna_pair_benchmark.py` now resolve artifacts from the manifest rather than
assuming flat filenames. The first reader is specifically for archived development
ON/OFF panels and explicitly permits omitted archive artifacts; it reports every
omission. It never fabricates an absent OFF run from demotion annotations. Its
legacy `viSUM_ON` / `viSUM_OFF` names mean PRIMARY detection; additional export
endpoints are separate. Fixed-count RNA regression checks are historical tests,
not expected scores for the changed alpha TE policy or a new benchmark.

Native comparative baselines and full publication figures still require an
explicit protocol and actual full runs. No benchmark score is claimed here.
