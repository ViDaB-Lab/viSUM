# Preliminary consolidated benchmark review — 9 September 2026

Implementation follow-up: the [VirSorter2 hallmark policy and saved-evidence replay](virsorter2-hallmark-followup.md) documents a rejected qualified-vote trial and the finalized **review-only** policy. The finalized policy preserves native predictions without routing/primary votes; all ten panels reproduce the original decisions, with no extra downstream input. Numbers in this original review remain valid for primary retention under the finalized policy.

## Verdict and scope

viSUM demonstrates a useful **sensitivity–specificity tradeoff**, not universal superiority over its constituent tools. With the optional RNA homology floor enabled, it retains 477/539 DNA positives and 788/874 RNA positives. geNomad retains 452/539 and 646/874, respectively. viSUM also retains more labelled negative controls. These are development benchmarks, not independent validation of a finished classifier.

This review evaluates `consolidated_review_20260909_115020.tar.gz`, particularly `rna_pair_floor_20260909_113555` and `rna_pair_OFF_20260909_114659`. The source snapshot reviewed locally is `4f8f377`; archive configuration and per-tool metadata are the run-specific authority. Numbers below describe the archived runs, **not the reduced default configuration**. No production decision policy was changed during this review.

The important new implementation finding is that viSUM excludes VirSorter2's unscored `lt2gene` hallmark predictions. The standalone comparison below restores these native predictions. This issue must be addressed explicitly before presenting a definitive integration benchmark.

## Counting and fairness

- Unit: original input sequence, counted once even when several extracted regions survive. A positive is detected if any primary region or intact sequence from it is retained. This is **parent-level detection**, not recovery of the entire viral genome or accuracy of extracted boundaries.
- viSUM calls: `database_candidates.tsv`; provisional/review-only outputs are not primary detections. Discovery-tool calls: configured pre-gate standardized viral evidence, with VirSorter2 native predictions reported separately and used as its headline comparator.
- Sensitivity = detected labelled positives / all labelled positives. FPR = retained labelled negatives / all labelled negatives. Inputs below a tool's supported length remain in this operational denominator; scope-specific performance is necessary before ranking specialist tools.
- Unsupported DNA/RNA combinations are not zero sensitivity. Deep6/VirBot are RNA-only here; DeepMicroClass2/GiantHunter are DNA-only. GiantHunter targets giant DNA viruses, so its all-DNA sensitivity is not a fair assessment of its intended niche.
- CheckV, TEsorter, VITAP and vConTACT3 are downstream refinement/interpretation tools in this workflow, not comparable standalone discovery methods.
- The RNA mobile-element panel is a separate challenge set, not part of the headline clean-RNA FPR.
- Exact counts are more defensible here than population precision or unqualified confidence intervals: sequences/segments can share parents, references and evolutionary relationships. There is no demonstrated training/reference exclusion across every tool. Species holdout against the viCAT nonviral database is not universal benchmark independence.

## DNA: detection and labelled-negative retention

Each cell is count / denominator (percent). Positive column is sensitivity; remaining columns are label-based FPR.

| Method | Viral | Cellular | Mitochondrial | Plasmid | Plastid | Retroelement |
| --- | --- | --- | --- | --- | --- | --- |
| geNomad | 452/539 (83.86%) | 3/500 (0.60%) | 1/500 (0.20%) | 18/500 (3.60%) | 1/500 (0.20%) | 15/390 (3.85%) |
| VirSorter2 native | 440/539 (81.63%) | 12/500 (2.40%) | 24/500 (4.80%) | 51/500 (10.20%) | 21/500 (4.20%) | 7/390 (1.79%) |
| Cenote-Taker3 | 453/539 (84.04%) | 27/500 (5.40%) | 29/500 (5.80%) | 91/500 (18.20%) | 7/500 (1.40%) | 6/390 (1.54%) |
| DeepMicroClass2 | 155/539 (28.76%) | 4/500 (0.80%) | 0/500 (0.00%) | 0/500 (0.00%) | 0/500 (0.00%) | 92/390 (23.59%) |
| GiantHunter | 24/539 (4.45%) | 122/500 (24.40%) | 7/500 (1.40%) | 15/500 (3.00%) | 19/500 (3.80%) | 13/390 (3.33%) |
| viCAT | 524/539 (97.22%) | 102/500 (20.40%) | 180/500 (36.00%) | 336/500 (67.20%) | 57/500 (11.40%) | 43/390 (11.03%) |
| viSUM | 477/539 (88.50%) | 15/500 (3.00%) | 11/500 (2.20%) | 67/500 (13.40%) | 1/500 (0.20%) | 7/390 (1.79%) |

Across the five DNA negative panels, viSUM retains 101/2,390 (4.23%), versus geNomad 38/2,390 (1.59%), VirSorter2 115/2,390 (4.81%), Cenote-Taker3 160/2,390 (6.69%), and viCAT 718/2,390 (30.04%). The plasmid panel is the largest viSUM burden. The table does not establish whether any particular plasmid is actually a phage-plasmid.

## RNA: detection and labelled-negative retention

| Method | Viral | Clean-confirmed | Clean-broad | Mobile challenge |
| --- | --- | --- | --- | --- |
| geNomad | 646/874 (73.91%) | 16/1,200 (1.33%) | 2/400 (0.50%) | 0/200 (0.00%) |
| VirSorter2 native | 429/874 (49.08%) | 12/1,200 (1.00%) | 3/400 (0.75%) | 1/200 (0.50%) |
| Cenote-Taker3 | 566/874 (64.76%) | 9/1,200 (0.75%) | 2/400 (0.50%) | 1/200 (0.50%) |
| Deep6 | 747/874 (85.47%) | 303/1,200 (25.25%) | 127/400 (31.75%) | 65/200 (32.50%) |
| VirBot | 749/874 (85.70%) | 0/1,200 (0.00%) | 1/400 (0.25%) | 0/200 (0.00%) |
| viCAT | 833/874 (95.31%) | 152/1,200 (12.67%) | 34/400 (8.50%) | 35/200 (17.50%) |
| viSUM, RNA rule ON | 788/874 (90.16%) | 20/1,200 (1.67%) | 7/400 (1.75%) | 0/200 (0.00%) |

Combined clean-RNA FPR: viSUM 27/1,600 (1.69%), geNomad 18/1,600 (1.13%), VirSorter2 15/1,600 (0.94%), Cenote-Taker3 11/1,600 (0.69%), Deep6 430/1,600 (26.88%), VirBot 1/1,600 (0.06%), viCAT 186/1,600 (11.63%). Zero observations do not imply a population FPR of zero.

VirBot is a particularly important comparator: it achieves 85.70% sensitivity with only one clean-negative call. viSUM gains 39 positives net, but retains 26 more clean negatives net. There is no basis for calling viSUM strictly better for every RNA use case.

### Native-output versus integrated-evidence discrepancy

The current VirSorter2 standardizer drops unscored `lt2gene` records. Its integrated counts are DNA viral 439; RNA viral 288; clean-confirmed 8; clean-broad 3; mobile 0. Native counts are 440, 429, 12, 3 and 1. All other DNA negative counts are unchanged.

VirSorter2 explicitly includes short hallmark-containing sequences in its identified-virus output; the code's description of them as merely diagnostic is incorrect. See the [upstream output specification](https://github.com/jiarong/VirSorter2#detailed-description-on-output-files). The 141 RNA-positive native calls excluded by the adapter are **not** 141 additional viSUM false negatives: many are already recovered by other tools. A saved-evidence replay is needed to measure the policy impact.

Native parent counts were also reconciled for geNomad and Cenote-Taker3, with no discrepancies. Deep6 and DeepMicroClass2 rows describe the configured score/class screens, not a threshold-free model argmax or an exhaustively audited author-default execution. Avoid labelling the whole table an author-default leaderboard.

## What viSUM adds—and loses

Relative to geNomad, on DNA positives: 437 shared detections, 40 viSUM-only, 15 geNomad-only. On RNA positives: 633 shared, 155 viSUM-only, 13 geNomad-only. Thus the gains are real, but viSUM does not preserve every geNomad positive.

Across DNA+RNA positives, viSUM retains 1,265/1,413 (89.53%) versus geNomad 1,098/1,413 (77.71%): +167 net inputs, or +11.82 percentage points. This pooled result depends on this panel mixture; retain separate DNA/RNA results as the headline.

The gain is not uniform by genome category. For example:

| Truth genome category | viSUM | geNomad |
| --- | --- | --- |
| dsDNA | 216/225 | 217/225 |
| dsDNA-RT | 29/30 | 23/30 |
| ssDNA(+) | 46/71 | 52/71 |
| ssDNA(-) | 27/45 | 20/45 |
| dsRNA | 177/210 | 121/210 |
| ssRNA(+) | 221/221 | 212/221 |
| ssRNA(-) | 189/202 | 147/202 |
| ssRNA-RT | 15/33 | 16/33 |

Categories are preserved verbatim from benchmark truth, including mixed/unspecified categories in the supplementary table. Segment-level detection is not equivalent to detecting independent virus species. Review reverse-transcribing and ssDNA(+) losses before broad claims about coverage.

## False positives versus putative proviral regions

**The supplied evidence cannot definitively divide these calls into biological true positives and true false positives.** A pipeline-generated boundary is not an independent truth label. The defensible split is intact retention versus extracted-region candidates, with an additional evidence tier.

| Negative panel | Retained parents | Parents with extracted region(s) | Intact-only parents | Extracted-region parents with a CT3 virion hallmark somewhere in that parent's CT3 calls |
| --- | --- | --- | --- | --- |
| Cellular DNA | 15 | 6 | 9 | 6 |
| Mitochondrial DNA | 11 | 11 | 0 | 0 |
| Plasmid DNA | 67 | 41 | 26 | 17 |
| Plastid DNA | 1 | 1 | 0 | 0 |
| Retroelement DNA | 7 | 0 | 7 | 0 |
| Clean RNA, pooled | 27 | 0 | 27 | 0 |
| Total | 128 | 59 | 69 | 23 |

The 59 DNA parents yield 62 extracted records. The hallmark audit is **parent-level**, not proof that each hallmark overlaps the final retained boundary. Coordinate-specific adjudication remains necessary.

1. **Higher-priority proviral candidates: 23 parents.** Six cellular and 17 plasmid parents have both an extracted region and CT3 virion-hallmark annotations. Examples include `CELL_NEG_000003` (major capsid protein, *Dehalogenimonas* chromosome), `CELL_NEG_000191` (portal protein, *Edaphobacter* chromosome), and `CELL_NEG_000486` (minor capsid protein, *Campylobacter* chromosome). These warrant contextual investigation; a single short structural-profile match is still insufficient confirmation.
2. **Extracted but weaker/ambiguous origin evidence: 36 parents.** These lack a CT3 virion-hallmark annotation in the inspected parent summaries. All 11 mitochondrial candidates are in this category. Their annotations are predominantly DNA polymerase B, with one RNA-dependent RNA polymerase case; the plastid case is DNA helicase. Repeated replication-protein annotations are not equivalent to distinct virion functions. These could reflect shared/mobile-element homology or other biological context; do not present them as rescued true proviruses.
3. **Intact-only labelled-negative retentions: 69 parents.** Forty-two DNA and 27 RNA. These remain false positives under the existing truth labels, without an extracted-region explanation in these outputs. This does not prove that every sequence is biologically nonviral.

If one merely removes the 59 region-candidate parents from the DNA numerator, the residual is 42/2,390 (1.76%). That is an **intact-only retention rate, not a corrected FPR**. It is still higher than geNomad's label-based 1.59%, and the two quantities are not equivalent. Excluding those inputs entirely would instead give 42/2,331 (1.80%); neither operation is justified as biological adjudication. A true relabelling would require rescoring every tool on the same revised truth set.

### How to adjudicate without circularity

Start with the 23 structural-hallmark candidates and all 12 organelle parents. Recover original genomic sequence plus flanks; map the retained boundary and individual hallmarks into source coordinates; inspect gene order and coherent packaging/capsid/replication architecture; compare cellular, organelle and mobile-element alternatives. Record evidence from resources not used to set the current decision wherever possible, with a blinded second reviewer. Use `confirmed_nonviral`, `supported_viral_or_proviral`, and `unresolved` labels, retaining original labels and the adjudication rationale. Expression from an integrated locus requires expression/context evidence; an annotated transcript alone does not demonstrate an active virus or provirus.

Publish unchanged-label statistics first. Report adjudicated statistics separately, with the denominator and all tools updated consistently. Do not use viSUM's own confidence field as the relabelling criterion.

## Are Deep6 and viCAT signals weak?

Not necessarily. Among retained clean-confirmed RNA negatives supported by Deep6, the 16 scores range 0.711–0.9995, median 0.874. Among the six supported clean-broad retentions, median is 0.946. High model scores occur on labelled nonviral input; they are not established posterior probabilities of viral origin in this use case.

viCAT supports 18 of the 20 retained clean-confirmed negatives. Their best-hit bit scores range 139–4,715, median 586.5; all 18 have **one viral-supported locus**. Five of seven clean-broad retentions have viCAT support, median best bit score 836, with one to two supported loci. These are not all marginal alignments, but a strong homologous protein match is not evidence of an entire virus. Retained source annotations include transporters, polymerase/primase, MutS, helicase, collagen and dynein. Shared domains, reference contamination and genuinely virus-derived genes are alternatives requiring case-level review.

This is a locus-level assessment, not a claim that every raw DIAMOND match is independent support. Multiple reference hits to one ORF do not constitute multiple viral genes.

## Optional RNA floor: verified effect

The rule applies only to RNA, without strong support, where the qualified discovery pair is exactly Deep6 + viCAT. It requires a viral-supported locus with protein length ≥100 aa, bit score ≥100, and query and subject coverage ≥50%. Failure moves the candidate out of primary retention into provisional output; it does not establish a nonviral truth label.

| Endpoint | OFF | ON |
| --- | --- | --- |
| RNA viral primary parents | 789/874 (90.27%) | 788/874 (90.16%) |
| Clean-confirmed RNA | 39/1,200 | 20/1,200 |
| Clean-broad RNA | 11/400 | 7/400 |
| Clean RNA pooled | 50/1,600 (3.13%) | 27/1,600 (1.69%) |
| RNA mobile challenge | 6/200 | 0/200 |

The six DNA panels have byte-identical ON/OFF `final_metadata`, `sequence_disposition`, `database_candidates`, and `provisional_metadata` tables. Clean-RNA false-positive retentions fall 46% relative, with one additional RNA-positive loss. The RNA-negative OFF counts in the core evaluator are reconstructed from the ON manifest's demoted IDs; they are not misrepresented as a separately packaged OFF negative run. The positive OFF run is packaged directly.

Recommendation: retain this as a documented opt-in profile for now. Its observed tradeoff is useful, but thresholds were developed using these data. Freeze the rule before testing untouched data; do not keep tuning this panel and calling the results independent validation.

## Required changes and release gates

### Before a usable first public pre-release

1. **Resolve the VirSorter2 adapter contract.** Preserve hallmark-only native evidence with an explicit evidence type and missing numerical score; do not invent a score or silently equate it with strong ML support. Correct the misleading code comment. Test whether such evidence can supply qualified support or should remain review-only. Replay saved evidence through affected downstream stages and compare all ten panels before choosing a retention policy. This review fixes comparator accounting only; the production adapter is unchanged.
2. **Ship honest documentation.** Use the new README, distinguish defaults from the all-tools benchmark, describe primary/provisional/region candidates, document RNA rule OFF by default, and keep the limitations above adjacent to performance claims. The RNA controls include annotated RNA and prokaryotic CDS proxies; they are not an experimental mRNA cohort, and `clean-confirmed` is the panel name, not independent biological certification.
3. **Complete legal and citation metadata.** Obtain the lab's license decision; add the license, real author/maintainer metadata, third-party tool/database attributions and `CITATION.cff`. Public visibility alone is not an open-source license; see [GitHub's licensing guidance](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository). Do not redistribute MetaVR-derived or other third-party databases without checking their terms.
4. **Pass a fresh Linux installation smoke test.** From a clean clone, follow only the README; test a small DNA and RNA input, zero-call handling, CSV batch input, resume, and the chosen full configuration with already provisioned databases. Record Nextflow, Java, Conda/Mamba, tool revisions and DB checksums. Cached scientific runs do not prove a fresh installation works.
5. **Make failures and reproducibility explicit.** Validate duplicate batch prefixes before tasks start (they can collide in published output paths); provide automated tests in CI; expose unsupported/disabled tools distinctly from a successful zero-call run. Audit the repository/history for secrets and large private artifacts. Extend ignore rules deliberately; do not stage the entire current untracked workspace.
6. **Tag an explicitly preliminary release.** After the above gates, create a versioned pre-release and release notes listing limitations, known integration behavior, tested configurations and benchmark access. Use the actual tag/commit in the fellowship reference; do not invent a DOI or describe this as peer-reviewed validation.

### Before a manuscript-level performance claim

Add independent, reference-aware held-out data; length/taxon/parent-stratified evaluation and uncertainty; native-output audits for every comparator; blinded negative adjudication; boundary accuracy on truth-annotated integrated regions; strict/exploratory taxonomy coverage and correctness; sensitivity to tool ablation; and uncached resource measurements under comparable hardware. The current analysis establishes detection/retention performance, not taxonomic correctness, proviral boundary accuracy or computational superiority.

**Next execution should not be another full discovery rerun.** First repair/test the narrowly identified adapter behavior, replay the saved evidence and affected downstream steps with the existing cache, and regenerate this report. Separately perform one clean-install smoke test. Only after freezing that release candidate should the next expensive run be an untouched validation panel.

## Audit artifacts and reproducibility

The review extracted 3,032 structured/code/log members from the archive; large FASTAs/raw data were not all materialized. Eighty available manifest-listed outputs passed SHA-256 checks; 70 listed outputs were not extracted or not packaged and were **not hash-verified** by this check. The full archive was not represented as exhaustively validated.

Core evaluation (Python standard library only):

```bash
python bin/evaluate_benchmark_review.py \
  --review-root /path/to/extracted/consolidated_review \
  --outdir /path/to/new/review_tables
```

The evaluator expects the packaged `rna_pair_floor_*` and `rna_pair_OFF_*` layout and all applicable discovery outputs; it fails on missing applicable tables rather than scoring them as zero. Original-input IDs are reconciled against the saved disposition, and the supplemental truth join covers the seven original benchmark panels. The RNA source-annotation supplement uses the previously supplied RNA-negative truth metadata.

Machine-readable supplement: [metrics](data/metrics.tsv), [paired tool overlaps](data/paired_comparison.tsv), [negative-parent categories](data/negative_parent_audit.tsv), [source/hallmark audit](data/negative_source_hallmark_audit.tsv), [positive subgroups](data/positive_subgroups.tsv), [native concordance](data/native_summary_concordance.tsv), [RNA score summaries](data/retained_rna_score_summary.tsv), and [DNA regression checks](data/regression.tsv).

These results were generated from saved evidence; no new Nextflow integration run was performed in this review environment. The spreadsheet-analysis workflow's reconciliation rules were used to preserve source labels, normalize denominators, and distinguish missing/unsupported data from zero calls.
