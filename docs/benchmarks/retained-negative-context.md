# Interpreting retained nominal-negative sequences

This follow-up summarizes coordinate-local review of the saved September 2026
benchmarks. It supplements the [original all-tool comparison](preliminary-2026-09.md)
without changing its input labels or reported detection counts.

## What the review found

The benchmarked viSUM configuration retained 101 of 2,390 nominal-negative DNA
inputs and 27 of 1,600 clean-RNA negative proxies. Each input is counted once,
even if more than one region was retained. Review examined the retained outputs
using saved gene annotations, viral/nonviral protein comparisons and candidate
region boundaries, rather than treating a sequence-level tool label as the
entire explanation.

The retained set includes several different situations:

- **Coherent virus-like content.** Seven plasmid inputs contained capsid,
  portal and terminase modules supporting phage-like candidates:
  `PLASMID_NEG_000059`, `000102`, `000113`, `000306`, `000326`, `000346` and
  `000489` (all IDs use the `PLASMID_NEG_` prefix).
- **Viral homologs with unresolved origin.** Some sequences match viral
  proteins but also have plausible cellular or mobile-element explanations.
  Replication-protein homology alone is less specific than a coherent module
  with multiple virion functions.
- **Mixed or potentially overbroad outputs.** A retained interval can contain
  viral evidence alongside substantial cellular or organellar content. This
  raises a boundary/output-purity question without erasing the viral evidence.
- **Weak or unresolved calls.** Some retentions lack comparably specific
  support. These remain candidates for review, not confirmed biological errors
  or confirmed viruses.

Phage–plasmids are established biological entities, combining phage and plasmid
features ([Pfeifer et al., 2021](https://doi.org/10.1093/nar/gkab064)). A plasmid
source label therefore does not exclude viral biology. The seven benchmark
examples above have supporting computational evidence; this review does not
establish infectivity, a complete lifecycle or experimental confirmation.

## How to report performance

Use **nominal-negative retention rate** for the fraction retained from the
original negative panels. It is also the apparent, label-based FPR:

- DNA: 101/2,390 = **4.23%**.
- Clean RNA: 27/1,600 = **1.69%**.

These are not independently established biological false-positive rates. Do
not automatically subtract every proviral call, viral homolog or plasmid
candidate from the numerator. Equally, do not label every retained nominal
negative a proven nonviral error. All tools retain the same original truth
labels for comparison; a future adjudicated truth set would require rescoring
every tool consistently.

Detection of viral content, correctness of the retained boundaries, taxonomic
accuracy and evidence confidence are distinct endpoints. The current headline
benchmarks assess input-level detection, not all four at once.

### Conditional exclusion of supported virus-like candidates

The README additionally reports two explicit sensitivity analyses of the negative
labels. The localized review found seven coherent capsid/portal/terminase modules
and eleven additional plasmid inputs with multiple CT3 virion-hallmark genes
inside a retained interval. The latter can represent partial viral modules;
multiple hallmark genes need not imply distinct functions or a complete virus.
These observations support a virus-like interpretation without establishing
whether each source is a phage-plasmid, a plasmid-associated prophage or a
mislabelled assembly.

| Scenario | Excluded retained inputs | Remaining DNA retentions | Remaining DNA negative inputs | Conditional rate |
| --- | ---: | ---: | ---: | ---: |
| Original labels | 0 | 101 | 2,390 | 4.23% |
| Exclude coherent-module candidates | 7 | 94 | 2,383 | 3.94% |
| Exclude all multiple-localized-virion-hallmark candidates | 18 | 83 | 2,372 | 3.50% |

The excluded source IDs are:

- **Coherent-module set (7):** `PLASMID_NEG_000059`, `PLASMID_NEG_000102`,
  `PLASMID_NEG_000113`, `PLASMID_NEG_000306`, `PLASMID_NEG_000326`,
  `PLASMID_NEG_000346`, `PLASMID_NEG_000489`.
- **Additional multiple-hallmark set (11):** `PLASMID_NEG_000017`,
  `PLASMID_NEG_000137`, `PLASMID_NEG_000158`, `PLASMID_NEG_000203`,
  `PLASMID_NEG_000271`, `PLASMID_NEG_000283`, `PLASMID_NEG_000359`,
  `PLASMID_NEG_000416`, `PLASMID_NEG_000424`, `PLASMID_NEG_000461`,
  `PLASMID_NEG_000486`.

The 18-input scenario includes the seven-input set. Each original input is counted
once and removed from both numerator and denominator: `(101 - k)/(2390 - k)`.
Neither scenario excludes RNA inputs, so clean-RNA retention stays at
27/1,600 (1.69%). No positive-panel counts, primary-retention rules or saved
outputs change. A boundary-only provirus call or an isolated viral homolog is
not sufficient for either exclusion criterion.

These are **conditional, evidence-adjusted FPR scenarios**, not an independently
adjudicated biological FPR. The evidence review was conducted on retained
nominal negatives, not a blinded reassessment of the entire negative panel.
They show how the interpretation changes if these explicitly identified
virus-like candidates are outside the intended nonviral controls. Retention of
such candidates can be a useful discovery outcome. Remaining retentions are
not thereby proven biological errors, and excluded candidates are not thereby
proven infectious viruses or accurately bounded proviruses.

## Experimental policy comparison

A saved-evidence replay tested localized mixed-content review, CT3
replication-only review, and their combination. Discovery outputs and selected
boundaries were held fixed; demoted records remained available as provisional
candidates rather than being declared nonviral.

| Policy | DNA positives /539 | RNA positives /874 | DNA nominal negatives /2,390 | Clean RNA negatives /1,600 |
| --- | ---: | ---: | ---: | ---: |
| Benchmarked policy | 477 | 788 | 101 | 27 |
| Mixed-content review | 477 | 788 | 88 | 27 |
| CT3 replication-only review | 476 | 788 | 83 | 24 |
| Both experimental rules | 476 | 788 | 72 | 24 |

All seven coherent-module plasmid candidates remained retained in every arm.
The replication-only option demoted one labelled DNA virus,
*Nitrosopumilus spindle-shaped virus 1*. These experiments did not change the
default retention policy. The experimental replay implementation remains local
development work and is not included in this documentation/default-settings update.

For this comparison, mixed-content review required at least three host-supported
loci exceeding viral-supported loci within the same source (class-aware viCAT
competition or exact-interval CheckV counts). A separate arm flagged at least
three CT3-annotated respiratory/photosynthetic genes exceeding distinct
virion/RdRp genes. Host evidence outside the output interval was not used as a
veto. The replication-only option withheld CT3 primary-origin votes and strong
boundary rescue when DNA-replication hallmarks were present without virion or
RdRp hallmarks. These are exploratory thresholds, not independently validated
operating points.

The all-tool benchmark configuration enabled viCAT and the optional RNA
homology floor. The published default configuration still requires those
options to be selected explicitly; do not describe the table as a benchmark
of the unmodified quickstart. These panels informed development and are not an
untouched external validation set.
