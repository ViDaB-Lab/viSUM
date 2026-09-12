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
