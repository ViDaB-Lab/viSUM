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

The README reports a broader virus-like-content discovery endpoint, rather than
restricting exclusions to the strongest plasmid modules. A retained DNA input
qualifies for the regional/hallmark set if it has either a localized CT3 virion
or RdRP hallmark, or an overlapping CheckV viral region with positive viral-gene
counts and local viCAT viral-supported protein loci. This is a retrospective
interpretation of saved evidence, not a new primary-retention rule.

| DNA source | Retained inputs | Inputs yielding extracted candidate regions | Regional/hallmark-supported inputs |
| --- | ---: | ---: | ---: |
| Cellular | 15 | 6 | 5 |
| Mitochondrial | 11 | 11 | 5 |
| Plasmid | 67 | 41 | 38 |
| Plastid | 1 | 1 | 0 |
| Retroelement | 7 | 0 | 0; seven endogenous-retroviral controls considered separately |
| Total | 101 | 59 | 48 |

The 59 inputs yielded 62 extracted records. The 48 supported inputs comprise
34 with localized virion/RdRP hallmarks and 14 additional inputs with CheckV
regional evidence plus viral homology. Of those 48, 36 yielded extracted regions
and 12 were retained intact. The 38 plasmid inputs include 28 with virion hallmarks
and 10 additional regional/homology-supported candidates. Eighteen of the 28 have
multiple virion-hallmark genes, including the seven coherent modules above.

All seven retained retroelement controls are source-labelled ERV1/ERVK-family
repeats. Their saved TEsorter fields also show `direct_hmm`, order `LTR`,
superfamily `Retrovirus`, evidence strength `qualified`, and status
`viral_like_mobile_element_with_viral_support`. Thus retroelement identification
was actually present in the pipeline outputs, not inferred only from source labels.
Counting these as intended discoveries is appropriate for an endpoint including
endogenous viral/mobile-element content. It is not evidence of seven infectious
retroviruses. Their saved `viral_entity_interpretation` remains `viral_contig`;
TEsorter's context must be read alongside that field, not mistaken for a claim
that viHARMONY assigned a dedicated endogenous-element entity label.

| Scenario | Excluded retained inputs | Remaining DNA retentions | Remaining DNA negative inputs | Conditional rate |
| --- | ---: | ---: | ---: | ---: |
| Original labels | 0 | 101 | 2,390 | 4.23% |
| Exclude regional/hallmark-supported candidates | 48 | 53 | 2,342 | 2.26% |
| Also exclude endogenous-retroviral controls identified by TEsorter | 55 | 46 | 2,335 | 1.97% |

The excluded source IDs are:

- **Cellular (5):** `CELL_NEG_000003`, `CELL_NEG_000077`, `CELL_NEG_000191`,
  `CELL_NEG_000202`, `CELL_NEG_000486`.
- **Mitochondrial (5):** `MITO_NEG_000141`, `MITO_NEG_000177`, `MITO_NEG_000272`,
  `MITO_NEG_000300`, `MITO_NEG_000489`.
- **Plasmid (38):** `PLASMID_NEG_000017`, `PLASMID_NEG_000059`, `PLASMID_NEG_000081`,
  `PLASMID_NEG_000088`, `PLASMID_NEG_000089`, `PLASMID_NEG_000102`, `PLASMID_NEG_000113`,
  `PLASMID_NEG_000136`, `PLASMID_NEG_000137`, `PLASMID_NEG_000142`, `PLASMID_NEG_000158`,
  `PLASMID_NEG_000189`, `PLASMID_NEG_000203`, `PLASMID_NEG_000251`, `PLASMID_NEG_000261`,
  `PLASMID_NEG_000271`, `PLASMID_NEG_000274`, `PLASMID_NEG_000279`, `PLASMID_NEG_000283`,
  `PLASMID_NEG_000286`, `PLASMID_NEG_000288`, `PLASMID_NEG_000306`, `PLASMID_NEG_000316`,
  `PLASMID_NEG_000326`, `PLASMID_NEG_000346`, `PLASMID_NEG_000353`, `PLASMID_NEG_000359`,
  `PLASMID_NEG_000397`, `PLASMID_NEG_000405`, `PLASMID_NEG_000416`, `PLASMID_NEG_000424`,
  `PLASMID_NEG_000432`, `PLASMID_NEG_000436`, `PLASMID_NEG_000461`, `PLASMID_NEG_000481`,
  `PLASMID_NEG_000486`, `PLASMID_NEG_000489`, `PLASMID_NEG_000499`.
- **Endogenous-retroviral controls (7):** `NEG_TE_LTR_0005`, `NEG_TE_LTR_0006`,
  `NEG_TE_LTR_0007`, `NEG_TE_LTR_0014`, `NEG_TE_LTR_0032`, `NEG_TE_LTR_0038`,
  `NEG_TE_LTR_0117`.

The 55-input scenario includes the 48-input set. Each original input is counted
once and removed from both numerator and denominator: `(101 - k)/(2390 - k)`.
Neither scenario excludes RNA inputs, so clean-RNA retention stays at
27/1,600 (1.69%). No positive-panel counts, primary-retention rules or saved
outputs change. An extracted-region label alone or an isolated viral homolog
without the specified hallmark/regional support is not sufficient for the
48-input set. RNA homologs are not automatically excluded.

These are **conditional, evidence-adjusted FPR scenarios**, not an independently
adjudicated biological FPR. The evidence review was conducted on retained
nominal negatives, not a blinded reassessment of the entire negative panel.
They show how the interpretation changes if these explicitly identified
virus-like candidates are outside the intended nonviral controls. Retention of
such candidates can be a useful discovery outcome. Remaining retentions are
not thereby proven biological errors, and excluded candidates are not thereby
proven infectious viruses or accurately bounded proviruses. In particular,
`PLASMID_NEG_000436` has regional/homology evidence but its retained interval
also contains substantial host content. Its exclusion is a viral-content
interpretation, not endorsement of the entire retained sequence as a clean
viral genome. These calculations assess discovery of content, not output purity.

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
