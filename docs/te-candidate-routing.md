# TE candidate routing

TEsorter contributes element identity, architecture and evidence-scope flags.
The flags distinguish strong TE support from partial, mixed or localized matches.
Read them alongside the accepted taxonomy and other origin evidence rather than
interpreting a single retrovirus-like match as a complete viral genome.

## Policy

- Strong, direct, complete and coherent Gypsy, Copia or Bel-Pao evidence routes
  to TE candidates when the accepted uncapped taxonomy is compatible or broader.
  The corresponding families are Metaviridae, Pseudoviridae and Belpaoviridae.
- Incompatible accepted taxonomy, or an explicit unresolved family/higher
  consensus conflict, prevents automatic TE-only routing. Such a record enters
  viral candidates only if it passes the ordinary viral-eligibility safeguards.
- Partial, mixed or ambiguous TE evidence annotates ordinary eligibility instead
  of acting as an early blanket veto. It supplies no new viral-support vote.
- Explicit retrovirus-like evidence remains compatible with viral candidacy.
  It does not establish endogenous versus exogenous status or infectivity.
- Nonselected individual taxonomy predictions remain audit warnings, not a
  substitute for accepted consensus. No taxonomy is invented from a TE label.

The actual TE interpretation is still used by the restricted RNA CheckV exception.
The RNA pair-locus floor, cellular/plasmid safeguards, their existing explicit
overrides and disputed-region boundary requirements still execute. Reporting
rank caps cannot hide an incompatible uncapped family from TE routing.

## Outputs

Under each sample's `viharmony/` folder, the relevant files are:

- `TE_candidates/<sample>.TE_candidates.fasta`
- `TE_candidates/<sample>.TE_candidates.summary.tsv`
- `<sample>.viral_candidates.fasta`
- `<sample>.viral_candidates.summary.tsv`
- `audit/<sample>.candidate_routing_audit.tsv.gz`

The two candidate FASTAs are disjoint. The combined viral FASTA contains retained
and provisional viral candidates with explicit labels. Legacy views in
`audit/legacy/` are not additional candidates. Review tables may overlap the
exports because warnings do not automatically remove a candidate. See the
[reporting bundle](reporting-bundle.md) for all columns and locations.

For TE-only records, `viral_decision` is `TE_candidate`. The detailed metadata and
always-written routing audit preserve `pre_TE_routing_decision` and
`ordinary_eligibility_decision`. The latter is the result after bypassing only the
early TE veto and retaining all later origin safeguards. TE routing never turns
an otherwise blocked conflict into a viral export. Old review annotations remain
available, including for records whose ordinary decision is nonviral.

Summaries include original/final/parent IDs, coordinates, sequence length, ordinary
eligibility, confidence, destination, TE flag/architecture/scope, TE order and
superfamily, viral-evidence scope, locus counts, cluster locations, supporting
tools and the full semicolon-delimited reported taxonomy. Missing values are
`NA`, distinct from zero. The audit also preserves uncapped taxonomy, domain
assignments, competing consensus and the original routing reason.

## Clarified flags

| Flag | Meaning |
| --- | --- |
| `TE_supported` | Strong, complete coherent TE architecture with compatible consensus. |
| `TE_viral_conflict` | Strong TE evidence competes with accepted or explicitly conflicting taxonomy. |
| `incomplete_TE_architecture` | TEsorter explicitly reports incomplete architecture; this does not imply weak TE matches. |
| `limited_TE_domain_support` | Domains are detected but complete architecture is not established. |
| `mixed_domain_assignments` | Different TE/retroviral assignments, not proof of a hybrid or contamination. |
| `ambiguous_mobile_element` | Ambiguous identity such as Maverick/Polinton; conflicting multiple domains are not required. |
| `localized_TE_domain_support` | Split-window TE evidence cannot automatically define the entire candidate. |
| `retrovirus_compatible` | Explicit viral-like domains; endogenous/exogenous status unresolved. |
| `weak_TE_signal` | No qualified direct TE-domain support. |

Architecture and window scope are separate fields, so identity flags do not hide
localized evidence. Different clades within one TE superfamily are not treated
as conflicting superfamilies. A missing TE hit is not evidence of a retrovirus.

## Viral-evidence scope

`viral_evidence_scope` distinguishes localized clusters alone, localized clusters
with record-level predictions, single-locus support, multi-locus support, record-
level predictions alone, and absence of qualified viral-origin support. A
record-level prediction is not proof that every base is viral or that a complete
genome exists. An audit's `applied_scope=region` alone does not establish that
the sequence was extracted; intact input contigs can receive that audit label.

Cluster locations are 1-based inclusive coordinates on the candidate record.
They are evidence locations, not newly inferred extraction boundaries, and need
not cover every viral-supported locus. For extracted children, they are not
silently relabeled as original-parent coordinates.

For example, a short viral-supported cluster inside a long scaffold can be
reported as localized evidence. It does not establish viral origin for the entire
scaffold. Likewise, an ambiguous mobile-element label is not evidence of a
confirmed hybrid. Inspect the candidate's evidence scope and supporting loci.
