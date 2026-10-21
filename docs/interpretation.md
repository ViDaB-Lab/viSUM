# Interpret viral, TE and taxonomy results

[Output files and columns](reporting-bundle.md) · [TE flags](te-candidate-routing.md)

## Start with the summary and matching sequences

`<sample>.viral_candidates.summary.tsv` describes the sequences exported in
`<sample>.viral_candidates.fasta`. The combined export includes retained and
provisional candidates. Check `candidate_status`, `viral_confidence`, conflicts
and `review_flags` rather than treating every FASTA record as equally supported.

- **Retained candidates** pass the applicable viral-retention rules. Evidence
  strength and conflicts still need interpretation.
- **Provisional candidates** remain accessible with a provisional label; they
  should not be presented as confirmed viral genomes.
- **TE candidates** are exported separately under `TE_candidates/`. Strong,
  coherent TE evidence can explain an otherwise virus-like prediction.
- **Review tables** collect reasons and warnings. A record may be flagged in a
  review table and still be exported as viral or TE. A review-only/excluded
  record has no additional review FASTA.

The confidence label describes support, not a probability of being viral.
No category establishes infectivity. An empty export does not necessarily mean
the input is biologically nonviral: tools have detection limits.

## Original sequences and extracted regions

Summary fields retain original IDs, normalized parent IDs, record type and
coordinates where applicable. Missing values are `NA`. A derived region is not
an additional independent original input. Check its parent and coordinate frame
before comparing counts or using the sequence downstream.

`viral_evidence_scope` separates localized homology from record-level predictions.
A predicted viral region or a cluster of viral genes does not prove that the
whole parent is viral. See [provirus boundaries](provirus-boundary-policy.md).

## Taxonomy and reporting depth

`reported_taxonomy` is the full semicolon-separated lineage. The default
`--harmonizer_max_rank family` masks deeper ranks with unclassified placeholders;
`<sample>.taxonomy.tsv` instead has separate rank columns and `NA` where unassigned.

You can permit a deeper supported call with `--harmonizer_max_rank genus` or
`species`. This only raises the reporting ceiling; it does not relax the
underlying evidence rules or manufacture an assignment. Possible ceilings are
`domain`, `realm`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`.

A call may stop above the configured ceiling when evidence is insufficient or
conflicting. A sequence's viral-origin evidence and its taxonomic resolution
are separate questions. Short or single-gene sequences need particular caution.

viSUM uses the shipped, pinned ICTV rank table and matching reconciliation bundle
to resolve supported taxonomy names. It does not automatically adopt a new ICTV
release during setup, analysis or resume. Original names, reconciled names and
unresolved mappings remain auditable. A successful name reconciliation is not
additional biological evidence for a sequence's classification.

The reporting audit preserves uncapped predictions even when the summary stops
at family. Deeper conflicts still matter; a reporting cap does not hide them
from the decision rules. `lineage_inferred` means a rank came from a resolved
lineage rather than a direct tool vote at that rank.

## CheckV quality and completeness

Read `quality/<sample>.checkv.tsv` alongside the summary. A value of
`Not-determined` means completeness was not established, not that a sequence is
nonviral. A direct `Complete` assessment gives `checkv_complete_call=true`;
determined high/medium/low categories give `false`, not proof of incompleteness.
Missing, undetermined or non-direct assessments give `NA`.

Check the assessment scope. Parent-only or proviral aggregate context does not
establish completeness of an extracted child. A complete call does not prove
infectivity or complete recovery of every segment of a segmented viral genome.
See [CheckV context](checkv-evidence-context.md).

## vCONTACT3 groups

`groups/<sample>.vcontact3.tsv` preserves the tool's predictions by searched
database and rank. Named predictions, novel groups, unplaced groups and
singletons are not interchangeable:

- A named tool prediction is evidence, not automatically the accepted consensus.
- A novel group is local gene-sharing context, not a formally named ICTV taxon.
- A singleton or missing match does not establish a new species or assigned realm.
- Database domain identifies which reference subset was searched, not a host call.

Raw predictions can be deeper than the final reporting cap. Consult taxonomy
audits to see whether a prediction was accepted, reconciled or rejected.

## TE evidence and flags

Strong, complete and coherent Gypsy/Copia/Bel-Pao evidence with compatible
consensus can route a sequence to the separate TE export. Partial or mixed TE
support does not by itself remove an otherwise eligible viral candidate; its
scope and flags remain visible. Retrovirus-compatible evidence alone cannot
distinguish endogenous from exogenous origin.

The viral and TE FASTAs are disjoint. Review views may overlap exports because
they are warnings, not a third set of candidate sequences. See [TE evidence](te-candidate-routing.md).

## Detailed evidence

`audit/` retains evaluated-record metadata, sequence maps, dispositions and
taxonomy/routing audits. `--harmonizer_audit full` adds row-level evidence/vote
audits; the summary is not a replacement for inspecting conflicting evidence.
Legacy views under `audit/legacy/` are compatibility outputs, not extra detections.
The manifest inventories the current bundle with checksums.
