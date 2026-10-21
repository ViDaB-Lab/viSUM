# viHARMONY reporting bundle

`reporting-bundle-v1` organizes the completed viHARMONY results. The summary,
quality, group and review views describe existing decisions; they do not
independently reclassify candidates. Start with the viral summary and matching
FASTA, then consult the detail tables below. See [interpreting results](interpretation.md).

## Layout

Under `<sample>_results/viharmony/`, the combined viral FASTA and summary plus
`<sample>.taxonomy.tsv` are at the top level. `TE_candidates/` contains the TE
FASTA and summary. `quality/` contains `<sample>.checkv.tsv`; `groups/` contains
`<sample>.vcontact3.tsv`; `review/` contains `<sample>.review.tsv` and
`by_reason/<sample>.review.<reason>.tsv` for reasons actually present.

`audit/` contains the unchanged detailed metadata, review queue, sequence map,
input dispositions, routing audit and enabled evidence/taxonomy/reconciliation/
reporting audits. `audit/legacy/` retains primary/provisional normalized and
original-ID FASTAs, provisional metadata, database-candidate FASTA/table and
the all-candidates FASTA. They are compatibility views, not additional candidates.
No review FASTA is written or declared as a workflow output.

All output filenames start with the sample prefix. Empty tables have headers;
empty viral/TE collections have empty FASTAs. Missing scalar values are `NA`,
distinct from measured zero and from explicit upstream `Not-determined`.

## Exact table schemas (in order)

Both candidate summaries (26 columns):

```text
sample_id final_sequence_id original_sequence_id normalized_parent_id
record_type coordinates length_nt candidate_status viral_confidence
strong_tools qualified_tools conflict_tools reported_taxonomy
reported_taxonomy_rank taxonomy_support_status taxonomy_supporting_tools
viral_evidence_scope TE_flag TE_architecture_status TE_evidence_scope
tesorter_orders tesorter_superfamilies checkv_quality checkv_complete_call
vcontact3_group_status review_flags
```

The original ID identifies the input source, including the source parent of
extracted records. Parent ID and coordinates are `NA` for an intact input;
the original self-parent mapping remains unchanged in the detailed audit.
Conflict tokens retain type (`cellular:tool;plasmid:tool`). Candidate status
copies `viral_decision`, not hypothetical eligibility before TE routing.
Viral confidence does not measure confidence in TE identity. Taxonomy support
copies the existing reported-rank confidence/tools, without inventing votes.
The full semicolon lineage keeps all nine ranks and unclassified placeholders.

Viral-only taxonomy table (10 columns):

```text
final_sequence_id domain realm kingdom phylum class order family genus species
```

Names omit prefixes and missing/capped ranks are `NA`. IDs exactly match the
combined viral FASTA. No abundance information or TE-only records are included.

CheckV table (24 columns):

```text
sample_id final_sequence_id checkv_sequence_id assessment_status
assessment_scope assessment_coordinates checkv_record_length_nt
checkv_quality miuvig_quality checkv_complete_call completeness_percent
completeness_method contamination_percent provirus proviral_length_nt
viral_genes host_genes complete_genome_prediction complete_genome_confidence
topology aai_confidence hmm_completeness_lower_percent
hmm_completeness_upper_percent warnings
```

There is a row per applicable saved CheckV evidence record; no matching record
gets an explicit placeholder row, not a fabricated completed/skipped run status.
The raw quality fields remain visible even on parent/aggregate context rows.
`checkv_complete_call=true` means a directly applicable assessment has the raw
quality category `Complete`. Determined high/medium/low categories give `false`
(not called complete, not proven incomplete). Missing/undetermined assessments
and non-direct context give `NA`. Summary quality also requires exact-candidate
scope and agreement if multiple direct assessments exist.

The current standardizer copies original CheckV aggregate values onto proviral
rows. An exact provirus ID is therefore labeled `proviral_aggregate_context`,
not treated as independent child completeness. Parent-only rows are labeled
`parent_context`; neither supplies a child's complete flag or fabricated
assessment coordinates. Complete calls do not establish viral identity,
infectivity, or complete representation of a segmented viral genome.

vCONTACT3 table (14 columns):

```text
sample_id final_sequence_id database_domain database_version assignment_status
rank raw_prediction prediction_type evidence_type network_support vog_support
standardization_decision host_domain vog_host_domain
```

Each candidate/database/rank prediction is a separate row, including subfamily.
Prediction types are `named_prediction`, `novel_group`, `unplaced_group` and
`unresolved`. Missing predictions remain explicit singleton/no-matching rows;
`singleton` is never presented as an assigned realm. Raw prediction spelling,
reference/VOG/novel evidence and zero-valued supports are preserved, not re-scored.
Database domain identifies the searched database, not a host assignment.
Novel groups are local analysis context, not globally stable ICTV taxa.
Named prediction alone does not mean a valid/accepted formal taxonomy vote;
the unchanged vote/reconciliation audits explain actual taxonomy selection.
The reporting cap applies to reported taxonomy, not the tool-prediction audit.

Review tables contain the 26 summary columns followed by `review_category`
and `review_reason`. Each original candidate/reason pair appears once in the
master table and in its reason-specific view. Categories distinguish flagged
viral candidates, flagged TE candidates, excluded candidate audits and
review-only records. View membership is non-exclusive and never changes exports.
The detailed original review queue and original reasons remain unchanged.

## Manifest and publication

`audit/<sample>.harmonizer_manifest.json` inventories every current output
except itself. Its `reporting_bundle.output_root` is `..`: output file paths
are relative to the enclosing `viharmony/` bundle root, not to `audit/`.
Checksums cover the actual relocated/generated files. Existing scientific
policy/reference metadata is unchanged. Native outputs retain relative folders
and publish with overwrite/failure checking enabled.
The harmonizer source fingerprint participates in Nextflow's cache identity.
Changing a reporting helper can therefore rerun viHARMONY without implying that
upstream search results must be rerun. Keep the original work/cache directories
and let Nextflow decide reuse; do not edit cached outputs.

Use a NEW output directory when migrating from an older layout. The publisher
does not remove old flat files or now-obsolete reason views from a previous
publication. The manifest is the authoritative current output inventory; this
update never deletes older published results. For a same-version run, a fresh
outdir can be combined with the appropriate native cache/workdir and `-resume`.
Do not assume cross-version beta.1 cache compatibility when upgrading.
