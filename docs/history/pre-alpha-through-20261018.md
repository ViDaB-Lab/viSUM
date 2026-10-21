# Historical development record; superseded by ../updates-prior-to-alpha-release.md

# Updates to evaluate before the alpha release

This is the release checklist for the isolated alpha-development checkout, including deferred decision/reporting proposals. Engineering work may proceed in isolation while student analyses stay on their frozen build. It is not approval to change scientific policy, migrate active runs, silently change their calls, or mix results from different revisions. This isolated copy supersedes the older 13 October planning status; historical evidence below is retained.

## Reporting bundle integration — 18 October 2026

The approved reporting layout and 26-column viral/TE summaries are integrated
only in this isolated checkout. It adds viral-only rank-separated taxonomy,
explicit-scope CheckV quality/completeness, per-database/per-rank vCONTACT3
predictions, and table-only review views split by original reason. Scientific
artifacts retain their values under audit/ and audit/legacy/; no review FASTA
is produced. Manifest paths now resolve from the enclosing bundle root.
See `reporting-bundle.md` for exact schemas and migration/publication limits.
Saved-input reporting-only replay passed all ten panels/3,068 candidates:
1,745 viral candidates and 40 TE candidates, identical to validated TE routing.
All 50 pre-existing non-run function/class ASTs and the complete TE helper AST
are unchanged. Original metadata, candidate decisions, FASTA content, mapping,
dispositions and full-mode audits match the previous validated package. All
218 saved package files and 1,252 protected primary-checkout/saved-input files
remain unchanged. One/eight-worker and full/compact/none reporting checks pass
on DNA retroelement and RNA viral panels; compact reconciliation audits retain
exactly the changed/unresolved subset of the full audit, as before.
The first reporting Linux review passed all three saved replays (1,211
candidates), exact nested publication, and unchanged six-task cached resume,
including repair of a stale summary and missing nested CheckV table. The
helper-only update incorrectly stayed cached: the SHA-256 was computed in a
script-local variable, not a dependency tracked by native task hashing. This
is corrected in isolation by referencing a fingerprint-valued `task.ext`
property in the script. The corrected native cache invalidation gate is still
pending; scientific Python sources, output schemas and all decision rules are
unchanged. No student source, database, cache or published result was changed.

The uploaded TE-routing Linux review has meanwhile closed that component's
native gate: 65 tests, all 3,068 saved candidates and 104,312 expected fields
passed; three native panels (1,211 candidates) and unchanged six-task cached
resume passed, including real published-file checks. TE routing can be frozen
as validated; this reporting-only update must preserve that scientific behavior.
Resource v2, ICTV main-entry-point and rank-cap publication Linux reports also
passed, superseding their historical pending labels below. These component
passes do not establish a fresh installation or final assembled alpha benchmark.

## TE routing update — 18 October 2026

The user-approved TE routing and clarified evidence-scope labels are integrated
only into this isolated checkout. Saved-input replay covers 3,068 candidates
across ten panels and matches the standalone prototype: 40 TE-only candidates,
97 additional flagged viral candidates, and all 1,648 previous viral exports
preserved. Taxonomy, original evidence and sequence checksums are unchanged.
One/four-worker and full/compact comparisons pass on the saved DNA retroelement
and RNA viral panels. The 53 TE-specific regression tests pass locally.

New sample-prefixed TE and combined viral FASTAs/summaries and an always-written
routing audit make destination and evidence scope explicit. Original decisions
and ordinary eligibility remain auditable; summaries/audits stream in the ordered
parent to avoid another dataset-sized routing-row collection. See
`te-candidate-routing.md` for policy, flags, output fields and limits.

This TE portion of reporting is approved and implemented. The broader folder
redesign/slim-summary/CheckV presentation work is not declared complete by this
update. Native Linux/Nextflow entry-point, publication and cache checks for this
exact candidate remain required before deployment. No student/server analysis,
database, search result or existing published output was changed.

## Previous release status — 17 October 2026

Scope: `C:/Users/alexv/.codex/worktrees/viharmony-performance/viSUM`, based on `d4b90f3c235eae41b03a24ea4fbfd62b5febdb1c` plus uncommitted candidate changes. No student/server checkout is updated by this checklist or resource-default change.

| Area | Status and evidence |
| --- | --- |
| viHARMONY indexing and bounded multiprocessing | Implemented in isolation. Saved-panel equivalence and completed `Tt_4_S2` one/eight-worker replay reproduced the frozen outputs; see historical measurements below. Default workers remain one, capped by `max_cpus`. |
| Native launch and resume | Implemented. Shipped `nextflow.config` includes `visum.config`; the optional wrapper only forwards to Nextflow. No custom resume/cache manager or preparation-workdir recovery flags. Keep both native task cache and work outputs. |
| Separate database setup | Implemented. `--setup` prepares persistent resources separately; analysis validates and consumes prepared paths. Setup validation may run again without rebuilding valid assets. Sample-specific viCAT ORF/reference-hit processing stays in analysis. |
| Cache regression tests | Linux fixture reports exercised completed/interrupted resume, worker changes, missing-output detection, and inline-preparation negative controls. Actual-entry-point tests exercised separate repeated setup and cached analysis. These synthetic tests do not prove real installation or biological accuracy. |
| Real-tool DNA/RNA smoke test | Completed using the frozen v3 test package with scoped runtime repairs. Final resume: 41 cached and 17 completed tasks, no failed/aborted tasks. All 41 cached hashes match prior successful tasks; all eight geNomad, VirSorter2 and viCAT searches were cached. Both viHARMONY manifests passed all 26 listed hashes; FASTA/metadata IDs reconcile and both samples have 1,000 unique input dispositions. |
| Defaults/readiness simplification | Implemented in the candidate: installation-relative paths resolved after overrides, shared setup/preflight structural contracts, incomplete-install rejection, duplicate-prefix checks, and native local CPU/memory budgets. Latest changes still need target-Linux validation. |
| Runtime repairs | Candidate includes task-local short TEsorter temporary roots and explicit vCONTACT3 PyArrow dependency plus in-memory Parquet check. Server smoke used targeted overrides/private-environment repair, not a fresh installation of these candidate files. |
| viCAT nonviral provenance | Existing 480-assembly inventory, build-matching records and hashes recovered under `assets/vicat/nonviral_classaware_v1_20260904/`. Unmatched 404-exclusion evidence remains separate. No database rebuild or redistribution clearance is claimed. |
| Resource defaults — this update | All 20 `*_time` defaults are `null`; removed the SLURM profile's 12-hour limit and both viCAT validation-only time limits. Existing per-tool/site overrides remain supported. viHARMONY requests 32 GB instead of 4 GB; overall default remains 64 GB and default workers remain one. |
| Reporting and scientific policy | Family-default `harmonizer_max_rank` implemented in isolation on 17 October: reporting-only mask after inference, retention and review. Uncapped consensus preserved in an always-written sample-prefixed reporting audit. Species-support and viral-retention rules unchanged. Option B summary/layout remain pending; full semicolon-separated taxonomy, original/parent labels, `NA` missing cells, sample-prefixed exports and no review FASTA are agreed design requirements. |
| ICTV nomenclature reconciliation | Standalone pinned MSL41.v1 resolver passed Linux validation: 74 dependency-free tests, 445 reproduced flagged baselines, exact detailed local/Linux agreement, and 160 unchanged saved inputs. Integrated only in isolated viHARMONY with original evidence preserved, a streaming reconciliation audit, and bundle/manifest staged as native task inputs. Final local regression/contract suite: 195 passed plus 113 subtests. Complete-harmonizer replay passed all 10 panels/3,068 candidates: identical viral decisions, primary/provisional FASTAs, original evidence audits, mapping/disposition, and all non-taxonomy metadata; all strict taxonomy outputs matched independent standalone expectations. Historical RNA settings and current defaults were both checked, with no retention differences in these saved panels. All 337 source-input hashes unchanged. Report: primary repository `outputs/ictv_reconciliation_integration_20261016/full_replay_final/report.json`. Native target-Linux integration/cache verification and deployment remain separate gates. No species-rule, retention-rule or reporting-rank-cap change. |
| Deployment/public release | Not done. Candidate changes remain isolated/uncommitted; clean installation, final validation, release provenance and packaging remain gates. Student rescue publication/obsolete workflow shutdown must be confirmed separately rather than inferred from an old note. |

### Reporting cap and agreed output design — 17 October 2026

Implemented only in the isolated checkout: `harmonizer_max_rank = 'family'`,
with validated rank choices, native module CLI forwarding, and a reporting-only
mask after unchanged taxonomy inference, viral adjudication and review-reason
generation. Both strict/exploratory final lineages are capped. Species support
rules are unchanged; `species` restores the previous reporting depth.
`<sample>.taxonomy_reporting_audit.tsv.gz` always preserves uncapped consensus
fields, including deeper conflicts and species-support explanations, even with
compact/none audits. Original evidence and full-mode taxonomy votes stay intact.
Reported confidence/support is attached to the displayed rank; an ICTV-path-only
ancestor is explicitly `lineage_inferred`, not a fabricated direct family vote.

Agreed but not yet implemented: Option B slim summary using `reported_taxonomy`
for the complete semicolon-separated lineage, `NA` for missing scalar values,
original/parent labels and coordinates, sample-prefixed exports, table-only
review views (no review FASTA), and separate quality and vCONTACT3 group tables.
Combined primary/provisional exports and the folder reorganization remain part
of that later output change. TEsorter interpretation/routing remains undecided.

This is not deployment. No student run, production source checkout, database,
upstream search, existing server cache or published student result is changed.
Native Linux execution/cache checks for the new reporting cap remain a gate.

Configured local regression suite: **549 passed, 250 subtests passed, 15 skipped,
one deselected**. The deselected check requires an extracted student-recovery
package and is unrelated to this reporting change; it is not counted as passed.
The first broad attempt lacked archived-fixture paths and hit the existing
Windows pytest-temporary-directory ACL. The configured rerun used a new private
test directory and the verified local archived package; no application logic
was changed to bypass those environment errors. The historical indexing guard
still compares every biological statement, allowing only the exact two final
reporting statements; frozen-output comparisons explicitly use species depth
and validate the newly added audit checksum before excluding that extra file.

Independent complete-harmonizer replay against the pre-cap source snapshot:
**10 saved panels / 3,068 candidates passed**. Family and genus lineages matched
an independent rank mask; all viral/non-taxonomy metadata, primary/provisional
and other existing FASTAs, mapping/disposition, review membership/reasons,
original evidence, taxonomy votes and reconciliation audit contents remained
identical. Every uncapped consensus field was recovered exactly from the new
reporting audit. `species` reproduced all pre-cap scientific outputs. Selected
DNA/RNA panels also matched serial/eight-worker outputs. All 334 saved-file
fingerprints and 41 pre-existing policy/helper function ASTs were unchanged.
Family masked deeper consensus in 1,495 candidates; genus masked species in
1,096 candidates (these counts include strict and/or exploratory changes, not
only primary retained sequences). Report: primary repository
`outputs/harmonizer_rank_cap_20261017/saved_replay/report.json`.
No search, database preparation, Nextflow scientific workflow or student input
was launched for this replay. This local result does not claim native Linux
execution/cache validation of the new reporting-cap setting.

### Reporting-cap native package — 18 October 2026

The small self-contained Linux gate contains two saved benchmark panels (221
DNA cellular candidates and 856 RNA viral candidates), the actual VIHARMONY
module/configuration, immutable expected-output fingerprints, and 17 dependency-
free regression tests. Native stages are the real family default, unchanged
resume, genus, species, and restoration of the original cached family tasks.
Only the harmonizer and a saved-file copy process are imported; no searches,
setup, deployment, production sources or student caches are used. Linux PASS is
not claimed by local-only package tests. Expected outputs inherit the passed
3,068-candidate pre-cap comparison, including exact retention and uncapped audits.

The exact final archive was extracted and locally validated: all 17 tests and
all six Python replays (two panels by family/genus/species; eight workers) passed,
with all packaged files unchanged. Status is deliberately
`LOCAL_ONLY_NATIVE_PENDING`, not Linux PASS. Package SHA-256:
`210c7851c75dc509e28ba6dde59daa4a68cfe79f7247ab898dc4e4993bb0df23`.
Archive and complete shell instructions are under the primary repository's
`outputs/rank_cap_linux_20261018/` directory.

This original gate verified task-work outputs and native hashes/cache status,
not published files. Nextflow's documented `publishDir` overwrite default is
false during resume; the publication follow-up below addresses that correctness
gap. Do not silently treat a task/cache-only PASS as proof of publication refresh.
Reference: https://docs.seqera.io/nextflow/reference/process/directives/publish-dir

### Linux rank-cap review and publication adjustment — 18 October 2026

The uploaded `rank_cap_linux_review.tar.gz` passed all 17 regression tests and
all five native stages for 1,077 candidates. An independent local review checked
the exact pinned package identity, traces, original/restored task hashes,
actual rank/worker commands, receipt manifests and non-audit output digests.
First family tasks completed; unchanged/restore tasks were cached with original
hashes; genus/species recomputed only VIHARMONY while saved-input staging stayed
cached. The only log warnings were the expected unmatched `DATABASE_SETUP:.*`
selector, because the reduced fixture imports no setup processes. The uploaded
review does not contain published outputs, so it cannot certify publication.
Independent review: primary `outputs/rank_cap_publication_20261018/review_summary.json`.

Implemented in isolation: only the VIHARMONY `publishDir` adds `overwrite: true`
and `failOnError: true` plus explanatory comments. Scientific scripts and config
are byte-identical to the passed package; the module differs only in publication
settings/comments. No student/server run, source, cache or published file changed.
Scope proof: primary `outputs/rank_cap_publication_20261018/scope_verification.json`.

The v2 native gate additionally verifies the real published manifest and every
emitted output against the current task, byte-for-byte, after every completed
stage, including cached family restoration. The verifier never repairs or copies
published files. New negative regressions reject internally consistent old
bundles, mixed new-manifest/old-table states, and stale species files during a
cached family restore. This test fixture still imports no search/setup programs.
The publisher adjustment is not yet target-Linux verified or deployed.

The exact v2 archive was extracted and locally checked: **22 dependency-free
tests, all six two-panel family/genus/species Python replays, and eight workflow
contracts passed**. Scientific configuration/scripts remained byte-identical to
the passed v1 package. Local package status is `LOCAL_ONLY_NATIVE_PENDING`;
the real Nextflow publisher still requires the new Linux gate. The final archive
is `visum_rank_cap_publication_linux_20261018.tar.gz`, SHA-256
`d2c73fd0a21bb88b6a7c453843b3cc789592aa72b001f6cb798065667affeccf`.
The archive, complete shell runner and local verification artifacts are in the
primary `outputs/rank_cap_publication_20261018/` directory. Server execution uses
a new private test directory, not a student cache or output directory.

### Integrated ICTV Linux follow-up — 17 October 2026

The uploaded `integrated_ictv_linux_review.tar.gz` establishes **86 regression
tests passed**, with no errors or skips, and a complete **3,068-candidate Linux
replay passed**. All **213 expected taxonomy changes** match the independently
verified expectations; viral decisions, retained/provisional sequences and
one/eight-worker scientific outputs remain identical. This closes the Linux
complete-harmonizer replay gate, not deployment or the native cache gate.

The native fixture stopped before scheduling tasks because its test runner put
the global Nextflow `-C` option after `run`. Corrected only the test command
builder and added six local regression checks, including the exact reported
warning. A small native-only continuation package reads checksum-verified prior
passed results and the unchanged original package. It does not repeat the 86
tests or complete replay. Its four native stages (first run, unchanged resume,
manifest change, equivalent bundle-byte change) remain **pending Linux
execution**. Local offline continuation/preflight and checksum-rejection checks
do not establish native execution. Scientific source/module/reference hashes
are unchanged by this repair. No student run, species rule, rank cap, retention
rule or deployment is changed.

### Integrated ICTV native completion and main-entry-point gate

Integrated reconciliation gate update (17 October): the supplied
`integrated_ictv_native_review.tar.gz` passed all six command-fix checks and all
four native module stages. Independent trace/log review confirmed initial
completion, all-four-task unchanged resume, and two cached staging tasks plus
two harmonizer reruns after either reference-manifest or equivalent bundle-byte
changes. All task exits were zero; supplied logs contained no WARN/ERROR lines.
Original sources/outputs remained unchanged. The Linux complete replay and
native module/cache gates are therefore complete, superseding the pending
native status above. Deployment remains unperformed.

The next narrow fixture invokes the unchanged `visum_nextflow.nf`, automatic
config and real VIHARMONY module. All other imported modules become explicitly
test-only saved-file channel adapters inside a NEW scratch project: it checks
DNA/RNA evidence joins, viCAT context/locus/projection slots, vCONTACT3 groups,
default pinned reference/manifest staging, unchanged resume, and fail-fast
missing-manifest behavior. It does not test actual upstream tools, preparation,
readiness checks or fresh installations. Older saved panels lack separate
viCAT context artifacts; a documented header-only neutral table supplies that
slot, with no invented prediction rows. No production workflow/source is
rewritten. Actual main-entry-point Linux execution remains pending the new
package; local contracts and direct saved-input replays do not establish it.
The final extracted entry-point package passed 17 local contracts. Four direct
saved-input replays (DNA/RNA in saved and reversed evidence order) matched all
16 output entries per comparison. Only audit row/column order and gzip header
timestamps are normalized; audit values/multiplicities, other output bytes and
policy/reference metadata must remain identical.

### Resource policy and evidence

- No viSUM-imposed wall-time limits for analysis, utility/standardization, database construction or existing-asset validation. Use native `null`, not a large nominal duration.
- Keep all existing per-tool time parameters. Explicit duration overrides and native site process selectors remain available; external scheduler/queue limits may still apply or require a time request. Test harnesses retain their own bounded timeouts.
- viHARMONY's **32 GB fixed reservation** gives headroom above approximately 20.00 and 24.16 GiB sampled summed process-tree RSS in two large eight-worker replays. Shared pages can be counted multiple times; this is not a measured maximum or a guarantee for larger inputs/evidence. Increase `--harmonizer_memory` deliberately where needed; memory is not automatically multiplied by workers.
- CPU defaults, biological thresholds, evidence selection, taxonomy/retention logic, output schemas, database contents and cache management are unchanged by this resource-only update. A native local reservation is not an OS memory cap or a cross-run resource budget.
- Do not apply candidate source/config changes to an active student analysis. A changed workflow revision is not assumed to be cache-compatible with old student runs.

### Local verification — 15 October 2026

- Focused resource/TEsorter/viHARMONY workflow contracts: **24 passed, 58 subtests passed**.
- Full isolated-checkout suite: **438 passed, 15 skipped, 196 subtests passed** in 144.68 seconds. Skips are not passes; this Windows host has no native Nextflow runtime, and Linux/runtime-specific checks remain separate gates.
- Before/after fingerprints confirm that all existing `bin/`, `envs/`, `assets/` and `workflows/` files, the workflow entry point, wrapper, automatic config and viHARMONY module are unchanged by this update. The only process-module changes are the two viCAT setup time declarations; restoring each old declaration reproduces its pre-update file hash.
- Config verification permits only the 20 time defaults, harmonizer memory, SLURM time removal and explanatory comments. No decision, tool-selection, thread, database or output-schema setting is changed. New regression contracts cover all time parameters, every module time directive, profile limits, viCAT validation/build overrides and the 32 GB harmonizer request.
- `git diff --check` passed. No student analysis, production checkout, database, cache or published output was accessed or changed by this update.

### Open alpha gates, in order

Resource-test follow-up: the supplied v1 Linux report passed all twelve resource
contracts and six configuration checks (default, local, local_safe, local_full,
vicat_build, slurm). Docker printed an unrecognized `docker.pull` warning before
valid JSON; the harness failed to parse it and stopped before Singularity or
any executable probes. The isolated follow-up removes that ignored setting,
retains warnings explicitly while decoding JSON, and regression-tests the exact
uploaded console bytes. The v2 package selects only Docker, Singularity and the
three unstarted local probe cases. Those native checks remain pending; there is
no reason to rerun completed biological analyses or prior passing checks.

1. **Complete target-Linux verification of this resource-only update.** Local regression tests and source fingerprints passed as recorded above. Run the latest candidate on Linux to verify resolved no-limit/default/override behavior; portable text contracts alone do not establish runtime behavior.
2. **Approve reporting presentation before implementing it.** Design the slim summary, candidate-versus-final taxonomy scope, single-method species qualification, CheckV quality display and understandable review categories. Preserve raw decisions/evidence. Any retention-rule proposal remains separate, explicitly reviewed and replay-tested.
3. **Validate fresh Linux environments and setup.** Include the Parquet gate and TEsorter smoke path, real database preparation/reuse, small DNA/RNA and batch analyses, zero-call/bad-input behavior, and interruption/resume. The completed smoke test with existing databases does not close this gate.
4. **Triage upstream warnings and packaging.** Review missing `vclust`, vCONTACT3 realm/rank coverage limits, TEsorter partial-codon and Deep6 retracing warnings. Separate installation defects from upstream limitations; do not alter scientific inputs/models merely to suppress warnings.
5. **Finalize database delivery.** Keep prepared-path and source-builder routes. Verify actual release bytes/checksums, preserve provenance limits and source notices, and resolve distribution review before hosting viCAT bundles. VITAP's existing lifecycle needs validation, not another redesign.
6. **Freeze and benchmark one candidate.** Consolidate only approved changes, freeze source/config/environment/reference identities and evaluation endpoints, then finish the interrupted geNomad-derived benchmark. Keep student builds separately attributable and reuse valid saved evidence/cache where compatible.
7. **Package alpha.** Reconcile stale guides, add repeatable automated release checks, select licensing/citation metadata, audit distribution for private/large artifacts, and publish a versioned preliminary tag and release notes. No tag/DOI is asserted here.

Rust migration, a new resume manager, and wholesale database redesign are not alpha requirements. VirSorter2 contribution/ablation remains a separate evidence-based evaluation, not permission to drop a slow tool.

### Reconciliation-specific reporting follow-up

Restored taxonomy evidence can change the winning lineage or expose previously
hidden cross-tool conflicts; it is not merely cosmetic renaming. Candidate
taxonomy is still independent of primary viral retention. For `CELL_NEG_000206`,
GiantHunter supplies qualified viral evidence, vCONTACT3 supplies taxonomy,
and DeepMicroClass2/viCAT supply cellular evidence. The saved origin decision
is `likely_nonviral`; reconciliation does not turn it into a retained virus.
The recovered class conflict can leave only domain `Viruses` with the existing
`multi_method_consensus` label. Future human-facing reports must clarify the
rank at which agreement exists, without implying agreement on deeper lineage
or verified viral identity. Do not change that confidence policy in this task.

## Historical engineering progress — 13 October 2026

The following records what was known on 13 October; use the current status above for subsequent engineering completion and remaining release gates.

The isolated performance candidate at base commit `d4b90f3c235eae41b03a24ea4fbfd62b5febdb1c` now indexes evidence, RNA loci, and parent dispositions and supports bounded candidate multiprocessing. Existing decision functions and extracted candidate rules were checked against the frozen build; 138 policy/performance/workflow tests passed. Eleven saved panels (4,103 inputs, 2,576 candidates) reproduced the frozen scientific outputs. The candidate is not yet deployed to the shared student installation.

The server replay of completed sample `Tt_4_S2` used the recorded task inputs, evidence order, settings, and Python 3.11.15 environment. Both the updated one-worker and eight-worker runs passed all 14 output comparisons against the published results. Elapsed time was 165.61 versus 28.01 seconds (5.91-fold speedup); sampled summed process-tree RSS was 1.10 versus 1.54 GiB. These are single-run measurements for that sample, not a general scaling estimate.

The same tested source finished `BW_10.2_S7` in 202.53 seconds and `KB_5.2_S9` in 237.59 seconds with eight workers, sequentially. Their 14 output integrity checks passed, and saved inputs and sources remained unchanged. The original runs were still incomplete, so these two samples do not have a direct baseline-output comparison. Their sampled summed process-tree RSS was approximately 20.00 and 24.16 GiB; shared pages can be counted multiple times. The existing 4 GB harmonizer request is not an adequate scheduling budget for these large runs.

Student rescue output generation is finished; stopping the obsolete whale workflow and publishing the rescued files into the usual sample folders have not yet been confirmed. Keep the completion reports with the student provenance. The other student workflow remains on its original installation; no pull into that active checkout has been performed as part of this work.

## Current evidence and limits

The Tt_4_S2 stranded-whale metatranscriptome is a difficult, reportedly degraded postmortem RNA sample. Its viHARMONY outputs account for 77,559 input sequences: 412 primary retained viral contigs, 7,040 provisional viral contigs, 1,574 ambiguous records, 8,255 `likely_nonviral` records, and 60,278 sequences excluded during discovery. Of the primary contigs, 318 are shorter than 1 kb. These counts describe pipeline dispositions, not sensitivity or biological truth.

The present `likely_nonviral` rule applies to an ambiguous discovery candidate with cellular or plasmid conflict, no strong viral tool, and at most one qualified viral tool. It does not require regional host dominance or compare calibrated probabilities. In Tt_4_S2, all 8,255 such records have one qualified viral tool; 8,239 have a cellular conflict, 48 have a plasmid conflict, and 32 have both. The label therefore means a conservative non-retention decision, not proof that no viral sequence is present.

Two records illustrate the different evidence patterns. `NODE_1_length_23084_cov_168.148538_g0_i0` has qualified Deep6 viral support, cellular calls from CheckV and viCAT, and 22 viCAT cellular-supported loci with no viral-supported loci. `NODE_11285_length_392_cov_1.811912_g11189_i0` has one qualified viCAT viral call and one viral-supported locus opposed by a Deep6 cellular call. Both receive `likely_nonviral`, although the second is more accurately described as unresolved conflicting evidence on a short fragment. Across the sample, 2,034 `likely_nonviral` records have viCAT's `predominantly_viral` pattern. This does not establish that they are false negatives; their locus-level homology, host context, and competing evidence need review.

Taxonomy is calculated independently of viral-origin disposition. The outputs contain 290 `likely_nonviral` records with a strict species field; that field must not be presented as a final viral species call. Among the 412 primary retained records, 41 have a strict species assignment. All 41 are single high-tier VITAP assignments, and 39 are on contigs shorter than 1 kb. Their current `strict_single_high` basis is visible in metadata, but a named species on a short fragment can be mistaken for a confirmed species identity.

vConTACT3's reference-supported, length-eligible ICTV assignments may contribute formal taxonomy votes. Its synthetic `novel_` and `unplaced_` network groups are auxiliary project-group context, not formal ICTV taxa. None of the 412 primary Tt_4_S2 records credits vConTACT3 at its final reported taxonomic rank. The long, comma-delimited group field is difficult to interpret and can combine multiple group memberships.

## Candidate changes to test after student analyses

1. **Separate exclusion from biological certainty.** Preserve non-retention of weak, conflicting candidates while testing whether `likely_nonviral` should be reserved for records with convincing exact cellular or plasmid context. Place single qualified viral calls opposed only by generic cellular classification in an explicit conflicting-origin review state, or use a more accurate non-retained label. Do not assume all 2,034 viCAT-predominantly-viral records should be promoted. Stratify them by fragment length, viCAT viral and cellular loci, protein alignment strength and coverage, CheckV findings, and other viral evidence first.
2. **Gate the displayed final taxonomy by disposition.** Preserve all taxonomy votes and candidate-level annotations for audit, but distinguish candidate taxonomy on nonprimary records from final taxonomy on retained viral records. A strict species field on a `likely_nonviral` row must not be read as a final virus identification. Test downstream tables and FASTA selection for this distinction.
3. **Qualify single-method species names.** Keep `strict_consensus` separate from `strict_single_high` in reports. Describe the latter as a VITAP-supported closest-reference species hypothesis until raw alignment length, reference coverage, specificity among related species, and independent corroboration have been audited. Do not impose an arbitrary nucleotide-length cutoff or demote these calls without a saved-results comparison.
4. **Make auxiliary context legible.** Report vConTACT3 formal rank votes separately from project-group memberships, including database domain, group rank, status, and ambiguity. Expose CheckV's raw quality category and completeness evidence separately from viHARMONY's combined `checkv_complete_viral_contig_support` flag. Review queues should distinguish primary retained flags, provisional candidates, ambiguous conflicts, and likely-nonviral audit records.
5. **Add one slim, human-facing sequence summary.** Keep the 61-column integrated metadata and all tool-level evidence unchanged for audit. Add a compact TSV with one row per evaluated final sequence, the source and final IDs, length, primary-retained flag, plain-language origin status and its basis, viral supporting tools, conflicting tools, best supported taxon and rank, taxonomy basis (`multi_method_consensus`, `single_reference_hypothesis`, `exploratory`, or `unclassified`), candidate-versus-final taxonomy scope, provirus status, and concise review flags. If CheckV quality is included, use its actual category or `unknown`; do not convert a false combined support flag into an incomplete-genome claim. Sort or document filters so primary candidates are easy to find. Avoid copying long vConTACT3 synthetic-group strings into this view; retain them in the detailed output. The proposed labels and column set require user review before implementation.
6. **Integrate the validated viHARMONY performance candidate.** Indexing and bounded multiprocessing are implemented and tested in isolation; see the engineering progress above. The candidate exposes `--harmonizer_cpus`, defaults to one worker, and caps its request by `max_cpus`. The isolated candidate now requests 32 GB with no imposed wall-time limit; this is an adjustable starting reservation, not a bound. Finish the student handoff, integrate only verified changes, and validate the latest resource settings and scheduler accounting on Linux before deployment. Preserve deterministic ordering, exact decisions, audit fields, and output parity. A Rust rewrite is not required for the present alpha milestone; reconsider further implementation changes only if profiling reveals a remaining practical bottleneck.

## Engineering order and deferred work

The current gate order above supersedes the original six-step engineering plan. Native setup/run/resume and bounded multiprocessing are implemented; do not repeat that redesign. Resource defaults are updated in isolation, with Linux runtime verification still required. Reporting proposals above remain unapproved for implementation, and scientific retention/taxonomy policy must not change silently.

**VirSorter2 contribution analysis (deferred):** audit elapsed time/execution phases and replay affected discovery/refinement/harmonization with and without saved evidence. Measure rescued detections, lost corroboration, changed boundaries/taxonomy, nominal-negative retention and intact-genome recovery separately for DNA and RNA. Saved evidence cannot reveal sequences an alternative discovery configuration never sent downstream; disclose that limit or run only the missing downstream tasks when needed. Do not remove VirSorter2 or change defaults merely because it is slow.

## Evidence required before changing policy

- Obtain Tt_4_S2 standardized and raw viCAT, Deep6, CheckV, VITAP, and vConTACT3 evidence for representative records from each pattern above. The current viHARMONY outputs alone cannot establish whether the conflicted fragments are truly viral or cellular.
- Replay each proposed change against saved results before running discovery tools again. Report effects on primary retention, provisional and ambiguous counts, DNA and RNA positive sensitivity, intact-positive recovery, and original-label nominal-negative retention. Review putative viral-like negatives separately; do not silently relabel them as true positives.
- Benchmark viHARMONY scalability separately from discovery: record sequence and evidence-row counts, wall time, CPU utilization, and peak memory on small and million-sequence inputs. Compare any optimized or parallel implementation against frozen outputs before adopting it.
- Finish and independently review the interrupted geNomad-paper-derived benchmark with its frozen inputs and policy provenance. The benchmark is a separate validation step; do not infer its outcome from this student sample.
- Run a clean setup and end-to-end smoke test after policy and reporting changes, then freeze the release commit, configuration, database versions, and documentation. Keep the student runs attributable to their original build.

## Source of truth

This checklist records decisions to revisit; it does not replace run-specific outputs or commits. The Tt_4_S2 figures above come from `Tt_4_S2.sequence_disposition.tsv`, `Tt_4_S2.final_metadata.tsv`, `Tt_4_S2.review_queue.tsv`, and `Tt_4_S2.harmonizer_manifest.json` supplied for review on 7 October 2026. The local development checkout examined on that date was `dev/review-visum-scripts` at `0924706` with uncommitted work; it is not asserted to be the exact server revision used for Tt_4_S2.

The 15 October smoke completion/cache findings come from the uploaded `completed_smoke_review.X70dtJ.tar.gz`: final run `resume_20261014_203229`, original session `529b89df-9cd1-44d3-98ba-efeda26a9b3d`, traces/logs, both samples' outputs and manifests. Earlier Linux fixture evidence includes `cache_test_review_20261014_v3.tar.gz` and the preceding rerun reports. Neither source proves deployment or a clean installation of the latest candidate. Local verification above applies only to this isolated resource-default update.
