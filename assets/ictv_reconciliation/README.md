# Pinned ICTV reconciliation — isolated alpha integration

Integrated into viHARMONY and wired as native Nextflow task inputs in the
isolated alpha checkout only. Not deployed to server/student checkouts. The
standalone resolver remains available for independent validation. No automatic
release discovery/update feature is enabled. Viral-retention and
species-acceptance rules are unchanged; reference-compatible taxonomy evidence
can now contribute where naming/version differences previously excluded it.

The compressed bundle targets **MSL41.v1**, derived from the official
[ICTVdatabase](https://github.com/ICTV-Virus-Knowledgebase/ICTVdatabase) snapshot
`e83bc2c73a746165a209ba88d76ef8a8fb8635d2`. Its sidecar manifest records source
Git blob identifiers, SHA-256 digests, the target release, and bundle digest.

ICTV data are copyright International Committee on Taxonomy of Viruses and
licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
See `LICENSE.txt`. Changes made here: project the target tree, index historical
names/identities and official merge/split relationships, and include historical
immediate child-taxon relocations. This is a viSUM-derived prototype, not an
ICTV-endorsed classifier.

## What it does

`bin/ictv_reconciliation.py` resolves a dictionary of unprefixed formal
rank/name pairs. It returns the original dictionary, normalized taxonomy,
resolution status, target release, source release if supplied, and historical
mapping/proposal identifiers. It does not receive or change scores or evidence
strengths. Ancestors can be reconstructed; unsupported descendants are never
filled merely because the reference contains them.

Supported statuses:

| Status | Meaning |
| --- | --- |
| `unchanged` | Supplied names remain current; missing ancestors may be filled. |
| `translated` | Official relationships give a compatible current placement. |
| `ambiguous` | Several destinations remain; only shared ancestry is returned. |
| `partial` | An unknown name stops deeper recovery, or an unsupported higher prediction is flagged while compatible lower placement is preserved (see below). |
| `incompatible` | Supplied mapped claims cannot share a target lineage; no normalized claim is returned. |
| `no_named_taxonomy` | No named input taxa; this is not a naming mismatch. |

The current prototype supports the eight ranks already used by viHARMONY.
Intermediate ranks in the official tree are preserved as ancestry containers,
not added to the pipeline schema.

## Important limitations

- Supply `source_release` when it is known. When absent, exact current names
  use their current meaning, while obsolete names use historical observations.
  A name reused or redefined across releases may need source metadata to
  disambiguate it; the prototype cannot recover missing source-version facts.
- Splits, historical child relocations, and rank changes are not ordinary
  spelling substitutions. An old realm with multiple successor branches is
  ambiguous without context. A family promoted to an order supplies order,
  not an invented new family.
- Valid nomenclature does not establish correct biological classification.
  This does not validate a gene's species-discriminating power, change
  confidence tiers, or recalibrate an older classifier's reference clusters.
- Truly incompatible assigned lineages remain flagged. Unknown names are not
  matched approximately. Broader recovery is deliberately conservative.
- A missing higher rank in ICTV is not an affirmative competing placement.
  If every surviving, jointly compatible lower target lacks a supplied higher
  rank, that higher prediction is omitted from normalized taxonomy and flagged
  as `target_higher_rank_unassigned`. Its original name, mapping history, and
  supporting lower target IDs remain in the audit. Example: Naldaviricetes and
  Baculoviridae can be preserved without inventing their realm. This does not
  validate the biological accuracy or evidence strength of the lower call.
- The omission exception does not guess through unknown names, use a mixture
  of assigned/unassigned split destinations to discard an ancestor, or silently
  prefer a lower claim against an explicitly assigned conflicting ancestor.
  All known claim pairs are checked first, so missing ranks in a lower taxon
  cannot conceal an explicit conflict between two upper predictions.
  Rank-changing higher mappings are conservatively excluded from this exception.
- Historical paths marked as reversals or taxon recreations are not treated
  as equivalence mappings by this prototype.
- This first prototype uses a pinned export, not the public ontology API or a
  SQL server. No per-sequence internet requests are made.

## Standalone checks

From the isolated checkout:

```sh
python -B -m unittest discover -s tests -p 'test_ictv_reconciliation.py' -v
python -B bin/ictv_reconciliation.py \
  --bundle assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz \
  --source-release 39 \
  --lineage-json '{"genus":"Cystovirus","species":"Cystovirus phiYY"}'
python -B tests/audit_ictv_unassigned_ranks.py \
  --bundle assets/ictv_reconciliation/MSL41.v1.e83bc2c73a74.json.gz \
  --output PATH_TO_NEW_UNASSIGNED_RANK_AUDIT.json
```

The builder consumes already downloaded, exact pinned source files, verifies
their Git blob identifiers, and refuses to overwrite an existing bundle:

```sh
python -B bin/build_ictv_reconciliation_bundle.py \
  --source-dir PATH_TO_PINNED_OFFICIAL_EXPORT \
  --target-release MSL41.v1 --output PATH_TO_NEW_BUNDLE.json.gz
```

Required source files: `taxonomy_node_export.utf8.txt`,
`taxonomy_node_merge_split.utf8.txt`, `taxonomy_toc.utf8.txt`,
`taxonomy_level.utf8.txt`, and the official `LICENSE.txt`.

`tests/validate_ictv_reconciliation_replay.py` compares saved full evidence
audits, first requiring reproduction of the saved taxonomy baseline. It uses
the existing taxonomy engine unchanged, replacing only rank/name fields on
copies of evidence rows. It writes separate reports, never writes into the
saved results tree, and checks source result/FASTA/workflow hashes. Retention
is not recomputed by this taxonomy-only harness.

The initial benchmark check uses `--affected-only`: it reproduces and compares
records with a saved `invalid_configured_ictv_msl` rejection. Unaffected records
are explicitly counted as not replayed. Omitting that flag requests the slower
exhaustive comparison; do not describe a focused report as an exhaustive one.

Automatic updates and any species-rule or reporting-cap changes remain
separate future work. The previously distributed standalone Linux package is
a frozen pre-integration validation snapshot, not a deployment package.

## viHARMONY integration

The workflow defaults `ictv_reconciliation` to this pinned bundle and stages
both it and its JSON manifest as `path` inputs to VIHARMONY. Native Nextflow
hashing therefore sees changes to these reference inputs; there is no custom
resume manager. A direct harmonizer CLI invocation enables the component with
`--ictv-reconciliation PATH`. Omitting that argument preserves legacy behavior
for controlled comparisons. Setup/run never contacts ICTV to change a release.

The bundle's complete species lineages must match the configured `ictv_csv`;
a mismatching target fails rather than silently mixing ICTV releases. The
harmonizer manifest records the bundle and manifest checksums, target release,
and source revision. One resolver is cached per worker process. Original
evidence rows, scores, strengths, classification, retention adjudication, and
species voting rules are not changed. Only copied taxonomy columns enter the
existing taxonomy engine.

`*.taxonomy_reconciliation_audit.tsv.gz` streams original and normalized
rank/name dictionaries, resolution status, flags, and historical mapping
provenance. Full mode records every named eligible call; compact/none still
record changed/unresolved calls. Original full evidence audits stay unchanged.
The audit is published through the existing `*_audit.tsv.gz` output channel.
No reporting-rank cap is introduced by this integration.

`tests/replay_integrated_ictv_reconciliation.py` replays the complete frozen
and integrated harmonizers on original saved benchmark inputs. It requires
published baseline reproduction, compares every non-taxonomy metadata field,
primary/provisional FASTAs, mappings/disposition, and original evidence audit,
and independently checks every candidate's strict taxonomy against standalone
reconciliation plus the unchanged voting engine. Taxonomy-dependent review
flags/group context may change and are not mislabeled as retention changes.

The release-wide unassigned-rank audit enumerates every formal target taxon
with missing higher ranks. For each, it verifies the unchanged valid lineage
and separate synthetic challenges inserting a known label at each missing
higher rank. Distractors compatible with assigned upper ancestry must be
flagged while preserving the lower placement. If no compatible distractor
exists, the injected assigned-ancestor conflict must remain rejected. These
are synthetic safety tests, not measured biological predictions.
