# Pinned ICTV reconciliation reference

viHARMONY uses this bundle to reconcile supported taxonomy names with the
configured ICTV reference. It does not download taxonomy updates during setup
or analysis. The reference is pinned so that a run does not silently change
when ICTV publishes a new release.

## Identity and attribution

The compressed bundle targets **MSL41.v1**, derived from the official
[ICTVdatabase](https://github.com/ICTV-Virus-Knowledgebase/ICTVdatabase) snapshot
`e83bc2c73a746165a209ba88d76ef8a8fb8635d2`. Its sidecar manifest records the
source Git blob identifiers, SHA-256 digests, target release and bundle digest.
Keep the bundle and its manifest together.

ICTV data are copyright International Committee on Taxonomy of Viruses and
licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
See [LICENSE.txt](LICENSE.txt). viSUM projects the target tree, indexes historical
names and official merge/split relationships, and includes historical immediate
child-taxon relocations. This derived resource is not an ICTV-endorsed classifier.

## Interpretation

| Status | Meaning |
| --- | --- |
| `unchanged` | Names remain current; supported missing ancestors may be filled. |
| `translated` | Official relationships support a compatible current placement. |
| `ambiguous` | Multiple destinations remain; only shared ancestry is returned. |
| `partial` | Only part of the supplied taxonomy can be reconciled. |
| `incompatible` | Supplied mapped claims cannot share a target lineage. |
| `no_named_taxonomy` | No named input taxa were provided. |

The resolver never invents unsupported descendants. A taxon legitimately
unassigned to higher ranks in ICTV retains those gaps; for example, a supported
Naldaviricetes placement does not require an invented realm. Incompatible
assigned ancestors remain conflicts, and unknown names are not fuzzy-matched.
Splits, reused names and rank changes can remain ambiguous when an upstream
program does not supply its source release.

Reconciliation changes compatible taxonomy names, not evidence strengths or
viral-retention rules. A valid name does not establish biological accuracy or
justify a species call from a nondiscriminating fragment. The separate
`--harmonizer_max_rank` setting controls the maximum reported depth and defaults
to family; deeper predictions remain in the audits.

## Runtime use

The workflow stages the bundle and manifest as viHARMONY inputs, so ordinary
Nextflow caching sees reference changes. Complete species lineages must match
the configured `ictv_csv`; an incompatible pairing fails rather than silently
mixing releases. Do not replace either reference under an active analysis.

The harmonizer manifest records reference identities. Taxonomy reconciliation
audits preserve original and normalized names, resolution status, flags and
mapping provenance. Original evidence remains available.

See [interpreting results](../../docs/interpretation.md) and
[output files](../../docs/reporting-bundle.md).
