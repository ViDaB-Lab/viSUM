# Beta.2 local release review

Version: `0.9.0-beta.2`. Intended public tag: `viSUMv0.9.0-beta.2`.
This is a local review copy; no tag, push, GitHub release, Zenodo publication or
server deployment is performed by its preparation.

## Source identity

The approved scientific implementation is taken from frozen RC3:
`visum-alpha-rc3-20261020.tar.gz`, SHA-256
`af47611650ef39c160ba8ab85d0ee48344da28e24ec1cd99fd07dce175854868`.
The public license and documentation context were reviewed against GitHub main
commit `be1872df6056f26742dbdf6f38d3961d50ba2da2`.

The review changes documentation, CLI help, release metadata and packaging only.
The scientific helpers, modules, environments and reference bytes must match
RC3. The workflow may differ only in `visumHelp()`; `visum.config` may differ only
in manifest author/version and comments. That boundary is checked explicitly.

The MIT license and existing beta.1 citation are preserved. The beta.1 DOI is
not claimed as a beta.2 DOI. Machine-readable citation metadata uses the existing
published author initials rather than inventing names, ORCIDs or a new DOI.

One-off student-specific recovery scripts and unreferenced legacy fix scripts/
environment copies are excluded. General replay helpers required by retained
tests remain. No student sequence data or installed database is included.

## Evidence and limits

Prior Linux component checks cover setup/cache behavior, harmonizer replay,
taxonomy reconciliation, reporting caps, TE routing and reporting publication.
The combined MetaVR source preparation/reuse check passed on Linux on 21 October
2026: 42 tests, no failures/errors/skips/warnings, approved metadata checksum,
1,409 substitutions, zero non-vOTU changes, unchanged source/cache identities,
and no network/external-build calls. See `release/metavr-linux-validation.json`.
Historical RC summaries are retained as dated evidence, not current status files.

The beta.2 package's own local tests and source/packaging checks are recorded in
`release/validation-summary.json` when the review archive is generated. A local
unit pass is not a native Linux end-to-end pass. Nextflow/Java execution is not
available on the local Windows review host.

Still required before stronger installation/performance claims:

- Full fresh Linux setup and reuse with the selected databases and captured builds.
- Small actual-entrypoint DNA/RNA/batch analysis and resume on the assembled beta.
- Full evaluation with the frozen truth/endpoints and external reference identities.
- Site-specific SLURM validation on the target cluster.

## Publication handoff

Review the candidate, then promote the selected tree to a clean release branch
based on current public main. Preserve any later public-only changes rather
than blindly copying an old checkout. Commit and validate the exact source to be
tagged, then create an annotated beta.2 tag and GitHub prerelease after explicit
approval. Rebuild release provenance against that actual commit; this local
archive deliberately identifies itself as an uncommitted review snapshot.

Archive beta.2 as a new version linked to the existing Zenodo record, not a
replacement for beta.1. Record its actual DOI when assigned. Distribution of a
large prebuilt viCAT database is separate and requires source-term review.
