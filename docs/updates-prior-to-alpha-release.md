# Alpha development checklist — historical ledger

For the current beta.2 review status, see [release review](development/beta2-release-review.md).

Updated 20 October 2026. Scope: isolated alpha candidate only. No student source,
running process, database, analysis cache or published result is migrated.
Older development notes are preserved in `history/pre-alpha-through-20261018.md`;
their pending labels are historical and do not override this checklist.

## Closed components

- [x] Indexed viHARMONY and bounded multiprocessing, saved-output equivalence.
- [x] Native Nextflow launch/resume, separate database preparation, cache tests.
- [x] Defaults/readiness simplification and resource-policy Linux checks: no
  shipped program wall-time limits; configurable 32 GB harmonizer reservation.
- [x] TEsorter short temporary directory and vCONTACT3 Parquet dependency repair
  in the isolated candidate (fresh environment installation remains untested).
- [x] Pinned integrated ICTV reconciliation, including native entrypoint check.
- [x] Family-default reporting cap; deeper predictions preserved in audits.
- [x] Approved TE routing and clarified localized/record-level evidence flags.
- [x] Reporting bundle/publication/cache gate CLOSED by Linux v2 report:
  23 regression checks; 1,211 saved candidates; 107 published files per stage;
  first execution 6 completed; unchanged resume 6 cached; helper-only change
  3 staging cached + 3 harmonizers completed; restored helper 6 original cached.
  Actual hash inputs include the dependency fingerprint. Retention, taxonomy,
  TE routing and expected reporting values matched at every stage.
- [x] Automatic nonviral source acquisition and tiny native Linux setup/reuse
  gate: 43 tests, no skips; seven synthetic proteins through the real nonviral
  build; unchanged database files/mtimes on resume; current published setup status;
  one real NCBI assembly with 5,572 classified proteins and offline source reuse.
  RC2 includes this addition; the original RC1 is unchanged. The full 480-assembly
  acquisition/build is still part of the fresh installation gate below.

## This consolidation pass

- Manifest-driven readers support flat historical and nested current artifacts,
  validate checksums and reject ambiguous/mixed layouts.
- Frozen endpoint definitions distinguish retained viral decisions, broader
  provisional-inclusive viral exports, TE exports and TE-only input fragments.
  Supported confidence alone never creates a primary call. Regions deduplicate
  to the original input fragment; original truth labels remain unchanged.
- Prepared-benchmark evaluator keeps source-routed and both unresolved route
  strata separate. Missing outputs fail instead of becoming negatives. It does
  not invent absent OFF comparisons or pool paired arms.
- A checksummed, reproducible private RC snapshot freezes available source,
  config, environment specifications and packaged reference assets. External
  database identities and resolved environment builds still need host capture.
- Consolidated the already approved geNomad composition parameter from the
  preparation branch: default `auto`, validated CLI forwarding and metadata;
  the fixed benchmark preset selects `metagenome`. This is configuration wiring,
  not a viHARMONY rule or a newly tuned benchmark threshold. Main-entrypoint
  Linux validation of the assembled candidate remains required.

See `release-candidate.md` and `benchmark-evaluation.md` for exact scope/limits.

## Outstanding gates before public alpha

- [x] MetaVR metadata compatibility patch: native Linux export of all 1,409
  changed UViGs; all 24,435,662 normalized metadata records EXACT_MATCH after
  private-copy patching; every non-vOTU field unchanged and both sources untouched.
  Versioned/checksummed patch integrated into the shared source comparison/setup
  path, with provenance and immutable offline reuse; isolated offline regressions.
  This is metadata compatibility only, NOT full-source approval or deployment.
- [x] Native Linux check of the integrated metadata preparation/reuse path using
  existing metadata, with no protein access, network requests or database build.
- [x] Complete original/portal protein comparison: all 369,319,989 ordered
  first-token IDs/exact sequences match; compressed source identities pinned.
  Combined with native patched-metadata equivalence, approve the exact source
  profile and checksum-bound native evidence. No biological policy changed.
- [x] Native Linux combined source-preparation/reuse check through the actual
  approved setup helper, using existing portal files with network and repeat
  protein scans forbidden. Passed on 21 October 2026: 42 tests, no failures,
  errors, skips or warnings; cache/provenance hashes and timestamps unchanged.
- [ ] Fresh Linux installation, full selected database setup/reuse (including the
  480-assembly nonviral reference if viCAT is selected), runtime/environment identities.
- [ ] Latest assembled main-entrypoint real DNA/RNA/batch smoke and resume tests,
  zero-call/invalid-input cases, dependency/warning review.
- [ ] Verify existing prepared geNomad-derived input identities on the server.
- [ ] Freeze external database identities, benchmark resource/config manifests,
  comparative native-tool protocol, training/reference overlap disclosure.
- [ ] Full benchmark: primary/broad detection, TE export rates, parent-versus-
  intact-fragment recovery, taxonomy/rank accuracy where supported, resource use.
- [ ] Public packaging/licensing/citation and database/source redistribution review.
- [ ] Final approved alpha tag/release (not created by this pass).

Deferred, not new alpha requirements: Rust rewrite, a custom cache manager,
wholesale database redesign, unreviewed species rules, or removing VirSorter2
solely because it is slow. Any later biological-policy change requires a new
approval, candidate identity and validation; do not tune on benchmark outcomes
and then present them as independent validation.
