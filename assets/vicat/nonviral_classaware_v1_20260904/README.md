# Frozen viCAT nonviral reference inventory

This directory documents the existing `nonviral_classaware_v1` database completed
on September 4, 2026. It is **not** a rebuilt database, new classifier policy,
downloadable database bundle, or claim of redistribution clearance. Post-RC1
setup uses `assembly_manifest.tsv` as its default automatic nonviral download
inventory. The manifest's bytes are unchanged; the historical audit summaries
remain documentation only and are not classification inputs. Analysis still
uses a separately prepared database and does not install/download references.

## Authoritative snapshot

- `assembly_manifest.tsv`: 480 versioned RefSeq assembly accessions with the
  cellular group, organism name and taxid recovered from every row of the
  SHA256-verified source catalog. Compatible with the existing source builder.
  Taxids are source annotations, not a claim of 480 distinct species.
- `provenance.json`: source hashes, recovered counts, historical script identity,
  build versions/parameters, database hashes and explicit verification limits.
- `build_evidence/`: original small build records preserved byte-for-byte from
  `nonviral_classaware_v1_review_20260904_172747.tar.gz`. The classification
  summary with **3,731,712 included proteins and zero exclusions** matches the
  source catalog and build. The database contains 2,708,478 representatives.

The builder and metadata helper match their recorded SHA256 values at Git
commit `d8805bd8e83a38aff0c3d97f2b2306a4b74f97f1`. This proves the identity of those
two files, not the entire checkout/environment used for the historical build.
The source protein sequences and feature-table files are not included here;
accession versions alone do not freeze subsequent annotation changes.

The recovery helper `tests/recover_vicat_nonviral_snapshot.py` reads the original
source catalog, build-review archive and standalone summary and emits proposed
text assets as JSON. It does not download, build, install or modify databases.
Its source hash guard deliberately restricts it to this historical snapshot.

## Unmatched evidence, not release/build input

`unmatched_evidence/classification_summary.tsv` preserves the separately supplied
summary with 3,731,308 included and 404 excluded proteins. It does not match this
database. Its ordering relative to the build is unknown; do not silently replace
the build-matching summary or infer that the database is corrupted.

The 404-protein difference comprises 64 unplaced cellular, 150 plasmid, 186
plastid, three mitochondrial and one shared-reference proteins. The complete
catalog includes special replicon mappings and 64 no-feature-mapping fallbacks.
No classification rules or biological interpretations were changed in this audit.

## Remaining provenance and release checks

Original assembly-selection criteria, download dates/commands and per-assembly
protein/feature-table hashes and annotation revisions were not recovered. Do not
invent them. Before hosting, inspect applicable notices for the exact sources
and verify the actual database files against the recorded checksums. Large
database artifacts absent from the review and current server bytes have not
been verified here. Keep the existing builder as a fallback; this inventory
does not guarantee that today's NCBI downloads will rebuild identical bytes.
