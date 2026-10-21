# Alpha release candidate 3

Historical RC3 scope record. The current local beta.2 review is described in
[beta.2 release review](development/beta2-release-review.md). The combined
MetaVR source-preparation/reuse check has since passed; this record preserves
the state at the original RC3 freeze rather than replacing its archive.

Identity: `visum-alpha-rc3-20261020`. This is a local, private validation candidate,
not a public release/tag, deployment, installation certification or benchmark result.
`SOURCE_MANIFEST.json` identifies every source/config/environment/reference byte;
`checksums.sha256` inventories packaged files. The Git base alone is insufficient
because approved changes remain uncommitted. No commit, push or student migration
is performed by packaging.

RC2 adds automatic nonviral source acquisition from the unchanged 480-assembly
inventory and explicit refresh of the nonviral setup-status TSV on resume.
The native Linux component gate passed: 43 tests with no skips; a real NCBI
assembly's 5,572 proteins classified; a seven-protein six-class synthetic database
built through the actual nonviral setup/module/builder and reused with unchanged
contents/modification times and fresh published status. The source server was
offline during native reuse; the live-source reuse check also blocked network.
The fixture substitutes only its download host, one-entry manifest and viral
dependency process. It is not proof of a full 480-assembly or full-tool install.
See `release/nonviral-linux-validation.json` for the reviewed evidence identity.
The original RC1 and component-test tarballs remain unchanged. No viral retention,
taxonomy, TE routing, clustering policy, truth label or benchmark endpoint changes
are included in RC2.

RC3 adds the robust bounded/resumable MetaVR downloader and approves the exact
portal source pair using completed native Linux evidence: all 369,319,989 ordered
protein IDs/sequences match, and all 24,435,662 builder-consumed metadata rows
match after the reviewed 1,409-assignment frozen-vOTU patch. Source, prepared
metadata, patch and evidence checksums are pinned together. The biological
builder/harmonizer rules remain unchanged. Existing manual/prebuilt setup works.
The original RC1/RC2 archives are unchanged. The combined actual source-preparation
and unchanged-reuse Linux gate is still pending; full database installation and
fresh-server smoke/benchmark gates are also NOT established by source equivalence.

The package includes current approved isolated runtime source, modules, environment
specifications, config, workflows, curated assets, documentation and source tests.
It excludes databases, installed software, Conda caches, analysis outputs, student
sequences, private worktrees, `.git`, editor state and legacy pipeline copies.
Linux shell files have LF line endings and executable modes. Both checkout and
packaged hashes are recorded when shell line endings differ. The archive uses
fixed entry order/timestamps/ownership and a fixed gzip timestamp; repeating the
build from identical inputs must produce the same SHA-256.

Packaged reference assets (ICTV CSV and reconciliation bundle/manifest) and all
environment YAMLs are frozen by checksum. These YAMLs are NOT solved dependency
lockfiles: exact Conda package builds must be recorded by the clean installation.
Large external tool databases are NOT included, silently frozen, or redistributed.
Before the benchmark, record their actual selected releases, preparation manifests
and checksums on the Linux host using `release/external-identities.template.json`.
Its null identities are deliberate outstanding gates, not reproducibility claims.
The already reviewed viCAT nonviral source inventory is included under assets;
that alone does not establish redistribution clearance or every remote DB byte.

## Closed component checks

- Indexed/bounded-multiprocessing harmonizer and documented resource defaults.
- Native Nextflow setup/run/resume design and cache regression components.
- Automatic nonviral acquisition and tiny native setup/reuse/publication gate
  (full-size reference acquisition/build and fresh full installation still pending).
- Pinned ICTV reconciliation, family-default reporting cap, and TE routing.
- Reporting Linux v2: 23 tests, 1,211 saved candidates, exact nested publication,
  six cached tasks on unchanged resume, helper-only invalidation of three
  harmonizers, and restoration of all six original cached tasks. Scientific
  results remained identical to the pinned expectations at every stage.

These component checks do not prove installation of this assembled package.

## Next gates, in order

1. Verify archive/source hashes; run self-contained source/evaluator tests from
   the extracted candidate. Review remaining upstream packaging/runtime warnings.
2. Fresh private Linux installation: separate code/env/cache/work/output paths;
   run documented `--setup` with explicitly supplied/buildable viCAT inputs and
   capture environment/database identities. Test database reuse and useful
   failures for invalid/incomplete assets. No use of student task caches.
3. Small real-tool DNA, RNA, and mixed-batch analyses through the main entrypoint;
   validate manifests, outputs, short-path TEsorter, Parquet support, zero-call
   inputs, malformed/duplicate inputs, interruption and native resume.
4. Verify existing prepared benchmark inputs. Review disk/time/resource budgets,
   capture full commands and external references, then launch the unchanged
   scientific preset in `release/benchmark-scientific.params.json` through native
   Nextflow. Supply installation/input/output/resource paths in a separate site
   or parameter file. No benchmark is launched by this release pass.
5. Run the frozen evaluator endpoints plus the comparative/publication protocol;
   report strata separately, preserve truth, and disclose training/reference
   overlap. Fix genuine defects and rerun affected checks before final tagging.
6. Resolve root-project licensing/citation and dependency/database distribution
   notices, reconcile final docs, then authorize a public alpha tag/release.

Do not reuse the historical `genomad_final_bundle` launcher as the current RC:
it snapshots older code, reuses installed environments, and imposed old viCAT time
limits. Its already prepared DATA remain valid if their verifier passes. The
current workflow is ordinary Nextflow with separate setup; no new resume manager
or Rust migration is required for this milestone.
