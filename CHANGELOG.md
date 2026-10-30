# Changelog

## 0.9.0-beta.2 — development candidate

This candidate promotes the approved RC3 implementation for the next testing
beta. It has not been tagged or published. Changes below are relative to beta.1;
the release-cleanup pass itself does not change RC3 biological decisions.

### Usage and installation

- Ordinary Nextflow entrypoint/config discovery and native resume.
- Database preparation is a separate `--setup` command; analysis reads prepared assets.
- Automatic, checksum-pinned MetaVR acquisition and frozen-vOTU metadata adaptation.
- Atomic named-file MMseqs membership export avoids `/dev/stdout` shard failures;
  failed export/compression does not publish an incomplete membership file.
- Opt-in `vicat_build_low_memory`: 8 viral-build threads, 160-GB task request,
  independent 96-GB DuckDB buffer setting, 192-GB local budget. Default build
  resources are unchanged. Full-input taxonomy parity passed at 113.8 GiB peak
  RSS; complete fresh source-build validation of this preset remains pending.
- Nonviral reference acquisition from the packaged 480-assembly list.
- Shared structural readiness checks and clearer failures for incomplete installations.
- TEsorter short temporary root, vCONTACT3 Parquet dependency/runtime check, and
  aligned Deep6 preparation/analysis environment.
- No shipped task wall-time limits; viHARMONY defaults to a 32-GB request.
- Bounded viHARMONY multiprocessing, configurable with `--harmonizer_cpus`.
- Compact parent-side evidence storage and ordered row-reference indexes reduce
  viHARMONY memory use without changing retention, taxonomy, TE routing or output
  content in saved-data validation. The new helper participates in Nextflow's
  cache fingerprint; existing completed searches need not be rerun for this change.

### Taxonomy, TE interpretation and output

- Pinned ICTV name reconciliation with original/reconciled evidence kept in audits.
- Family-default reporting ceiling with configurable maximum rank and uncapped audits.
- Reviewed TE-specific routing and distinct partial/mixed/localized evidence flags.
- Combined labelled viral export, separate TE export, rank-separated taxonomy,
  scoped CheckV quality, readable vCONTACT3 context and reason-specific review tables.
- Sample-prefixed filenames and explicit missing-value placeholders; no review FASTA.
- Current outputs are republished on resume; the manifest defines current inventory.

### Documentation and packaging

- User guides focus on setup, single/batch runs, resume, configuration and interpretation.
- The runtime branch excludes test suites, benchmark/replay tools, historical
  development notes and release-validation records. Earlier committed versions
  remain available in Git history.
- Existing public MIT license, motivation and beta.1 citation are preserved.
- Version/help metadata identifies beta.2; no beta.2 DOI is invented.
- Only helpers used by setup and analysis are shipped, alongside user guides,
  required reference inputs, provenance and license information.

### Migration and limitations

Use a separate installation and new analysis work/results directories when moving
from beta.1. Reuse compatible prepared databases after validation; do not replace
code or mutate environments under active runs. Results and TE destinations may
legitimately differ from beta.1, and old filenames are not the new primary interface.

Full fresh installation, assembled native real-tool smoke/resume, site-specific
SLURM use and the complete benchmark remain separate checks. The large viral
viCAT build still has a configurable 350-GB allocation, not a measured minimum.

## 0.9.0-beta.1

Original published testing beta. Its historical results, tag and Zenodo record
remain unchanged: [viSUMv0.9.0-beta.1](https://github.com/ViDaB-Lab/viSUM/releases/tag/viSUMv0.9.0-beta.1),
[doi:10.5281/zenodo.22774539](https://doi.org/10.5281/zenodo.22774539).
