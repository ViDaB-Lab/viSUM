# Prepare and reuse databases

[Documentation](README.md) · [viCAT databases](vicat-database-preparation.md)

## Default setup

From the viSUM installation directory:

```bash
nextflow run visum_nextflow.nf --setup
```

Setup prepares enabled database-bearing programs sequentially, including the
DNA- and RNA-specific tools. Do not supply FASTA, prefix or molecule-type options.
It creates the required preparation environments and validates the assets.
Ordinary analysis only reads prepared assets; it never runs global installers.

Locations are already configured:

| Contents | Default location |
| --- | --- |
| Managed databases | `<installation>/databases/` |
| Software/model bundles | `<installation>/tools/` |
| Conda environments | `<installation>/.conda/` |
| Setup metadata | `<outdir>/database_setup/` (`outdir` defaults to `results`) |

Override `--dbdir` or `--tooldir` only when using another storage location.
For an alternative Conda cache, set `conda.cacheDir` in your site config.
Use the same asset configuration for setup and analysis.

## Include optional viCAT

All default programs are enabled except viCAT. To prepare those defaults plus
the viCAT viral and nonviral databases:

```bash
nextflow run visum_nextflow.nf --setup --run_vicat true -profile vicat_build
```

The approved MetaVR sources and the packaged nonviral assembly selection can be
downloaded automatically. No manually collected source files are required for
that route. **Check available memory and disk first:** the current viral build
requests 350 GB, and the profile budgets up to 32 CPUs/400 GB. This is a
configured allocation, not a measured minimum RAM requirement. The large build
does not describe the resources needed to use an already-prepared database.
See [viCAT preparation](vicat-database-preparation.md) for inputs and alternatives.

For the opt-in smaller preset, replace `-profile vicat_build` with
`-profile vicat_build_low_memory`: 8 viral-build threads, 160-GB task request,
96-GB DuckDB buffer setting and a 192-GB local budget. Full-input taxonomy
parity passed; complete fresh source-build validation of this preset is pending.
This is not a 96-GB or 128-GB machine-support claim. Both presets reuse valid
databases without reconstruction.

## Select tools and reuse existing assets

Use the normal `--run_<tool>` switches in setup, for example
`--run_virsorter2 false` to omit VirSorter2. If disabled during setup, prepare its
assets before enabling it in analysis.

To use an existing database outside the managed root, pass its per-tool path,
such as `--genomad_db /shared/genomad/genomad_db` or
`--vicat_db /shared/vicat/viral`. viCAT also needs `--vicat_nonviral_db`.
The same paths must be available from every execution host, including cluster nodes.

Rerunning the setup command validates and reuses valid installations. You do not
need a separate recovery launcher or `-resume` to reuse a completed database.
Preparation checks still execute even with `-resume`; an invalid database is not
silently overwritten. Interrupted work may be reused when its source/build
identity matches. Keep the existing private build files until a failure is understood.

## Updates are explicit

Setup does not automatically replace a valid database with a newer release.
VITAP and vCONTACT3 have explicit update options:
`--vitap_update_database true` and `--vcontact3_update_database true`.
Perform updates in setup, not analysis. For reproducibility, keep previously
used references and use a new destination/release label for changed inputs.
See [VITAP selection](vitap-database-lifecycle.md).

The ICTV rank table and reconciliation bundle shipped with viSUM are pinned
separately. A VITAP or vCONTACT3 reference update does not automatically update
that target taxonomy or its reconciliation bundle.

## What successful setup means

Setup checks required structure, identities and native readability where
implemented. It is not proof that every reference is biologically correct or
that all analysis runtimes work on a new machine. Run a small DNA/RNA analysis
before committing a large dataset.

TEsorter has no standalone database preparation task. Its packaged references
are checked when its analysis environment runs. This warning alone is not a setup
failure. viSUM imposes no task wall-time limits, but cluster queues may impose
or require one. See [resources](resources.md) and [troubleshooting](troubleshooting.md).
