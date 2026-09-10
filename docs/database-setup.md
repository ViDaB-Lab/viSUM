# Preparation-only setup

Run from the repository on Linux with the same config and database paths that
will be used for analysis:

```bash
bash ./visum -c visum.config --setup \
  --dbdir /shared/visum/databases \
  --tooldir /shared/visum/tools \
  --outdir setup_results
```

`--setup` prepares enabled database-bearing tools using the existing `--run_*`
flags. There is no FASTA, sample normalization, discovery, classification,
refinement or harmonization. Supplying input/prefix/type/prefix_many/indir is
rejected. It does not run viCAT's sample-specific reference-subset/hit preparation.
At least one database-bearing tool must be enabled.

Defaults are the existing config defaults: optional viCAT, VITAP and vConTACT3
are not implicitly enabled. DNA-only and RNA-only preparations both run if their
tool flags are enabled, since setup has no sample type. Disable unwanted tools
with the same `--run_<tool> false` settings used in analysis.

For a full-profile setup add `--run_vicat true --run_vitap true
--run_vcontact3 true` and supply both viCAT databases or their build inputs as
described in [viCAT preparation](vicat-database-preparation.md). Large local
builds require explicitly adequate `--max_memory`, disk and time limits; the
default 48 GB ceiling is not sufficient for the configured 350 GB viral build.

Preparations run sequentially to avoid overlapping their memory-intensive builds.
They use the same prepare modules, destinations, installers and validation as
analysis, and write per-tool records under `<outdir>/database_setup/`. A failed
preparation stops the dependency chain. Existing valid assets are reused.
Setup tasks are not restored from the Nextflow task cache, even with `-resume`:
the persistent installation is validated again without automatically rebuilding.

Update flags remain explicit. Setup does not enable VITAP/vConTACT3 updates or
check for newer releases by itself. See [VITAP selection](vitap-database-lifecycle.md).
No database is published/uploaded by this workflow.

## Limits of setup success

TEsorter has no separate database preparation process; setup prints a warning
when it is enabled. Its packaged references and runtime still require the small
analysis smoke test. Preparing the other databases also does not necessarily
create every distinct analysis Conda environment (notably Deep6).

The Conda cache remains configured separately from `dbdir` and `tooldir`. For
shared installations, set `conda.cacheDir` to a shared absolute path in your own
config. Use the same config for subsequent analyses. Paths should remain stable.

## Stronger completeness validation

geNomad validation checks the version, marker metadata, taxonomy files and all
three required MMseqs databases, including headers, indexes, nonempty data and
readability via `mmseqs dbtype`.

vConTACT3 validation follows the selected release's manifest for the requested
domains (`both` means prokaryotes and eukaryotes). It checks reference sequences,
gene/genome tables, all five runtime identity levels, their MMseqs databases and
cluster mappings, and the VOGDB profile database. Paths must stay within the
release directory. Missing domains, mismatched release keys, missing/empty files
and unsafe paths are rejected. Split MMseqs data files are supported.

These are structural/readability checks, not proof of biological accuracy or
cryptographic authenticity of every asset. Full-size upstream database tests
remain necessary.

Source contracts: [geNomad v1.12.0](https://github.com/apcamargo/genomad/blob/v1.12.0/genomad/database.py)
and [pinned vConTACT3](https://bitbucket.org/MAVERICLab/vcontact3/src/d57ae81e0be5ab57a8a14692a4cb8ab255ce2d52/vcontact3/databases.py).

## Routing regression test (no database downloads)

On Linux with Nextflow on PATH:

```bash
python -m unittest discover -s tests -p 'test_setup_nextflow_integration.py' -v
```

This uses temporary fake prepare implementations and verifies the trace contains
exactly the twelve selected preparation processes and no analyses. It tests
Nextflow routing, not the real installers. Then test a real clean setup, repeat
setup to verify reuse, and run small DNA/RNA analyses against those installations.

## Isolated Linux installation test

```bash
NEXTFLOW_BIN=/absolute/path/to/nextflow bash bin/smoke_database_setup.sh
```

This runs the routing test first, then installs the ten non-viCAT preparations
(including VITAP and vConTACT3) in a new `visum-clean-setup.*` directory beside
the repository. Databases, tool installations, Conda environments, Nextflow work
and outputs are isolated from the production installation. It performs a second
setup invocation against the same assets and checks both traces contain exactly
those ten preparations, with no analysis processes. Inspect the reuse log and
per-tool provenance to confirm existing assets were reused, not rebuilt.

Defaults are 16 CPUs and a 128 GB per-task memory ceiling; override with
`SETUP_CPUS` and `SETUP_MEMORY`. Set `SETUP_TEST_PARENT` to an existing directory
on a filesystem with room for these large downloads and builds. This is a real,
potentially long and disk-intensive installation test, not a small fixture test.
Conda must be on PATH. viCAT source builds require separate explicit inputs and
are excluded; TEsorter's analysis runtime is also not tested. On failure, preserve
the printed test directory and logs for diagnosis. The runner never deletes it.
