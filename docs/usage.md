# Install and use viSUM

[Documentation](README.md) · [Setup](database-setup.md) · [Outputs](reporting-bundle.md)

## Requirements and installation

Use Linux with Java, Nextflow, Conda and Mamba available on `PATH`. viSUM uses
Nextflow to create each tool's Conda environment. Native Windows execution is
not supported by this release's validation. The development runtime is Nextflow
26.04.6; use that version when reproducing the tested configuration.

The current beta.2 review archive can be extracted directly:

```bash
tar -xzf viSUM-v0.9.0-beta.2-review.tar.gz
cd viSUM-v0.9.0-beta.2
nextflow run visum_nextflow.nf --help
```

After the beta.2 tag is published, the equivalent repository installation is:

```bash
git clone --branch viSUMv0.9.0-beta.2 https://github.com/ViDaB-Lab/viSUM.git
cd viSUM
```

That tag is not created by the local review package. Keep the installation in
a stable directory with enough disk space for references and environments.

## Prepare databases once

```bash
nextflow run visum_nextflow.nf --setup
```

This prepares the enabled default programs, including DNA- and RNA-specific
assets. It does not analyze samples. Databases default to the installation's
`databases/`, software/model bundles to `tools/`, and Conda environments to
`.conda/`. No directory overrides are necessary for those defaults.

viCAT is off by default. To prepare its references too, see the large-build
resource warning in [viCAT setup](vicat-database-preparation.md), then use:

```bash
nextflow run visum_nextflow.nf --setup --run_vicat true -profile vicat_build
```

Use the same enabled tools and database paths during setup and analysis.
Existing valid databases are reused. Analysis does not download, build or
automatically update databases. TEsorter has no separate database setup task;
its packaged references are used when its analysis environment is created.

## Run a single sample

Input must be an assembled nucleotide FASTA, not raw reads. Use `dna` for DNA
assemblies and `rna` for transcript assemblies. RNA mode can detect evidence of
DNA and RNA viruses represented in the transcripts.

```bash
nextflow run visum_nextflow.nf \
  --input /data/contigs.fasta --prefix sample01 --type dna \
  --outdir results/sample01
```

For RNA, change `--type dna` to `--type rna`. Sample prefixes may contain letters,
numbers, periods, underscores and hyphens. Use unique prefixes in a batch.
FASTA headers and sequences must be nonempty. The first whitespace-delimited
header token is retained as the original ID; viSUM assigns its own unique
normalized IDs. Use unique original IDs to keep downstream label matching clear.

Enable prepared viCAT references during analysis with `--run_vicat true`.
Disable any unwanted program with its `--run_<tool> false` flag; its evidence
will then be absent from the final integration.

Nextflow uses the shipped `nextflow.config` and `visum.config` automatically.
The optional `bash ./visum ...` command forwards arguments to the same workflow;
it is not a separate installation or resume system.

## Run several samples

Create a comma-separated `samples.csv` with the header shown below. FASTA paths
are resolved against `--indir`:

```csv
prefix,type,fasta
dna_sample,dna,dna_contigs.fasta
rna_sample,rna,rna_contigs.fasta
```

```bash
nextflow run visum_nextflow.nf \
  --prefix_many samples.csv --indir /data/assemblies \
  --outdir results/batch01
```

Do not combine single-sample input options with batch input options. The
read-only preflight checks sample rows, duplicate prefixes, missing/empty input
files and required assets. FASTA contents are checked during normalization.

## Resume an interrupted analysis

Repeat the original analysis command from its original launch directory,
adding `-resume`:

```bash
nextflow run visum_nextflow.nf \
  --input /data/contigs.fasta --prefix sample01 --type dna \
  --outdir results/sample01 -resume
```

Preserve both `.nextflow/` in the launch directory and the task `work/` directory
(or the path selected with `-work-dir`). Results alone are not a resumable cache.
Keep inputs, references and settings stable; Nextflow decides which completed
tasks can be reused.

If you have launched multiple runs or setup commands from that directory,
select the intended analysis explicitly:

```bash
nextflow log
# Repeat the analysis command with: -resume RUN_NAME
```

You may launch a workflow by its absolute path from a dedicated analysis
directory, for example `nextflow run /opt/viSUM/visum_nextflow.nf ...`. This keeps
each analysis's history and work separate. Do not change the installation used
by an active run. Start a new work/output directory when upgrading versions;
cross-version reuse of beta.1 task caches has not been established.

## Set resources or use a shared installation

Defaults are an 8-CPU/64-GB local scheduling budget. Override them only when
needed, for example `--max_cpus 16 --max_memory '128 GB'`. Each tool still has
its own request; a smaller budget does not make a large task fit.

For parallel viHARMONY, add `--harmonizer_cpus 8` when that CPU budget is
available. Its default memory request is 32 GB and does not automatically scale
with workers; large inputs may need more memory. See [resource settings and
SLURM examples](resources.md).

For shared assets, set `--dbdir /shared/visum/databases` and, if needed,
`--tooldir /shared/visum/tools`. Conda cache location is a separate
`conda.cacheDir` setting in a site configuration. Supply the same site config
with `-c /path/site.config` for setup and analysis. These are optional overrides,
not required arguments.

You can also keep analysis parameters in a JSON file; copy and edit
[example-analysis.json](example-analysis.json), then run:

```bash
nextflow run visum_nextflow.nf -params-file /data/sample01.json
```

## Read the results

For the single-sample example, final files are under
`results/sample01/sample01_results/viharmony/`. Tool-specific raw/standardized
outputs remain alongside `viharmony/` in the sample result directory.

Start with `<sample>.viral_candidates.summary.tsv` and its matching FASTA, then
use `<sample>.taxonomy.tsv` for a rank-separated taxonomy table. Retained and
provisional candidates are distinguished by `candidate_status`; the FASTA
contains both. TE candidates are separate, and review files contain tables only.

See [interpreting results](interpretation.md) for the biological meaning of the
labels and [reporting bundle](reporting-bundle.md) for folders and exact columns.
Run reports and trace are under `<outdir>/reports/`; you do not need to supply
custom `-with-trace` or `-with-report` paths for the configured reports.

## Get help

Run `nextflow run visum_nextflow.nf --help`, consult
[troubleshooting](troubleshooting.md), or report an issue with the version,
redacted command, relevant task logs and expected behavior. Do not include
confidential sequences or credentials.
