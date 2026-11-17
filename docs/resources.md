# CPU, memory and execution settings

[Usage](usage.md) · [Setup](database-setup.md)

## Local execution

The default local scheduling budget is 8 CPUs and 64 GB. For setup with viral
viCAT enabled, the default memory budget is 192 GB. Set `--max_cpus` and
`--max_memory` to the resources available for this invocation, for example:

```text
--max_cpus 16 --max_memory '128 GB'
```

Tool-specific `--*_cpus` requests are capped at `max_cpus`. Each task retains
its own `--*_memory` request; viSUM does not shrink that request to fit a smaller
budget. Nextflow's local executor schedules against the invocation's budget.
These settings are not OS-enforced limits and do not reserve resources against
other users or independent Nextflow runs.

| Profile | CPU budget | Memory budget | Purpose |
| --- | --- | --- | --- |
| `local` / default | 8 | 64 GB | Ordinary setup and analysis |
| `local` / default with viCAT setup | 8 | 192 GB | Lower-memory viCAT source build |
| `local_safe` | 8 | 64 GB | Explicit conservative local preset |
| `local_full` | 32 | 128 GB | Larger dedicated local budget |
| `vicat_build_high_power` | 32 | 400 GB | Opt-in larger viCAT source build |

The default viral build task requests 8 CPUs and 160 GB RAM, with a separate
96 GB DuckDB buffer setting. The high-power profile requests 32 CPUs and
350 GB RAM. Its DuckDB buffer follows the task request unless overridden.
The buffer setting covers DuckDB-managed memory. Peak process RAM can be
higher. Override the build settings with `--vicat_build_cpus`,
`--vicat_build_memory` or `--vicat_taxonomy_memory`. These build resources
are separate from viCAT sample-analysis requests.

The full-input taxonomy replay used 113.8 GiB peak RAM and took about 2h 5m,
with all output values unchanged. The historical run used 310.2 GiB and took
23m 40s, but used different threads/runtime settings; this is not a controlled
performance comparison. Full-sized MMseqs stages historically exceeded 100 GiB.
Do not interpret the 96-GB DuckDB setting as a 96-GB machine requirement or
claim that the entire build fits on a 128-GB server. Complete fresh-build
validation of the lower-memory preset remains pending. Allow substantial
scratch space for DuckDB spill and the temporary uncompressed membership TSV.

For SLURM, use `-profile slurm` with your site's settings. Combine it with
`vicat_build_high_power` for the larger task request. The local scheduling
budget applies to local execution.

## viHARMONY workers

The default is one worker with a 32-GB task request. To use eight workers:

```text
--harmonizer_cpus 8
```

The parallel and serial implementations use the same decision rules. More
workers can require more RAM; `harmonizer_memory` is not multiplied automatically.
For example, when available, request `--harmonizer_memory '48 GB'` within a
`--max_memory '64 GB'` or larger local budget. These examples do not guarantee
that every dataset fits.

viHARMONY uses compact evidence storage, but its memory demand still grows with
the inputs and accumulated outputs. One saved RNA dataset with 1,759,634 final
records and eight workers used 68.7 GiB peak parent RAM, compared with 104.7 GiB
before this storage change; elapsed time was similar. These are measurements on
one dataset, not a memory bound or a general performance guarantee. Worker and
other process memory require additional headroom. For a run of that scale, a
128-GB task request and scheduling budget are a reasonable starting point if
the machine can provide them:

```text
--harmonizer_cpus 8 --harmonizer_memory '128 GB' --max_memory '128 GB'
```

The default remains a 32-GB request, not a promise that every run fits in 32 GB.
Monitor your own run and adjust requests, particularly on SLURM where exceeding
the allocated memory may terminate a task.

## Wall-time limits

viSUM does not impose task wall-time limits by default, including setup.
Queues may impose a maximum or require an explicit request. To request a limit
for an individual tool, use its `--*_time` setting, such as
`--harmonizer_time '48h'`. Tool-specific CPU, memory and time names are listed in
`visum.config` and in command-line help.

## SLURM

Use the `slurm` profile and a site configuration with your actual partition and
account. The shipped default queue name is a placeholder, not a portable choice.
For example, edit this `site.config` for your cluster:

```groovy
process {
  queue = 'YOUR_PARTITION'
  clusterOptions = '--account=YOUR_ACCOUNT'
  time = '48h'  // only if appropriate/required by your queue
}
```

Then run the same workflow and inputs with:

```bash
nextflow run /shared/viSUM/visum_nextflow.nf -profile slurm -c /shared/site.config \
  --input /shared/assembly.fasta --prefix sample01 --type dna \
  --outdir /shared/results/sample01
```

Use paths accessible from all compute nodes for the code, inputs, databases,
environments, work and results. Nextflow submits per-task resource requests;
SLURM controls allocation. Local `executor` budgets are not a cluster-wide RAM
reservation. Site-specific SLURM execution must be verified on your cluster.
Do not run the workflow as a large computation on a login node against local
execution defaults. Database setup can also use SLURM with suitable queue resources.

## Site overrides

Prefer `-c site.config` to editing workflow source. For example, share environments
with `conda.cacheDir = '/shared/visum/conda'`. Use the same site configuration
for setup and analysis. Keep it with your run's command and reference identities.
