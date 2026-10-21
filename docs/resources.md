# CPU, memory and execution settings

[Usage](usage.md) · [Setup](database-setup.md)

## Local execution

The default local scheduling budget is 8 CPUs and 64 GB. Set `--max_cpus` and
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
| `local` / default | 8 | 64 GB | Default local settings |
| `local_safe` | 8 | 64 GB | Explicit conservative local preset |
| `local_full` | 32 | 128 GB | Larger dedicated local budget |
| `vicat_build` | 32 | 400 GB | Large viCAT source construction |

CLI overrides still apply. A profile is a configuration preset, not evidence
that the machine has those resources. The viral viCAT build currently requests
350 GB; this is a configurable allocation, not an established minimum. Its
reference-building allocation is separate from sample-analysis requests.

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
