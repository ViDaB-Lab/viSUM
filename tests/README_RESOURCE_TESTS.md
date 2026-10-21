# Resource-default verification v2 — unfinished checks only

This package is a frozen copy of the isolated alpha candidate, not a production
update. Run only `tests/run_resource_defaults_linux.sh`. It never executes the
scientific workflow, any analysis binary, a database installer, or a resume.
No completed smoke/student cache is needed or read.

The runner reads the previously used Python 3.11 environment and copies the
installed Nextflow 26.04.6 launcher/framework into a fresh private directory.
It installs nothing, downloads nothing, runs sequential Nextflow invocations,
and retains logs/observations on success or failure. Strict parser v2 and
offline mode are explicit. No SLURM jobs are submitted.

## What is checked

1. Verify all package file hashes. The twelve resource contracts and six
   successful configuration checks from v1 are not rerun on the server.
2. Resolve shipped configuration with native `nextflow config -o json` for
   **docker and singularity only**. Default, local, local_safe, local_full,
   vicat_build and slurm already passed in the supplied v1 report.
   All twenty time defaults must be null; viHARMONY defaults must
   be 32 GB/one worker, and profile memory budgets must remain correct. This
   step only reads configuration, including for scheduler/container profiles.
3. Generate 27 tiny local resource probes by copying actual module time
   expressions, plus an inherited utility-time probe. Analysis/install scripts
   are replaced entirely by `printf` of resolved task resource settings.
4. Run defaults: every resolved task time must be null; the viHARMONY probe
   must request 32 GB and one CPU.
5. Run typed parameter overrides: all twenty time parameters are set to three
   hours, harmonizer memory to 48 GB, requested workers to eight with a two-CPU
   budget. The probes must resolve to those times/memory and two capped workers.
6. Run site selectors: utility gets five hours; harmonizer gets seven hours,
   48 GB and two CPUs; other probes remain unlimited. All traces must complete.

Each invocation uses at most two reserved task CPU slots; the JVM has two active
processors and a 1 GiB maximum heap. The 32/48 GB values are scheduling requests
on trivial echo-like tasks, **not allocated test data or a RAM stress test**.
Other probes request 64 MB because their memory sizing is outside this check.
These reservations/JVM settings are not hard OS limits on the full process tree.

The three-minute fixture watchdog is deliberate and does not change shipped
pipeline time defaults. On timeout it terminates only its freshly spawned test
process group, never processes found by a username/name search.

At completion the runner prints the absolute path to
`resource_defaults_review_20261015_v2.tar.gz`. Upload that file. It contains the
snapshot identity, report, configuration resolutions, probe graph, observations,
trace and logs, but excludes copied runtime jars and all biological data.

## Scope and current verification

The v1 Docker check emitted `WARN: Unrecognized config option 'docker.pull'`
before valid JSON. Its harness then failed by decoding the entire console as
JSON; no resource probes started. The isolated Docker profile now removes that
ignored option without changing its enabled state or user mapping.

The v2 parser accepts a single JSON object with explicit `WARN:` lines before
or after it, records every warning in the report and `config.warnings.json`,
and retains original console logs. It rejects error/noise prefixes, malformed
JSON and additional documents rather than blindly searching for a later `{`.
An unsupported `docker.pull` key is still rejected by configuration validation.
An exact 9,088-byte copy of the uploaded Docker console output is preserved in
`fixtures/resource_defaults/docker_warning_20261015.json`, with its source hash.
Regression tests reproduce the old failure and verify the corrected parser.

Portable generator/validator tests include negative controls for missing time
parameters, restored time limits, incorrect memory/workers, missing observations
and corrupted resource values. They do not substitute for Linux execution.
This Windows development host has no Nextflow/Java runtime; native execution of
this new resource harness must be verified on Linux. Existing full-workflow and
cache tests are separate and are intentionally not rerun by this runner.

Passing these probes validates resource resolution, not memory adequacy for
every input, scientific equivalence, scheduler-site acceptance, or fresh tool/DB
installation. Decision code remains unchanged from the validated candidate.
