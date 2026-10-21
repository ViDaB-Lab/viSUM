# Troubleshooting

[Usage](usage.md) · [Setup](database-setup.md) · [Resources](resources.md)

## Find the relevant error

Read the Nextflow error and the failed task's work directory. That directory
contains `.command.sh`, `.command.out`, `.command.err` and `.command.log` where
available. The launch directory has `.nextflow.log`; `<outdir>/reports/trace.txt`
records task status. Preserve these files before changing or cleaning anything.

## Missing or incomplete databases

Analysis never installs missing assets. Run setup with the same tool selection,
database paths and site config as analysis. Verify paths from the execution host,
not just the machine where you launched Nextflow. TEsorter's lack of a standalone
setup task is expected; a nonzero exit or failed runtime import is not.

Do not remove `.visum_installing` markers or edit receipts/checksums to force an
incomplete installation to pass. Rerun matching setup or inspect the reported
failure. Use a new destination for genuinely changed source/reference identities.

## Task memory exceeds available resources

A scheduling budget cannot create RAM. viCAT's full viral reference build is
configured to request 350 GB, whereas ordinary analysis uses separate requests.
Check the selected profile, per-tool memory and actual machine/queue capacity.
Do not merely lower a setting and assume the same full build will fit.

## Resume does not reuse expected tasks

Keep `.nextflow/` from the original launch directory and the original work path.
Repeat the same analysis command and select the correct run with `-resume RUN_NAME`
if needed. Check `nextflow log` to distinguish setup from analysis sessions.
Changed inputs, source, references or task parameters can invalidate a cache.
Do not manually edit cached outputs or copy another run's database/cache into it.

## Outputs look stale after a changed run

Use the harmonizer manifest as the current inventory. Current emitted files
replace same-named published copies, but publishing does not delete obsolete
files from an older layout or optional view. Use a new output directory when
changing versions or when you want to preserve an older report.

## TEsorter reports `AF_UNIX path too long`

This version gives TEsorter a short `/tmp` temporary root for multiprocessing.
Ensure `/tmp` is writable on the execution node. Check the task's generated
script and software version if an older environment/config override is in use.
Do not replace code or environments underneath an active analysis.

## vCONTACT3 reports a missing Parquet engine

The supplied environment includes the expected pandas/pyarrow combination and
an in-memory round-trip check. Confirm Nextflow is using this version's YAML and
the intended Conda cache. Do not mutate a shared environment while another run
uses it. Preserve the error and environment identity before repairing a private
environment or choosing a new cache.

## No viral sequences were exported

Empty FASTAs and header-only tables are valid zero-call outputs. Verify the
workflow completed, then inspect the summaries, dispositions, review tables and
enabled programs. An empty export is not proof that the input contains no viruses.

## Report an issue

Provide the viSUM version/commit, Linux/runtime versions, input type, redacted
command, failed task name and relevant logs. Explain expected versus observed
behavior. Do not upload credentials, confidential sequences or private sample
identifiers; use a small non-sensitive example when possible.
