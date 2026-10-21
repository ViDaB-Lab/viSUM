# Native cache-reuse regression tests

These tests use real Nextflow and unchanged copies of six viSUM process modules.
geNomad, VirSorter2, DIAMOND and viHARMONY executables are tiny substitutes;
normalization runs normally. No databases are downloaded or biological analyses
performed. Each actual substitute execution writes an independent event file.
Assertions check both those events and Nextflow trace statuses.

## Running

Use a separate Linux/Unix development checkout with Python 3.11+ and a working
Nextflow/Java installation. The Nextflow runtime must already be available: the
fixture sets offline mode and disables Conda and containers. Use the same
Nextflow version as the supported installation when establishing the baseline.
Tests reserve two local task CPU slots in total. The Nextflow JVM has a 1 GiB
maximum Java heap and processor count of two; these are not a hard quota on
total process-tree memory or CPU usage.
The fixture supplies only its copied viSUM config and test overrides using
Nextflow's [exclusive configuration option](https://docs.seqera.io/nextflow/cli#hard-configuration-override),
so an installation's default config cannot select a different executor.

From that isolated checkout:

```bash
python3 -B tests/test_cache_reuse_fixture.py
VISUM_CACHE_TEST_ARTIFACTS="$(mktemp -d /tmp/visum-cache-tests.XXXXXXXX)" \
  python3 -B tests/test_cache_reuse_nextflow_integration.py
```

The direct native runner exits unsuccessfully if its runtime is unavailable.
Ordinary test discovery skips the five native tests on unsupported systems;
skips are **not** evidence of cache reuse. The artifact directory retains the
per-invocation console/log/trace files, scratch modules, work/cache and execution
events. Without that environment variable, temporary fixtures are removed.

## Contracts

1. DNA and RNA samples complete, then every task is cached on resume. Calls from
   different working directories use the fixture's fixed launch directory. A
   read-only asset check writes a separate report without changing search inputs.
2. Interrupt the disposable Nextflow process while both viHARMONY fixtures wait;
   resume finishes them without reexecuting any of the four completed searches.
3. Change only viHARMONY workers from one to two: completed work stays cached
   under Nextflow 26.04.6. The interrupted-task test separately verifies that new
   viHARMONY executions receive two workers, without repeating searches.
4. Remove one fixture-owned cached search output: the affected search must rerun.
   This negative control checks that the test detects execution, not just text.
5. Keep the production prepare declarations, including `cache false`, but replace
   installers in scratch copies with tiny symlink/metadata producers. Repeating
   this inline preparation recreates task-local database paths and must expose reexecution of viCAT's viral and nonviral
   searches. This characterizes the existing cache hazard; it is not the desired
   behavior for alpha. The stable-path branch is the contrasting target contract.

## Scope and limits

This is a focused module-level graph, **not** an end-to-end test of
`visum_nextflow.nf`, the current launcher, real database validators, or scientific
decisions. It deliberately does not fix those components or manipulate cache
records. Its viHARMONY files contain test sentinels, not scientific results;
identical hashes here do not establish biological equivalence.

The separate full-workflow tests below cover the actual entry point,
missing-asset/input preflight, zero-call routing and setup idempotence.
Keep viHARMONY decision-equivalence tests as a separate requirement.

## Actual entry-point tests after the native refactor

`test_native_workflow_integration.py` now exercises the actual
`visum_nextflow.nf`, shipped automatic `nextflow.config`, read-only preflight,
standardizers, discovery gate, refinement and real viHARMONY on tiny DNA/RNA
zero-call inputs. Only geNomad's external binary and the setup MMseqs binary are
substituted. Two separate setup invocations reuse the prepared fixture database;
the original analysis resumes with every scientific task cached. Missing assets
and duplicate prefixes must fail before any scientific task starts.

```bash
python3 -B tests/test_prepared_assets.py
python3 -B tests/test_native_workflow_contracts.py
VISUM_CACHE_TEST_ARTIFACTS="$(mktemp -d /tmp/visum-native-tests.XXXXXXXX)" \
  python3 -B tests/test_native_workflow_integration.py
```

This is still not an all-tools biological benchmark or a clean installation of
real databases. The native-refactor tests need Linux execution before deployment.

Never point the fixture at student run/cache/output directories. It accepts no
server, database, input or work-directory override; all assets and cache state are
generated beneath a fresh temporary directory.

## Initial verification — 2026-10-13

On the Windows development host, all 11 portable fixture checks passed. Combined
with the selected viHARMONY, setup and resume-contract regression tests, 123 tests
passed and six Linux/Nextflow-only tests were skipped (five new cache tests and
one existing setup test). Python syntax checks also passed. The viHARMONY and
evidence-schema source hashes stayed unchanged from the validated multiprocessing
build at the start of this work.

**Native Nextflow execution remains unverified.** No local WSL distribution or
Nextflow runtime was available. Establish the Linux baseline before treating
these cache contracts as passed or changing production setup/run/resume behavior.
No student installation, analysis, database, results or cache was accessed.

## First Linux run and fixture correction — 2026-10-13

The supplied server report used Nextflow 26.04.6: portable checks passed (10
passed, one expected runtime-availability skip); native tests ran in 46 seconds
with one pass and four failures. The inline-preparation characterization passed.

The other tests accidentally selected that same preparation branch: the harness
passed `--fixture_inline_prepare false`, but the v2 parser preserved this as the
nonempty string `"false"`, which the fixture's `if` treated as true. The cold
stable-path test scheduled 15 tasks instead of the expected 12, including all
three preparation tasks. Resume traces and execution events confirmed fresh
preparations followed by repeated geNomad and viCAT searches. They also showed
normalization and VirSorter2 reusing their caches. This is not a valid result
for the intended stable-path branch.

The corrected harness uses a JSON `-params-file` with a real Boolean and numeric
CPU values. The fixture rejects non-Boolean flags before scheduling; the native
runner also checks that stable runs have no preparation tasks and inline runs
have exactly the three expected preparations. Portable tests cover typed JSON
round trips and rejection of string Boolean values. This follows the documented
[Nextflow 26.04 parameter typing behavior](https://training.nextflow.io/latest/nfcore_use/02_configure_execution/#12-set-parameter-values).

That fixture correction did not change production code. The subsequent Linux
report passed completed resume, interrupted resume, missing-output detection and
the inline-preparation control. Its worker-change expectation was incorrect:
changing only `task.cpus` does not invalidate completed Nextflow tasks. The test
now expects completed viHARMONY work to stay cached and checks that an unfinished
task receives the new worker allocation. The revised tests and new native
full-workflow tests still require Linux validation.
