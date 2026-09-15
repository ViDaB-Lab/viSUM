# viSUM documentation

The [main README](../README.md) introduces the four modules, programs and preliminary performance. Detailed operation and evidence interpretation live here.

## Using viSUM

- [Usage guide](usage.md): installation, single and batch inputs, resources, database locations, output files and troubleshooting.
- [Preparation-only setup](database-setup.md): prepare databases without analysis inputs, installation checks and smoke-test scope.
- [viCAT database preparation](vicat-database-preparation.md): viral and nonviral resources, local source builds and supplied databases.
- [VITAP database lifecycle](vitap-database-lifecycle.md): database selection, reuse and explicit updates.

VITAP and vConTACT3 are enabled by default for taxonomy refinement. viCAT remains opt-in because it requires separately supplied or built viral and nonviral resources.

## Benchmarks and interpretation

- [Preliminary all-tool comparison](benchmarks/preliminary-2026-09.md): positive sensitivity, nominal-negative retention and subgroup results.
- [Retained-negative context](benchmarks/retained-negative-context.md): viral homologs, coherent phage-like modules, mixed content and limits of source labels.
- [VirSorter2 evidence-policy follow-up](benchmarks/virsorter2-hallmark-followup.md): finalized review-only treatment of unscored native hallmark calls.
- [Experimental policy comparison](benchmarks/retained-negative-context.md#experimental-policy-comparison): saved-result tests of mixed-content and replication-only review; neither changes the default retention policy.

Benchmark documents are dated development records, not a substitute for run-specific configuration and provenance.
