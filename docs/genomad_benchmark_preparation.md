# Prepared geNomad-derived benchmark: reuse, do not restart

The approved 29 September preparation design is unchanged: 1,445,955 unique input
fragments; 1,960,294 planned route executions; 128 planned batches. Source-supported
RNA parents run as RNA, genomic chromosome/plasmid controls as DNA, and unresolved
viral fragments retain paired DNA/RNA executions. No predictions assign routing.

The server preparation directory recorded in the plan is:
`/home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/alpha_inputs_dual_v1`.
Its current presence/contents have NOT been reverified from this local machine.
Before a benchmark, run the read-only verifier from the extracted RC:

```bash
python3 -B /path/to/candidate/bin/verify_genomad_preparation.py \
  /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/alpha_inputs_dual_v1
```

Keep `truth.routed.tsv`, `execution_map.tsv`, `manifest.csv`,
`preparation_summary.json`, reviewed routing, author cluster/split metadata and
all checksums. Do not regenerate inputs or choose routes based on new outcomes.
Fragment execution strata must remain separate; see [evaluation contract](benchmark-evaluation.md).

The standalone `prepare_genomad_benchmark.py`, `finalize_genomad_routing.py`,
`fingerprint_genomad_parents.py`, `build_genomad_routing.py`,
`audit_genomad_routing_evidence.py` and verifier are carried forward with identical
Python ASTs to the approved preparation snapshot (line endings may differ).
They are available if a genuine preparation
repair is required; this consolidation neither rebuilds nor launches anything.

The old `genomad_final_bundle` and its launch/preparation wrappers are historical
snapshots, NOT the alpha candidate. Do not reuse their pipeline code, installed
environment-cache assumptions, or 24-hour viCAT limits. The preserved historical
design/commands are in `history/genomad-preparation-20260929.md`, not current launch
instructions. Use ordinary Nextflow on this candidate with separate setup and
explicit prepared-input/database/resource parameters, after the fresh-install gates.
The scientific benchmark preset is frozen in `release/benchmark-scientific.params.json`.
Native baselines must record actual versions/configs; this is a geNomad-derived
evaluation, not an automatic reproduction of the paper's training/cross-validation.
