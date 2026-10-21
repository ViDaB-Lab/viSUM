# HISTORICAL PREPARATION PLAN — not the current RC launcher

# Preparing the geNomad-author dataset for viSUM

Preparation scripts do not launch analysis or change retention thresholds.
Python 3.10+ and the standard library are sufficient. The separate launch script
uses the existing server installation and an isolated, checksummed pipeline copy.

## Approved final design (2026-09-29)

Keep all 1,445,955 original fragments. Reviewed RNA source assignments run as RNA;
genomic chromosome/plasmid controls run as DNA. Unresolved viral parents, including
conflicting records and DNA annotations lacking independent genome-type support,
run **both** routes under explicit `--unresolved-policy both` authorization.
The finalized map has 78,169 RNA parents (89,225 fragments). The remaining 162,242
viral parents (514,339 fragments) receive paired executions. This produces
1,960,294 executions, not 1,960,294 independent biological observations.

Candidate RNA records that are merely mRNA or lack genome-type support also remain
unresolved unless they match the published RNA reference exactly. No benchmark
predictions or outcomes influence routing. `finalize_genomad_routing.py` records
every parent-level decision and the source checksums.

The three evaluation strata are `source_routed`, `unresolved_dna`, and
`unresolved_rna`. Never pool the paired strata, select the better result, or score
their union as a single source-routed sensitivity. Keep the same original IDs for
pairing; the unique execution key is `(benchmark_batch, sequence_id)`.

## Server preparation command

Transfer `genomad_final_bundle.tar.gz` to the benchmark root, then run this child
shell. It never launches discovery tools and will not close the parent terminal
if a check fails. Extraction goes into a new temporary directory.

```bash
bash <<'BASH'
set -euo pipefail
BENCH=/home/alex/programs/visumproject/benchmark_data/genomad_paper_v1
PACKAGE="$BENCH/genomad_final_bundle.tar.gz"
test -s "$PACKAGE"
DEST=$(mktemp -d "$BENCH/final_bundle_XXXXXXXX")
tar -xzf "$PACKAGE" -C "$DEST"
bash "$DEST/genomad_final_bundle/pipeline/bin/prepare_genomad_final.sh" \
  "$BENCH" "$BENCH/alpha_inputs_dual_v1"
printf '\nBundle for the later launch: %s\n' "$DEST/genomad_final_bundle"
BASH
```

The archive freezes the committed pipeline plus the explicitly listed composition
and preparation changes. It excludes unrelated, uncommitted viHARMONY experiments.
`bundle_provenance.json` and `bundle_checksums.sha256` identify the actual files.
No Git pull, checkout, overwrite of the server checkout, or push is needed.

After successful preparation, the **separate explicit launch** is:

```bash
bash /actual/extracted/genomad_final_bundle/pipeline/bin/launch_genomad_final.sh \
  /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/alpha_inputs_dual_v1 \
  /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/alpha_run_dual_v1
```

Replace only the extracted bundle path with the path printed by preparation.
The launcher defaults to `/home/alex/programs/visumproject/nextflow`. It reuses the
original project's databases, tools and Conda cache, checks the known complex
database directories, records its exact command, and refuses an existing run
directory. Environment overrides are `NEXTFLOW_BIN`, `VISUM_ORIGINAL_PROJECT`,
`VICAT_DB`, `NONVIRAL_DB`, `VITAP_DB`, `VCONTACT3_DB`, and `PYTHON_BIN`.
Directory existence does not establish database completeness; pipeline preparation
tasks still perform their own checks. Runtime checks require the Linux server.

## What the supplied metadata establishes

The uploaded `metadata.tar.gz` contains the original author instructions and a
prepared fragment truth table. Its SHA-256 is
`efe8d774f33758c01a1b3693e6d88df8f895651c11b2309dd041c9a62e39f7be`.

The table contains 1,445,955 fragments: 733,534 chromosome, 108,857 plasmid and
603,564 virus fragments, totaling 9,367,521,658 bases. The viral fragments derive
from 240,411 parent IDs. Many are metagenomic identifiers, not GenBank accessions.
Neither the truth table nor the four-column reference-cluster table specifies
viral DNA/RNA molecule type. A nucleotide FASTA alphabet cannot distinguish a
DNA genome from an RNA genome represented with T instead of U.

The six `author_weight_split_*` fields are described by the authors as **training
weights**. We preserve their exact strings; we do not designate them as evaluation
weights or infer held-out membership from zero weight. The reference-cluster
table explicitly supplies the author split assignment, which is joined separately.
Running current pretrained software on this collection is not an independent
training-unseen benchmark or an exact cross-validation reproduction.

Source: supplied `source/benchmark_data/README.md` and
`source/reference_sequences/README.md`, from https://doi.org/10.5281/zenodo.8049246.

## 1. Audit the metadata

From the viSUM repository on the server:

```bash
python bin/prepare_genomad_benchmark.py inspect \
  --metadata /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1 \
  --out /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/alpha_metadata_audit
```

Alternatively, `--metadata` accepts the uploaded metadata tarball. Only the four
named metadata files are read; archive paths are not extracted. Use a **new output
directory** for each attempt. Existing directories are never overwritten.

Outputs include a disk-backed audit index, reconciled counts, checksums,
`routing_required.tsv`, and `routing_used.tsv`. An inspection can succeed while
reporting `METADATA_VALID_ROUTING_INCOMPLETE`; that is **not** permission to launch.
Malformed rows, duplicated IDs, unexpected classes, nonfinite/negative weights,
retained author exclusions, absent cluster members and inconsistent counts fail.

## 2. Supply source-supported viral routing

### Recover sequence identity when published names differ

The uploaded metadata does not contain original parent sequences. The published
RVMT RNA metadata uses identifiers that did not directly match the benchmark's
parent IDs in the local audit. Do not infer molecule type from that failed join.
The original author ZIP contains the parent reference FASTA, which can be read
without extracting its 5.4-GB compressed member:

```bash
python bin/fingerprint_genomad_parents.py \
  --truth /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/prepared/truth_metadata.tsv \
  --archive /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/downloads/genomad_supplementary_data_code.zip \
  --out /path/to/new_parent_fingerprint_directory
```

This exports headers, lengths, original sequence hashes and strand-independent
hashes for exact-content joins to source collections. Canonical hashes normalize
case and U/T only. They do not establish identity for trimmed, mutated or partial
sequences. The script checks all required parents are present exactly once,
consumes the compressed member to EOF for archive integrity checks, and packages
the small result as `routing_evidence.tar.gz`. A missing-parent result has status
`INCOMPLETE` and exit code 2. Neither a successful fingerprint export nor a source
nonmatch assigns a DNA/RNA route or authorizes launch.

Published RNA source: https://zenodo.org/records/7368133 and
https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/.

Create a tab-separated map with these columns:

```text
parent_sequence_id  input_type  evidence_source  evidence_record
```

The spaces above illustrate columns; the actual separator must be a TAB. Values
for `input_type` must be `dna` or `rna`. Each exact viral parent appears at most
once and needs a nonempty source citation and record identifier. The pipeline's
own predictions are not acceptable independent routing evidence. Do not default
unmapped viral sequences to DNA, infer molecule type from identifier prefixes,
or propagate one member's route across an entire reference cluster.

`build_genomad_routing.py` can join saved NCBI GenBank/INSDSeq XML to exact parent
accessions using accession-version and molecule-type fields:

```bash
python bin/build_genomad_routing.py \
  --required /path/to/alpha_metadata_audit/routing_required.tsv \
  --genbank-xml /path/to/saved_records.xml \
  --source 'NCBI nuccore; saved metadata snapshot DATE' \
  --out /path/to/new_routing_audit
```

This is an offline optional helper, not a complete solution for metagenomic IDs.
It records input hashes and emits `routing.tsv` plus the remaining unresolved
IDs. Versioned parents require an exact accession-version match. Versionless
parents can match the explicit primary accession, with the returned version
recorded as provenance. Conflicting metadata fails. Source-associated metadata
is still needed for non-accession metagenomic sequences when a single route is
required. Under the approved dual policy, those records remain explicitly
unresolved and get both executions. The XML helper itself makes no network calls.

**Record molecule type is not always viral genome type.** A GenBank DNA record
can represent an RNA-virus cDNA clone or integrated copy. Treat the XML helper's
map as candidate evidence requiring review, not an automatically approved route.
An exact match to a published RNA-virus reference can conflict with that field;
preserve both observations and resolve the biological interpretation explicitly.

`audit_genomad_routing_evidence.py` performs an offline reconciliation of the
fingerprint package, saved exact RNA-reference matches, and cached NCBI JSON
summaries. It verifies fingerprint integrity, requires exact accession-version
and length matches, retains every evidence source, and separates conflicts and
missing records. Its output is deliberately named `routing.candidate.tsv` and its
summary always states `launch_ready: false`; it does not adjudicate cDNA records
or approve a benchmark launch. The supplementary ten tests cover exact joins,
version/length mismatches, invalid RNA identities, missing routes and conflicts.

```bash
python bin/audit_genomad_routing_evidence.py \
  --fingerprints /path/to/routing_evidence.tar.gz \
  --rna-matches /path/to/rna_exact_matches.json \
  --ncbi-cache /path/to/ncbi_summary_cache \
  --out /path/to/new_routing_evidence_audit
```

The chromosome and plasmid source classes route to DNA based on their benchmark
source class. Those are genomic DNA fragments, not an RNA-negative transcript panel.

## 3. Build and reconcile the FASTAs

Use the reviewed, potentially partial routing map with the approved dual policy:

```bash
python bin/prepare_genomad_benchmark.py build \
  --metadata /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1 \
  --routing /path/to/routing.reviewed.tsv --unresolved-policy both \
  --fasta /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/prepared/genomad_paper_fragments.fna.gz \
  --shards 32 --seed visum-alpha-gnp-v1 \
  --out /home/alex/programs/visumproject/benchmark_data/genomad_paper_v1/alpha_inputs_v1
```

This command streams the FASTA and requires exactly one matching record for every
truth ID, with the expected length and a valid DNA/IUPAC representation. Original
IDs and sequence case/content are preserved. All fragments of a parent share a
deterministic SHA-256 hash shard. Shards are grouped by DNA/RNA input type, **not**
positive/negative labels; fragment order follows the source FASTA. `--shards` is
the number of possible shards per type, not a strict sequence/base limit or a
cross-validation partition. Different parents in a reference cluster can share
or differ in execution shards without changing their recorded author split.

Without `--unresolved-policy both`, the default still refuses incomplete routing.
Missing, duplicated,
unexpected or length-mismatched FASTA records prevent manifest creation. Partial
failed outputs are kept for diagnosis and carry a `FAILED` summary; they must not
be used. Resolve the issue and choose a new output directory.

Successful output includes:

- `manifest.csv`: viSUM `prefix,type,fasta` header with relative FASTA paths;
  supply its directory as `--indir` when launching later.
- `truth.routed.tsv`: all original columns plus parent, author cluster split,
  source-supported input type (or `unresolved`), batch(s), and sequence SHA-256;
  exactly one row per original fragment.
- `execution_map.tsv`: each actual batch/sequence execution, with its route,
  evaluation stratum and hash. An unresolved fragment has exactly two rows.
- `batch_composition.tsv`: counts and bases by source class within each batch.
- `preparation_summary.json`, `routing_used.tsv`, and `checksums.sha256`.

## Launch boundary and evaluation

`INPUTS_PREPARED_NOT_LAUNCHED` means sequence preparation passed, not that Nextflow
or tool installation was tested. Full FASTA reconciliation must run on the server;
the local metadata package contains no complete sequence FASTA.

The wrapper now exposes validated `--genomad_composition auto|metagenome|virome`.
The normal pipeline default remains `auto`. This benchmark explicitly uses
`--genomad_score_calibration true --genomad_composition metagenome`, default
geNomad filters, and `--genomad_cleanup false` to retain native evidence.
The fixed preset avoids composition estimates derived from positive-only route
batches. It does not claim geNomad's estimated FDR equals this benchmark's empirical
FPR. Supported syntax was checked against geNomad 1.12.0's official CLI source:
https://github.com/apcamargo/genomad/blob/v1.12.0/genomad/cli.py.

`verify_genomad_preparation.py INPUT_DIRECTORY` independently checks checksums,
unique truth IDs, both paired arms, execution-map identity, manifest coverage,
FASTA lengths and sequence hashes. It refuses altered or incomplete inputs.

Freeze the actual pipeline commit and configuration; preserve database provenance.
Review batch sizes and resource requirements before scheduling this 9.37-Gb input.
Retain native per-tool outputs for comparisons, rather than assuming adapters
preserve every native prediction. Keep original-label FPR, parent-level detection,
and intact-input-fragment recovery separate. These are metagenome-sized fragments,
not uniformly complete genomes. Mark previously examined pilot fragments during
evaluation; this utility neither excludes nor relabels them. RNA-negative controls
remain a separate benchmark, not supplied by these genomic chromosome/plasmid rows.

## Tests

```bash
python -m unittest discover -s tests -p test_prepare_genomad_benchmark.py -v
```

Only `launch_genomad_final.sh` launches analysis. Preparation, routing review,
verification and packaging scripts never do. No script commits or pushes changes.
