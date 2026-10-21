# viCAT viral and nonviral databases

[Database setup](database-setup.md) · [Resources](resources.md)

viCAT requires both a viral reference and a class-aware nonviral reference.
Viral homology alone does not replace comparison against cellular, organellar
and plasmid proteins. viCAT is optional; use `--run_vicat true` in both setup
and analysis when you want its evidence.

## Automatic preparation

From the installation directory, after checking available resources:

```bash
nextflow run visum_nextflow.nf --setup --run_vicat true -profile vicat_build
```

Default destinations are:

```text
databases/vicat/metavr5/
databases/vicat/nonviral_classaware_v1/
```

Setup downloads the approved sources when the destinations are absent, builds
the two references sequentially, and validates them before use. Valid existing
databases are reused without reconstruction. No prebuilt public viCAT database
download is configured in this version.

### Viral sources

Setup automatically acquires the authors' `METAVR.faa.bgz` protein file and
`METAVR_main_table.parquet` metadata from the MetaVR portal. The shipped manifest
pins their sizes and SHA-256 checksums; users do not need the old browser-download
filenames. Allow space for approximately 50 GB of compressed protein sources,
1.09 GB of Parquet metadata, converted metadata, build scratch and the finished
reference. Those source sizes are not estimates of the final database size.

Interrupted transfers retain their saved prefix and validate HTTP identity
before resuming. Completed sources are checksum-verified for reuse. Unknown or
changed bytes are rejected, not silently substituted.

The metadata is converted to the existing compact format and a pinned
compatibility patch restores 1,409 vOTU assignments from the original reference
snapshot. It does not change taxonomy/method fields or declare the portal's
newer grouping incorrect. This preserves the selected reference identity.
Sources/receipts remain under `metavr5.sources/` beside the database; completed
database provenance records the inputs and patch identity.

### Nonviral sources

Setup uses the packaged list of **480 versioned RefSeq assemblies** to acquire
proteins and feature tables automatically. A separate NCBI Datasets installation
is not required. Downloads validate upstream checksums and record the exact
acquired bytes. An unavailable assembly causes an error rather than being
silently skipped or replaced by a later assembly version.

References are classified into cellular chromosome, unplaced cellular, plasmid,
plastid, mitochondrial and shared-nonviral groups. Each group is clustered
separately at 90% protein identity and 80% bidirectional coverage. The builder
requires all six groups to be represented. Source/receipt files remain under
`<nonviral-destination>.build/sources`, and finished provenance records the
selected assemblies and source fingerprints.

A later download of the same assembly accession need not have identical
annotation bytes. Record the finished database identity when reporting results.

## Resources

The viral build requests **32 CPUs and 350 GB RAM**, capped by `max_cpus` for
CPU allocation. The `vicat_build` profile provides a 32-CPU/400-GB local budget.
The nonviral build requests 16 CPUs/64 GB by default. These are configurable
allocations, not measured minimums or guarantees; lowering them is not known
to make the full build fit on a smaller machine. Verify RAM and disk before
starting. viSUM sets no build wall-time limit; queue policies still apply.

Use of an already-built database does not require the full build allocation.
Analysis has its own resource settings. See [resources](resources.md).

## Use existing databases

If both databases are already installed elsewhere, supply their paths during
setup and analysis:

```text
--run_vicat true
--vicat_db /shared/vicat/viral
--vicat_nonviral_db /shared/vicat/nonviral
```

Prepared references are validated, not modified. A valid viral reference
contains `IMGVR5_UViG_representatives.dmnd` and
`IMGVR5_UViG.vicat_taxonomy_lookup.parquet`; the older development layout with
`diamond/` and `taxonomy/` subdirectories is also accepted. Do not rename files
or replace their provenance to make a different reference appear compatible.

## Supply original local sources instead of downloading

The manual viral route accepts the original checksum-pinned protein and compact
metadata pair. It does not accept the portal Parquet file as the compact input:

```text
--vicat_metavr_proteins /inputs/IMGVR5_UViG.faa.gz
--vicat_metavr_metadata /inputs/IMGVR5_UViG.tsv.gz
```

Do not also set `--vicat_db`. Set `--vicat_managed_db` only to override the build
destination. Alternative source releases need deliberate identity changes and
validation; do not change checksum pins just to bypass an error.

For local nonviral NCBI sources, supply these three options:

```text
--vicat_nonviral_manifest /inputs/assemblies.tsv
--vicat_nonviral_package_root /inputs/protein_package
--vicat_nonviral_metadata_root /inputs/feature_tables
```

The tab-separated manifest requires `assembly_accession` and `cellular_group`.
Expected local files are:

```text
<package_root>/ncbi_dataset/data/<assembly_accession>/protein.faa
<metadata_root>/<cellular_group>/<assembly_accession>/<assembly_accession>.feature_table.txt.gz
```

Alternatively, set `--vicat_nonviral_classified_dir` to a directory containing
the six class-named `.faa.gz` files, `vicat_nonviral_reference_metadata.tsv.gz`
and `vicat_nonviral_classification_summary.tsv` in the classifier's format.
Do not combine the raw-source and preclassified routes.

`--vicat_metavr_auto_download false` and `--vicat_nonviral_auto_download false`
disable their respective automatic source routes. Existing completed databases
remain usable; missing source files will not be fetched by those disabled routes.

## Interrupted builds and updates

Rerun setup with the same sources, destination and settings. Managed installation
locks prevent overlapping builders, and matching checkpoints may be reused.
Do not remove installation markers or incomplete data just to bypass validation.
Keep failed logs and build files for diagnosis.

Use new destination paths to build changed references. A source option is not
an instruction to overwrite a completed database. Keep versioned databases for
analyses that must remain reproducible. See [troubleshooting](troubleshooting.md).

## Attribution

Cite MetaVR and the selected nonviral sources in addition to viSUM and the
underlying tools. Public source access is not a blanket redistribution license;
verify applicable terms before hosting a derived reference bundle. This software
archive does not redistribute the large databases.
