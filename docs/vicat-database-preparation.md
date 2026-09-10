# viCAT database preparation

viCAT requires **two independent databases**: MetaVR-derived viral proteins and
taxonomy, plus class-aware nonviral references. Neither substitutes for the other.
Both prepare processes validate their outputs before downstream consumers use them.
Their heavyweight builds are serialized (viral, then nonviral).

Preparation is available through `--setup` without sample inputs, or within an
analysis run. Add the options below to either command. See [setup](database-setup.md).
No public prebuilt download URL is configured or implied.

## Reuse completed databases

```text
--run_vicat true
--vicat_db /shared/references/vicat/viral
--vicat_nonviral_db /shared/references/vicat/nonviral
```

Existing databases are validated without modification. Invalid existing targets
are rejected, not replaced. The legacy viral development layout remains accepted.
These paths also support future downloaded prebuilt bundles after extraction.

## Build both databases from local inputs

```text
--run_vicat true
--vicat_metavr_proteins /inputs/MetaVR/IMGVR5_UViG.faa.gz
--vicat_metavr_metadata /inputs/MetaVR/IMGVR5_UViG.tsv.gz
--vicat_managed_db /shared/references/vicat/viral
--vicat_nonviral_db /shared/references/vicat/nonviral
--vicat_nonviral_manifest /inputs/nonviral/assemblies.tsv
--vicat_nonviral_package_root /inputs/nonviral/protein_package
--vicat_nonviral_metadata_root /inputs/nonviral/feature_tables
```

Do not also set `--vicat_db` when supplying MetaVR sources. The viral checksums
and expected record counts in `visum.config` pin the tested source release; a
different release requires deliberately updating and validating those settings.

The nonviral source manifest is tab-separated and requires `assembly_accession`
and `cellular_group`. Optional organism/taxid fields are retained by the classifier.
Each assembly needs these exact local paths:

```text
<package_root>/ncbi_dataset/data/<assembly_accession>/protein.faa
<metadata_root>/<cellular_group>/<assembly_accession>/<assembly_accession>.feature_table.txt.gz
```

The workflow does not invent an assembly selection or download these NCBI inputs.
Use the frozen reference inventory for benchmark reproduction. Changing the
nonviral reference inventory changes the scientific resource and requires new
performance validation.

The classifier separates chromosomes, unplaced cellular sequences, plasmids,
plastids, mitochondria and shared nonviral proteins. Ambiguous/unsupported records
are audited as exclusions. The existing builder clusters each class separately
at 90% identity and 80% bidirectional coverage, then creates DIAMOND and metadata
artifacts. The current builder requires all six classes to contain sequences.

### Already-classified inputs

Instead of the three NCBI arguments, supply:

```text
--vicat_nonviral_classified_dir /inputs/nonviral/classified
```

This directory must contain `cellular_chromosome.faa.gz`,
`cellular_unplaced.faa.gz`, `plasmid.faa.gz`, `plastid.faa.gz`,
`mitochondrial.faa.gz`, `shared_nonviral.faa.gz`,
`vicat_nonviral_reference_metadata.tsv.gz`, and
`vicat_nonviral_classification_summary.tsv`, in the existing classifier's format.
The raw and classified input routes are mutually exclusive.

## Reuse, recovery and resources

- Builds occur only when their target database is absent. Source parameters are
  bootstrap inputs, not update commands. To rebuild with different sources, use
  new viral/nonviral destination paths; do not overwrite a frozen reference.
- Preparation validation is not skipped by Nextflow `-resume`. This can require
  substantial I/O for checksums and metadata counts, but does not rebuild a valid
  database. Other analysis task caching remains available.
- Managed builds use destination locks. User source files are not deleted.
- Viral checkpoint recovery checks pinned source hashes and a builder/parameter
  fingerprint. Old checkpoint directories lacking this stamp are rejected;
  inspect them and use a new destination/work path rather than deleting blindly.
- Nonviral raw classification records input/output fingerprints. Clustering
  retains its existing input/parameter/builder fingerprint checks. Interrupted
  work remains under `<nonviral-destination>.build`.
- Raw-source builds retain the assembly manifest and source fingerprint inventory
  in the nonviral package. Preclassified inputs should be accompanied by their
  own upstream provenance before redistribution.
- Viral defaults: 32 CPUs, 350 GB memory, 12 hours. Nonviral defaults: 16 CPUs,
  64 GB memory, 24 hours. These are configurable resource requests, not measured
  guarantees that every input collection will fit or finish. Check disk space,
  available resources and scheduler limits before building.
- No analysis decision thresholds are changed by this preparation wiring.

## Hosting policy

Keep the builder even if prebuilt bundles are hosted later. Publish versioned
packages with checksums, source accessions/releases, build provenance and notices.
MetaVR redistribution terms still need to be established for the viral bundle.

For NCBI-derived nonviral references, NCBI states that it places no restrictions
on use/distribution of molecular data, but explicitly notes possible third-party
rights and cannot grant blanket permission for them. Audit the exact source
inventory and retain applicable notices before hosting. This is not a statement
that an arbitrary nonviral bundle is cleared for redistribution.

Source: https://www.ncbi.nlm.nih.gov/home/about/policies/

Full-size builds and Nextflow execution on the target Linux host remain part of
the clean-install verification; local unit tests are not a substitute.
