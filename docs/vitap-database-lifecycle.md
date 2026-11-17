# Select, reuse or update a VITAP database

[Database setup](database-setup.md)

VITAP is enabled by default. Prepare its assets with the ordinary `--setup`
command. On first installation, setup resolves the current ICTV VMR source and
uses VITAP's preparation path, which may download additional required references.

## Reuse

Normal setup reuses a valid managed `current` database without checking for a
new VMR. Analysis reads the selected completed database and never downloads or
updates it. To use an existing external installation, set `--vitap_db PATH` in
both setup and analysis; that database is not modified.

## Select a particular VMR

During setup, supply a local workbook/CSV with `--vitap_vmr PATH`. A release label
can be inferred from the filename or supplied with `--vitap_db_label LABEL`:

```bash
nextflow run visum_nextflow.nf --setup \
  --vitap_vmr /references/VMR-source.xlsx --vitap_db_label my-vmr-release
```

This explicitly selects that VMR even if managed `current` already exists.
Building may still need other network downloads. A label alone does not specify
the contents of an ICTV release; pair it with the intended VMR.

## Check for an update

```bash
nextflow run visum_nextflow.nf --setup --vitap_update_database true
```

This resolves the supplied VMR, or ICTV's current endpoint if none was supplied.
It is not a forced rebuild of a matching completed release. Other enabled tools
are also checked by setup. Update flags are not accepted in ordinary analysis.

## Preserve reference identity

A requested release is reused only when its recorded VMR checksum matches.
Conflicting/missing provenance fails rather than relabeling existing files.
Sources are stored under `sources/<SHA256>/<original-filename>`. Matching
interrupted work can be resumed; mismatching work requires inspection or a new
label. Existing completed legacy databases can be supplied explicitly.

Retain concrete release paths for analyses you need to reproduce. A matching
VMR checksum identifies that VMR, not every UniRef90/taxonomy source or software
environment used in construction. The VMR source, generated references and
runtime identities should all be recorded.

Updating VITAP does not automatically replace viSUM's separately pinned ICTV
rank table or reconciliation bundle. Do not change reference assets used by an
active analysis; prepare a new selection separately.
