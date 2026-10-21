# viCAT nonviral assembly inventory

`assembly_manifest.tsv` is the default input list for automatic viCAT nonviral
source acquisition. It contains 480 versioned RefSeq assembly accessions, their
cellular groups, organism names and source taxids. Taxids are source annotations,
not a claim that the list represents 480 distinct species.

Inventory SHA-256:
`7039630fd08289d871c5c714ff3eaf4fe96a4f164063b7fa73dc6e40f6cf2a29`.

The list was recovered from the complete source catalog for the historical
`nonviral_classaware_v1` database. No assembly selections or classification rules
were changed. The source catalog SHA-256 is
`dbc7cf8fffecad85e65d4670300691f2cc30f3a4e9ac9848c73e9640429dc9c8`.

Setup downloads proteins and feature tables for these assemblies, records source
checksums, classifies the references and builds the nonviral search database.
Analysis uses that prepared database; it does not download sources.

This directory is an accession inventory, not a prebuilt database or permission
to redistribute source sequences. Third-party data retain their applicable terms.
Versioned assembly accessions do not freeze later annotation changes, so a new
download is not guaranteed to reproduce the bytes of the historical database.
Retain setup's downloaded-source receipts and database build records with a run.

See [viCAT database setup](../../../docs/vicat-database-preparation.md).
