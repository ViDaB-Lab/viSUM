# VITAP database selection

Normal reuse of a valid managed `current` database does not contact ICTV.
There is no automatic update check on each analysis.

## Existing interfaces

- `--vitap_db PATH`: validate and use this completed database without building
  or downloading. User-supplied databases are not modified.
- No explicit VMR or label, with `--vitap_update_database false`: validate and
  reuse local `current`. If none exists, the existing first-install path remains
  available when automatic downloading is enabled.
- `--vitap_vmr PATH`: explicitly select source contents, even when local `current`
  exists. No ICTV VMR download is needed, although building may download other
  required references.
- `--vitap_update_database true`: explicitly resolve a supplied VMR, or fetch
  ICTV's current VMR when no file was supplied. This is not a forced rebuild.
- `--vitap_db_label LABEL`: identify a separate release/build directory. A label
  does not independently identify an ICTV download; pair it with a specific VMR
  when selecting particular reference contents.

## Fingerprint and preservation guarantees

A requested release can be reused only when its recorded VMR SHA256 matches.
Missing or conflicting provenance causes an error before changing `current`.
Use a distinct label to build separately; never relabel an existing database's
provenance to claim it was built from a different source.

Source files are retained under `sources/<SHA256>/<original-filename>`, so files
with identical names but different contents do not overwrite one another.
Interrupted builds record their VMR fingerprint and reject mismatching or
unrecorded build work. Legacy interrupted builds require inspection or a new
label; they are not deleted automatically. Legacy completed databases remain
usable through ordinary local reuse or an explicit `--vitap_db` path.

A matching VMR fingerprint establishes VMR source identity, not full identity of
UniRef90, taxonomy downloads, software environments or every generated file.

## Separate setup/update workflow

Use `--setup` with `--run_vitap true` to run preparation without analysis inputs.
The same existing prepare process handles both modes. To explicitly check for an
update, also set `--vitap_update_database true`; setup alone does not enable updates.
Other enabled tools are also prepared; use the usual `--run_*` flags to select
them. See [setup](database-setup.md). Full clean-install verification remains pending.
