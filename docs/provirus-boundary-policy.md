# Provirus boundary authority

By default, only geNomad or corroborated CheckV boundaries can change sequence
coordinates. CT3 discovery calls, gene annotations and proposed regions remain
available as evidence; a CT3-only region does not authorize trimming.

Boundary selection uses geNomad first. A CheckV boundary requires overlapping
CT3 support or qualifying viCAT regional support (the configured overlap
threshold remains 0.5). viCAT is advisory, not a boundary authority. Unsupported
regions remain recorded in the boundary audit, and their parent sequence remains
unchanged at refinement. Final retention is still decided by viHARMONY: preserving
a parent at refinement does not guarantee that it enters the primary database.

`--allow_ct3_only_refinement true` remains an explicit experimental opt-in.
The pipeline default and benchmark launcher set it to `false`.

When comparing policies, report both parent-level detection and intact sequence
recovery. A correctly detected viral parent can still be incorrectly shortened.
Replay refinement from saved discovery evidence first. If sequence IDs or contents
change, regenerate affected downstream evidence before claiming final sensitivity,
FPR or taxonomy accuracy. Do not attach the former subregion's taxonomy evidence
to the newly restored full parent as though it had been rerun.
