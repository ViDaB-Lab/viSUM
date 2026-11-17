# Candidate provirus boundaries

[Interpreting results](interpretation.md)

By default, geNomad boundaries or corroborated CheckV boundaries can change
sequence coordinates. geNomad is considered first. CheckV boundaries need
overlapping Cenote-Taker3 support or qualifying viCAT regional support; the
configured viCAT overlap threshold is 0.5.

Cenote-Taker3-only regions are advisory by default. viCAT also provides regional
support but does not supply an independent extraction boundary. Unsupported
proposed regions stay in the boundary audit while the parent remains unchanged
at refinement. A preserved parent still has to pass viHARMONY's origin rules to
be exported.

`--allow_ct3_only_refinement true` is an explicit experimental opt-in. The
default is false. Use the default unless you deliberately want different boundary
authority and understand the effect on downstream sequences and evidence.

Inspect record type, parent ID and coordinates before using a derived sequence.
A computational provirus call does not prove integration experimentally. If
boundaries or sequence contents change, old subregion evidence must not simply
be treated as evidence for a newly restored parent; affected downstream tasks
need to run on the appropriate sequences. Keep the task cache and let native
Nextflow resolve dependencies rather than editing cached files.
