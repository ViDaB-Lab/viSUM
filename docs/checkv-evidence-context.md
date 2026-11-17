# CheckV quality and supporting evidence

[Interpreting results](interpretation.md) · [Output columns](reporting-bundle.md)

Completeness estimates and gene-content evidence answer different questions.
`checkv_quality=Not-determined` means completeness is unresolved; it does not by
itself mean nonviral. Viral/host gene counts and warnings remain available even
when no completeness estimate is possible.

- No informative viral or host genes: weak, unclassified evidence.
- Host genes without viral genes: qualified cellular evidence.
- Viral genes with undetermined quality: weak viral evidence, retaining host context.

Check the quality table's assessment scope before applying a value to an
extracted candidate. Parent-only and proviral aggregate context cannot establish
a child's completeness. Direct `Complete` assessments give a true complete flag;
determined high/medium/low categories give false, and undetermined/missing or
non-direct context gives `NA`. See the exact field definitions in the reporting guide.

## RNA candidates with CheckV-only cellular conflict

An intact RNA candidate can remain retained despite CheckV-only cellular
conflict when Deep6, VirBot and viCAT each provide qualified/strong support,
there is no plasmid/mobile conflict, and no viCAT locus is cellular or ambiguous.
At least one viral-supported protein must have length >=100 aa, bit score >=100,
identity >=50%, query/subject coverage >=80%, and a coding interval spanning
>=70% of the input. Other safeguards still apply. DNA and extracted regions are
not eligible, and missing supporting loci cannot satisfy this exception.

The CheckV conflict is not erased. Detailed metadata records
`rna_checkv_conflict_exception` and `rna_checkv_exception_supporting_loci`; the
harmonizer manifest records thresholds and affected IDs. These rules do not
establish that the sequence is a complete or infectious virus.
