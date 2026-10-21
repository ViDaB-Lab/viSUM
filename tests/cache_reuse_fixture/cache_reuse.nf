nextflow.enable.dsl=2

// Focused cache contract, NOT a second biological workflow. Analysis modules
// are copied verbatim from viSUM. Only their external executables are fakes.
include { NORMALIZE_FASTA } from './modules/local/normalize_fasta'
include { RUN_GENOMAD } from './modules/local/run_genomad'
include { RUN_VIRSORTER2 } from './modules/local/run_virsorter2'
include { RUN_VICAT_DIAMOND } from './modules/local/run_vicat_diamond'
include { RUN_VICAT_NONVIRAL_DIAMOND } from './modules/local/run_vicat_nonviral_diamond'
include { VIHARMONY } from './modules/local/viharmony'
include { PREPARE_GENOMAD_DATABASE } from './modules/local/genomad_database'
include { PREPARE_VICAT_DATABASE } from './modules/local/vicat_database'
include { PREPARE_VICAT_NONVIRAL_DATABASE } from './modules/local/vicat_nonviral_database'

workflow {
    // Nextflow 26.04 CLI values are strings: never use their Groovy truthiness.
    // The harness supplies a genuine JSON Boolean through -params-file.
    if (!(params.fixture_inline_prepare instanceof Boolean)) {
        error 'fixture_inline_prepare must be a Boolean supplied through -params-file'
    }
    log.info("CACHE_FIXTURE inline_prepare=${params.fixture_inline_prepare}")
    def data = file("${projectDir}/data")
    def assets = file("${projectDir}/databases")
    NORMALIZE_FASTA(Channel.of(
        tuple('sample_dna', 'dna', file("${data}/dna.fasta")),
        tuple('sample_rna', 'rna', file("${data}/rna.fasta"))
    ))

    if (params.fixture_inline_prepare) {
        // Preserve prepare inputs, outputs and cache directives. Test-only
        // bodies create symlinks/metadata without installing anything.
        PREPARE_GENOMAD_DATABASE(Channel.value(tuple("${assets}/genomad", 'user-supplied', false)))
        PREPARE_VICAT_DATABASE(Channel.value(tuple(
            "${assets}/vicat", 'user-supplied', '', '', '', '', '', '', 0, 0, 0, '', '', 0.9, 0.8, 50
        )))
        PREPARE_VICAT_NONVIRAL_DATABASE(Channel.value(tuple(
            "${assets}/nonviral", 'user-supplied', '', '', '', '', file("${data}/helper.txt")
        )))
        genomad = PREPARE_GENOMAD_DATABASE.out.database.map { db, metadata -> db }.first()
        vicat = PREPARE_VICAT_DATABASE.out.database.map { db, metadata -> db }.first()
        nonviral = PREPARE_VICAT_NONVIRAL_DATABASE.out.database.map { db, metadata -> db }.first()
    } else {
        genomad = Channel.value(file("${assets}/genomad"))
        vicat = Channel.value(file("${assets}/vicat"))
        nonviral = Channel.value(file("${assets}/nonviral"))
    }

    RUN_GENOMAD(NORMALIZE_FASTA.out.normalized_records, genomad)
    RUN_VIRSORTER2(NORMALIZE_FASTA.out.normalized_records, Channel.value(file("${assets}/virsorter2")))
    orfs = NORMALIZE_FASTA.out.normalized_records.map { prefix, type, fasta, header ->
        tuple(prefix, type, file("${data}/${prefix}.orfs.faa"), file("${data}/orf_map.tsv"), header)
    }
    RUN_VICAT_DIAMOND(orfs, vicat)
    RUN_VICAT_NONVIRAL_DIAMOND(orfs, nonviral)

    // Actual output paths are inputs downstream: invalidation must propagate.
    completed = RUN_GENOMAD.out.results.map { row -> tuple(row[0], row[3]) }
        .join(RUN_VIRSORTER2.out.results.map { row -> tuple(row[0], row[3]) })
        .join(RUN_VICAT_DIAMOND.out.results.map { row -> tuple(row[0], row[4]) })
        .join(RUN_VICAT_NONVIRAL_DIAMOND.out.results.map { row -> tuple(row[0], row[4]) })
    harmony = NORMALIZE_FASTA.out.normalized_records.join(completed)
        .map { prefix, type, fasta, header, geno, vs2, viral, nonviralHits ->
            tuple(prefix, type, fasta, header, file("${data}/gate.tsv"),
                file("${data}/refined.fasta"), file("${data}/regions.tsv"),
                4, [geno, vs2, viral, nonviralHits], 0, file("${data}/groups.tsv"),
                file("${data}/msl.csv"), 'compact', 1000)
        }
    VIHARMONY(harmony)
}
