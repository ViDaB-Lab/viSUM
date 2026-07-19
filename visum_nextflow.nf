#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NORMALIZE_FASTA } from './modules/local/normalize_fasta'
include { PREPARE_GENOMAD_DATABASE } from './modules/local/genomad_database'
include { RUN_GENOMAD } from './modules/local/run_genomad'
include { STANDARDIZE_GENOMAD } from './modules/local/standardize_genomad'
include { PREPARE_VIRSORTER2_DATABASE } from './modules/local/virsorter2_database'
include { RUN_VIRSORTER2 } from './modules/local/run_virsorter2'
include { STANDARDIZE_VIRSORTER2 } from './modules/local/standardize_virsorter2'


def validatePrefix(rawPrefix, source) {
    if( rawPrefix == null || rawPrefix.toString().isEmpty() ) {
        throw new IllegalArgumentException("Missing sample prefix (${source}).")
    }

    def prefix = rawPrefix.toString()
    def allowedPattern = '[A-Za-z0-9][A-Za-z0-9._-]*'

    if( !prefix.matches(allowedPattern) ) {
        def invalidChars = prefix
            .replaceAll('[A-Za-z0-9._-]', '')
            .toList()
            .unique()
            .join(' ')

        def suggestion = prefix.replaceAll('[^A-Za-z0-9._-]', '_')
        if( !suggestion.matches('[A-Za-z0-9].*') ) {
            suggestion = "sample${suggestion}"
        }

        def detail = invalidChars
            ? " Invalid character(s): '${invalidChars}'."
            : ' The first character must be a letter or number.'

        throw new IllegalArgumentException(
            "Invalid sample prefix '${prefix}' (${source}).${detail} " +
            'Use only letters, numbers, periods, underscores, and hyphens; ' +
            "the first character must be a letter or number. Suggested prefix: '${suggestion}'."
        )
    }

    return prefix
}


workflow {

    def singleMode = params.input != null
    def multiMode  = params.prefix_many != null

    if( singleMode && multiMode ) {
        error 'Choose one input mode: --input/--prefix/--type or --prefix_many/--indir.'
    }

    if( !singleMode && !multiMode ) {
        error 'Provide inputs with --input/--prefix/--type or --prefix_many/--indir.'
    }

    if( singleMode ) {
        if( !params.prefix ) {
            error 'Single mode requires --prefix.'
        }
        if( !params.type ) {
            error 'Single mode requires --type (dna or rna).'
        }

        def prefix = validatePrefix(params.prefix, '--prefix')
        def type = params.type.toString().trim().toLowerCase()
        if( !(type in ['dna', 'rna']) ) {
            error "Invalid --type '${params.type}'. Use dna or rna."
        }

        def fasta = file(params.input)
        if( !fasta.exists() ) {
            error "Input FASTA not found: ${params.input}"
        }

        ch_samples = Channel.of(tuple(prefix, type, fasta))
    }

    if( multiMode ) {
        if( !params.indir ) {
            error 'Multi mode requires --indir.'
        }

        def baseDir = file(params.indir)
        if( !baseDir.exists() || !baseDir.isDirectory() ) {
            error "--indir must be an existing directory: ${params.indir}"
        }

        def mapfile = file(params.prefix_many)
        if( !mapfile.exists() ) {
            error "prefix_many CSV not found: ${params.prefix_many}"
        }

        ch_samples = Channel
            .fromPath(mapfile)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                def prefix = row.prefix?.toString()
                def type = row.type?.toString()?.trim()?.toLowerCase()
                def relativeFasta = row.fasta?.toString()?.trim()

                if( !prefix || !type || !relativeFasta ) {
                    error "CSV rows require prefix,type,fasta values. Bad row: ${row}"
                }

                prefix = validatePrefix(prefix, 'prefix_many row')
                if( !(type in ['dna', 'rna']) ) {
                    error "Invalid type for prefix '${prefix}': '${type}'. Use dna or rna."
                }

                def fasta = file("${baseDir}/${relativeFasta}")
                if( !fasta.exists() ) {
                    error "FASTA not found for prefix '${prefix}': ${fasta}"
                }

                tuple(prefix, type, fasta)
            }
    }

    NORMALIZE_FASTA(ch_samples)

    NORMALIZE_FASTA.out.normalized_records.view { prefix, type, fasta, headerMap ->
        "NORMALIZED sample=${prefix} type=${type} fasta=${fasta.name} map=${headerMap.name}"
    }

    if( params.run_genomad ) {
        def userSuppliedDatabase = params.genomad_db != null
        def genomadDatabasePath = userSuppliedDatabase
            ? file(params.genomad_db).toString()
            : file("${params.dbdir}/genomad/genomad_db").toString()
        def genomadDatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_genomad_database_request = Channel.of(
            tuple(
                genomadDatabasePath,
                genomadDatabaseSource,
                params.genomad_auto_download
            )
        )

        PREPARE_GENOMAD_DATABASE(ch_genomad_database_request)

        PREPARE_GENOMAD_DATABASE.out.database.view { database, metadata ->
            "GENOMAD_DB database=${database} metadata=${metadata.name}"
        }

        ch_genomad_database = PREPARE_GENOMAD_DATABASE.out.database
            .map { database, metadata -> database }
            .first()

        RUN_GENOMAD(
            NORMALIZE_FASTA.out.normalized_records,
            ch_genomad_database
        )

        RUN_GENOMAD.out.results.view { prefix, type, review, virusSummary, virusFasta, virusGenes, virusProteins, plasmidSummary, metadata ->
            "GENOMAD sample=${prefix} type=${type} virus_summary=${virusSummary.name} output=${review.name}"
        }

        ch_genomad_for_standardizer = RUN_GENOMAD.out.results.map {
            prefix, type, review, virusSummary, virusFasta, virusGenes, virusProteins, plasmidSummary, metadata ->
                tuple(prefix, type, virusSummary, plasmidSummary, metadata)
        }

        ch_genomad_standardizer_input = ch_genomad_for_standardizer
            .join(
                NORMALIZE_FASTA.out.normalized_records.map { prefix, type, fasta, headerMap ->
                    tuple(prefix, headerMap)
                }
            )

        STANDARDIZE_GENOMAD(ch_genomad_standardizer_input)

        STANDARDIZE_GENOMAD.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
    }

    if( params.run_virsorter2 ) {
        def userSuppliedDatabase = params.virsorter2_db != null
        def virsorter2DatabasePath = userSuppliedDatabase
            ? file(params.virsorter2_db).toString()
            : file("${params.dbdir}/virsorter2/db").toString()
        def virsorter2DatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_virsorter2_database_request = Channel.of(
            tuple(
                virsorter2DatabasePath,
                virsorter2DatabaseSource,
                params.virsorter2_auto_download
            )
        )

        PREPARE_VIRSORTER2_DATABASE(ch_virsorter2_database_request)

        PREPARE_VIRSORTER2_DATABASE.out.database.view { database, metadata ->
            "VIRSORTER2_DB database=${database} metadata=${metadata.name}"
        }

        ch_virsorter2_database = PREPARE_VIRSORTER2_DATABASE.out.database
            .map { database, metadata -> database }
            .first()

        RUN_VIRSORTER2(
            NORMALIZE_FASTA.out.normalized_records,
            ch_virsorter2_database
        )

        RUN_VIRSORTER2.out.results.view {
            prefix, type, review, score, boundary, viralFasta, metadata ->
                "VIRSORTER2 sample=${prefix} type=${type} score=${score.name} output=${review.name}"
        }

        ch_virsorter2_for_standardizer = RUN_VIRSORTER2.out.results.map {
            prefix, type, review, score, boundary, viralFasta, metadata ->
                tuple(prefix, type, score, boundary, metadata)
        }

        ch_virsorter2_standardizer_input = ch_virsorter2_for_standardizer
            .join(
                NORMALIZE_FASTA.out.normalized_records.map {
                    prefix, type, fasta, headerMap -> tuple(prefix, headerMap)
                }
            )

        STANDARDIZE_VIRSORTER2(ch_virsorter2_standardizer_input)

        STANDARDIZE_VIRSORTER2.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
    }
}
