#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NORMALIZE_FASTA } from './modules/local/normalize_fasta'
include { PREPARE_GENOMAD_DATABASE } from './modules/local/genomad_database'
include { RUN_GENOMAD } from './modules/local/run_genomad'
include { STANDARDIZE_GENOMAD } from './modules/local/standardize_genomad'
include { PREPARE_VIRSORTER2_DATABASE } from './modules/local/virsorter2_database'
include { RUN_VIRSORTER2 } from './modules/local/run_virsorter2'
include { STANDARDIZE_VIRSORTER2 } from './modules/local/standardize_virsorter2'
include { PREPARE_CENOTETAKER3_DATABASE } from './modules/local/cenotetaker3_database'
include { RUN_CENOTETAKER3 } from './modules/local/run_cenotetaker3'
include { STANDARDIZE_CENOTETAKER3 } from './modules/local/standardize_cenotetaker3'
include { PREPARE_DEEP6_DATABASE } from './modules/local/deep6_database'
include { RUN_DEEP6 } from './modules/local/run_deep6'


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


def parseBooleanParameter(rawValue, parameterName) {
    if( rawValue instanceof Boolean ) {
        return rawValue
    }

    if( rawValue == null ) {
        throw new IllegalArgumentException(
            "Missing boolean value for ${parameterName}. Use true or false."
        )
    }

    def normalizedValue = rawValue.toString().trim().toLowerCase()
    if( normalizedValue == 'true' ) {
        return true
    }
    if( normalizedValue == 'false' ) {
        return false
    }

    throw new IllegalArgumentException(
        "Invalid value '${rawValue}' for ${parameterName}. Use true or false."
    )
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

    def runGenomad = parseBooleanParameter(params.run_genomad, '--run_genomad')
    def runVirsorter2 = parseBooleanParameter(params.run_virsorter2, '--run_virsorter2')
    def runCenotetaker3 = parseBooleanParameter(params.run_cenotetaker3, '--run_cenotetaker3')
    def runDeep6 = parseBooleanParameter(params.run_deep6, '--run_deep6')

    println(
        "viSUM program selection: " +
        "geNomad=${runGenomad}, " +
        "VirSorter2=${runVirsorter2}, " +
        "Cenote-Taker3=${runCenotetaker3}, " +
        "Deep6=${runDeep6}"
    )

    NORMALIZE_FASTA(ch_samples)

    NORMALIZE_FASTA.out.normalized_records.view { prefix, type, fasta, headerMap ->
        "NORMALIZED sample=${prefix} type=${type} fasta=${fasta.name} map=${headerMap.name}"
    }

    if( runGenomad ) {
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

    if( runVirsorter2 ) {
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

    if( runCenotetaker3 ) {
        def userSuppliedDatabase = params.ct3_db != null
        def cenotetaker3DatabasePath = userSuppliedDatabase
            ? file(params.ct3_db).toString()
            : file("${params.dbdir}/cenotetaker3/ct3_DBs").toString()
        def cenotetaker3DatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_cenotetaker3_database_request = Channel.of(
            tuple(
                cenotetaker3DatabasePath,
                cenotetaker3DatabaseSource,
                params.ct3_auto_download,
                params.ct3_hmm_db_version
            )
        )

        PREPARE_CENOTETAKER3_DATABASE(ch_cenotetaker3_database_request)

        PREPARE_CENOTETAKER3_DATABASE.out.database.view { database, metadata ->
            "CENOTETAKER3_DB database=${database} metadata=${metadata.name}"
        }

        ch_cenotetaker3_database = PREPARE_CENOTETAKER3_DATABASE.out.database
            .map { database, metadata -> database }
            .first()

        RUN_CENOTETAKER3(
            NORMALIZE_FASTA.out.normalized_records,
            ch_cenotetaker3_database
        )

        RUN_CENOTETAKER3.out.results.view {
            prefix, type, summary, virusFasta, virusProteins, pruneSummary,
            geneAnnotations, runArguments, log, metadata ->
                "CENOTETAKER3 sample=${prefix} type=${type} summary=${summary.name}"
        }

        ch_cenotetaker3_for_standardizer = RUN_CENOTETAKER3.out.results.map {
            prefix, type, summary, virusFasta, virusProteins, pruneSummary,
            geneAnnotations, runArguments, log, metadata ->
                tuple(
                    prefix,
                    type,
                    summary,
                    virusFasta,
                    pruneSummary,
                    geneAnnotations,
                    metadata
                )
        }

        ch_cenotetaker3_standardizer_input = ch_cenotetaker3_for_standardizer
            .join(
                NORMALIZE_FASTA.out.normalized_records.map {
                    prefix, type, fasta, headerMap -> tuple(prefix, headerMap)
                }
            )

        STANDARDIZE_CENOTETAKER3(ch_cenotetaker3_standardizer_input)

        STANDARDIZE_CENOTETAKER3.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
    }

    if( runDeep6 ) {
        def deep6MinimumLength
        try {
            deep6MinimumLength = params.deep6_minlen as Integer
        }
        catch( Exception ignored ) {
            error "Invalid --deep6_minlen '${params.deep6_minlen}'. Use an integer of at least 250."
        }
        if( deep6MinimumLength < 250 ) {
            error "Invalid --deep6_minlen '${params.deep6_minlen}'. Deep6 requires at least 250 nt."
        }

        def userSuppliedInstallation = params.deep6_dir != null
        def deep6InstallationPath = userSuppliedInstallation
            ? file(params.deep6_dir).toString()
            : file("${params.dbdir}/deep6/Deep6").toString()
        def deep6InstallationSource = userSuppliedInstallation
            ? 'user-supplied'
            : 'viSUM-managed'

        def userSuppliedModels = params.deep6_model != null
        def deep6ModelPath = userSuppliedModels
            ? file(params.deep6_model).toString()
            : file("${deep6InstallationPath}/Models").toString()
        def deep6ModelSource = userSuppliedModels
            ? 'user-supplied'
            : 'bundled-with-installation'

        ch_deep6_samples = NORMALIZE_FASTA.out.normalized_records.filter {
            prefix, type, fasta, headerMap -> type == 'rna'
        }

        // Trigger one shared setup only when at least one RNA sample is present.
        ch_deep6_database_request = ch_deep6_samples
            .take(1)
            .map { prefix, type, fasta, headerMap ->
                tuple(
                    deep6InstallationPath,
                    deep6InstallationSource,
                    deep6ModelPath,
                    deep6ModelSource,
                    params.deep6_auto_download,
                    params.deep6_repository,
                    params.deep6_revision
                )
            }

        PREPARE_DEEP6_DATABASE(ch_deep6_database_request)

        PREPARE_DEEP6_DATABASE.out.database.view {
            installation, models, metadata ->
                "DEEP6_DB installation=${installation} models=${models} metadata=${metadata.name}"
        }

        ch_deep6_bundle = PREPARE_DEEP6_DATABASE.out.database
            .map { installation, models, metadata ->
                tuple(installation, models, metadata)
            }
            .first()

        RUN_DEEP6(ch_deep6_samples, ch_deep6_bundle)

        RUN_DEEP6.out.results.view { prefix, type, scores, logFile, metadata ->
            "DEEP6 sample=${prefix} type=${type} scores=${scores.name}"
        }
    }
}
