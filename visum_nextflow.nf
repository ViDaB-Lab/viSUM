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
include { STANDARDIZE_DEEP6 } from './modules/local/standardize_deep6'
include { PREPARE_DEEPMICROCLASS2 } from './modules/local/deepmicroclass2_installation'
include { RUN_DEEPMICROCLASS2 } from './modules/local/run_deepmicroclass2'
include { STANDARDIZE_DEEPMICROCLASS2 } from './modules/local/standardize_deepmicroclass2'
include { PREPARE_VIRBOT_DATABASE } from './modules/local/virbot_installation'
include { RUN_VIRBOT } from './modules/local/run_virbot'
include { STANDARDIZE_VIRBOT } from './modules/local/standardize_virbot'
include { PREPARE_GIANTHUNTER_DATABASE } from './modules/local/gianthunter_database'
include { RUN_GIANTHUNTER } from './modules/local/run_gianthunter'
include { STANDARDIZE_GIANTHUNTER } from './modules/local/standardize_gianthunter'
include { PREPARE_VICAT_DATABASE } from './modules/local/vicat_database'
include { PREDICT_VICAT_ORFS } from './modules/local/predict_vicat_orfs'
include { RUN_VICAT_DIAMOND } from './modules/local/run_vicat_diamond'
include { STANDARDIZE_VICAT } from './modules/local/standardize_vicat'
include { DISCOVERY_GATE } from './modules/local/discovery_gate'
include { PREPARE_CHECKV_DATABASE } from './modules/local/checkv_database'
include { RUN_CHECKV } from './modules/local/run_checkv'
include { STANDARDIZE_CHECKV } from './modules/local/standardize_checkv'
include { REFINE_PROVIRAL_REGIONS } from './modules/local/refine_proviral_regions'
include { PREPARE_VITAP_DATABASE } from './modules/local/vitap_database'
include { RUN_VITAP } from './modules/local/run_vitap'
include { STANDARDIZE_VITAP } from './modules/local/standardize_vitap'
include { PREPARE_VCONTACT3_DATABASE } from './modules/local/vcontact3_database'


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


def parsePositiveIntegerParameter(rawValue, parameterName) {
    Integer parsedValue
    try {
        parsedValue = rawValue as Integer
    }
    catch( Exception ignored ) {
        throw new IllegalArgumentException(
            "Invalid value '${rawValue}' for ${parameterName}. Use a positive integer."
        )
    }

    if( parsedValue < 1 ) {
        throw new IllegalArgumentException(
            "Invalid value '${rawValue}' for ${parameterName}. Use a positive integer."
        )
    }
    return parsedValue
}


workflow {

    if( params.threads != null ) {
        error '--threads is no longer used because it could not distinguish a total workflow budget from per-tool threads. Use --max_cpus and, when needed, a tool-specific --*_cpus option.'
    }
    if( params.vicat_build_threads != null ) {
        error '--vicat_build_threads was renamed to --vicat_build_cpus.'
    }
    if( params.vicat_orf_threads != null ) {
        error '--vicat_orf_threads was renamed to --vicat_orf_cpus.'
    }
    if( params.vicat_threads != null ) {
        error '--vicat_threads was renamed to --vicat_cpus.'
    }

    def maxCpus = parsePositiveIntegerParameter(params.max_cpus, '--max_cpus')
    def analysisCpuParameters = [
        '--genomad_cpus': params.genomad_cpus,
        '--virsorter2_cpus': params.virsorter2_cpus,
        '--ct3_cpus': params.ct3_cpus,
        '--deep6_cpus': params.deep6_cpus,
        '--deepmicroclass2_cpus': params.deepmicroclass2_cpus,
        '--virbot_cpus': params.virbot_cpus,
        '--gianthunter_cpus': params.gianthunter_cpus,
        '--vicat_orf_cpus': params.vicat_orf_cpus,
        '--vicat_cpus': params.vicat_cpus,
        '--checkv_cpus': params.checkv_cpus,
        '--vitap_cpus': params.vitap_cpus,
        '--vcontact3_cpus': params.vcontact3_cpus,
    ]
    analysisCpuParameters.each { parameterName, rawValue ->
        parsePositiveIntegerParameter(rawValue, parameterName)
    }

    println(
        "viSUM resource budget: max_cpus=${maxCpus}, max_memory=${params.max_memory}; " +
        "preferred tool CPUs (each capped at max_cpus)=" + analysisCpuParameters.collect { name, value ->
            "${name.substring(2)}=${value}"
        }.join(', ')
    )

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
    def runVirbot = parseBooleanParameter(params.run_virbot, '--run_virbot')
    def runGianthunter = parseBooleanParameter(
        params.run_gianthunter,
        '--run_gianthunter'
    )
    def runVicat = parseBooleanParameter(params.run_vicat, '--run_vicat')
    def runDeepmicroclass2 = parseBooleanParameter(
        params.run_deepmicroclass2,
        '--run_deepmicroclass2'
    )
    def runCheckv = parseBooleanParameter(params.run_checkv, '--run_checkv')
    def runVitap = parseBooleanParameter(params.run_vitap, '--run_vitap')
    def runVcontact3 = parseBooleanParameter(
        params.run_vcontact3,
        '--run_vcontact3'
    )
    def vitapIncludeLowConfidence = parseBooleanParameter(
        params.vitap_include_low_confidence,
        '--vitap_include_low_confidence'
    )
    def allowCt3OnlyRefinement = parseBooleanParameter(
        params.allow_ct3_only_refinement,
        '--allow_ct3_only_refinement'
    )

    // Every standardizer emits sparse, threshold-qualified evidence. These
    // channels are merged and grouped by sample for the discovery gate.
    ch_discovery_evidence = Channel.empty()
    // Only boundary-capable tools contribute to post-discovery provirus
    // refinement. Their normal evidence tables are reused directly.
    ch_provirus_evidence = Channel.empty()

    println(
        "viSUM program selection: " +
        "geNomad=${runGenomad}, " +
        "VirSorter2=${runVirsorter2}, " +
        "Cenote-Taker3=${runCenotetaker3}, " +
        "Deep6=${runDeep6}, " +
        "VirBot=${runVirbot}, " +
        "DeepMicroClass2=${runDeepmicroclass2}, " +
        "GiantHunter=${runGianthunter}, " +
        "viCAT=${runVicat}, " +
        "CheckV=${runCheckv}, " +
        "VITAP=${runVitap}, " +
        "vConTACT3=${runVcontact3}"
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
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_GENOMAD.out.evidence
        )
        ch_provirus_evidence = ch_provirus_evidence.mix(
            STANDARDIZE_GENOMAD.out.evidence
        )
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
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_VIRSORTER2.out.evidence
        )
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
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_CENOTETAKER3.out.evidence
        )
        ch_provirus_evidence = ch_provirus_evidence.mix(
            STANDARDIZE_CENOTETAKER3.out.evidence
        )
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

        def deep6MinimumScore
        try {
            deep6MinimumScore = params.deep6_min_score as Double
        }
        catch( Exception ignored ) {
            error "Invalid --deep6_min_score '${params.deep6_min_score}'. Use a number from 0 to 1."
        }
        if( !Double.isFinite(deep6MinimumScore) || deep6MinimumScore < 0.0 || deep6MinimumScore > 1.0 ) {
            error "Invalid --deep6_min_score '${params.deep6_min_score}'. Use a number from 0 to 1."
        }

        def deep6MedianMultiplier
        try {
            deep6MedianMultiplier = params.deep6_median_multiplier as Double
        }
        catch( Exception ignored ) {
            error "Invalid --deep6_median_multiplier '${params.deep6_median_multiplier}'. Use a number of at least 1."
        }
        if( !Double.isFinite(deep6MedianMultiplier) || deep6MedianMultiplier < 1.0 ) {
            error "Invalid --deep6_median_multiplier '${params.deep6_median_multiplier}'. Use a number of at least 1."
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

        ch_deep6_for_standardizer = RUN_DEEP6.out.results.map {
            prefix, type, scores, logFile, metadata ->
                tuple(prefix, type, scores, metadata)
        }

        ch_deep6_standardizer_input = ch_deep6_for_standardizer
            .join(
                NORMALIZE_FASTA.out.normalized_records.map {
                    prefix, type, fasta, headerMap -> tuple(prefix, headerMap)
                }
            )

        STANDARDIZE_DEEP6(ch_deep6_standardizer_input)

        STANDARDIZE_DEEP6.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_DEEP6.out.evidence
        )
    }

    if( runDeepmicroclass2 ) {
        def deepmicroclass2Model = params.deepmicroclass2_model?.toString()?.trim()
        if( !(deepmicroclass2Model in ['8class', 'high_precision', '300bp']) ) {
            error "Invalid --deepmicroclass2_model '${params.deepmicroclass2_model}'. Use 8class, high_precision, or 300bp."
        }
        def deepmicroclass2ModelDescription = [
            '8class': 'fastest; length-matched non-overlapping windows; minimum 500 nt',
            'high_precision': 'overlapping length-matched windows; includes 300-499 nt contigs',
            '300bp': 'advanced; applies the 300-nt model to every contig of at least 300 nt'
        ][deepmicroclass2Model]
        println(
            "DeepMicroClass2 model: ${deepmicroclass2Model} " +
            "(${deepmicroclass2ModelDescription})"
        )

        def userSuppliedInstallation = params.deepmicroclass2_dir != null
        def deepmicroclass2InstallationPath = userSuppliedInstallation
            ? file(params.deepmicroclass2_dir).toString()
            : file("${params.tooldir}/deepmicroclass2/DeepMicroClass2").toString()
        def deepmicroclass2InstallationSource = userSuppliedInstallation
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_deepmicroclass2_samples = NORMALIZE_FASTA.out.normalized_records.filter {
            prefix, type, fasta, headerMap -> type == 'dna'
        }

        // Trigger one shared installation check only when at least one DNA
        // sample is available for DeepMicroClass2.
        ch_deepmicroclass2_installation_request = ch_deepmicroclass2_samples
            .take(1)
            .map { prefix, type, fasta, headerMap ->
                tuple(
                    deepmicroclass2InstallationPath,
                    deepmicroclass2InstallationSource,
                    params.deepmicroclass2_auto_download,
                    params.deepmicroclass2_repository,
                    params.deepmicroclass2_revision
                )
            }

        PREPARE_DEEPMICROCLASS2(ch_deepmicroclass2_installation_request)

        PREPARE_DEEPMICROCLASS2.out.installation.view {
            installation, metadata ->
                "DEEPMICROCLASS2_INSTALL installation=${installation} metadata=${metadata.name}"
        }

        ch_deepmicroclass2_bundle = PREPARE_DEEPMICROCLASS2.out.installation
            .map { installation, metadata -> tuple(installation, metadata) }
            .first()

        RUN_DEEPMICROCLASS2(
            ch_deepmicroclass2_samples,
            ch_deepmicroclass2_bundle
        )

        RUN_DEEPMICROCLASS2.out.results.view {
            prefix, type, scores, logFile, metadata ->
                "DEEPMICROCLASS2 sample=${prefix} type=${type} scores=${scores.name}"
        }

        ch_deepmicroclass2_for_standardizer = RUN_DEEPMICROCLASS2.out.results.map {
            prefix, type, scores, logFile, metadata ->
                tuple(prefix, type, scores, metadata)
        }

        ch_deepmicroclass2_standardizer_input = ch_deepmicroclass2_for_standardizer
            .join(
                NORMALIZE_FASTA.out.normalized_records.map {
                    prefix, type, fasta, headerMap -> tuple(prefix, headerMap)
                }
            )

        STANDARDIZE_DEEPMICROCLASS2(ch_deepmicroclass2_standardizer_input)

        STANDARDIZE_DEEPMICROCLASS2.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_DEEPMICROCLASS2.out.evidence
        )
    }

    if( runVirbot ) {
        def virbotSensitive = parseBooleanParameter(
            params.virbot_sensitive,
            '--virbot_sensitive'
        )
        def virbotTaxaMode = params.virbot_taxa?.toString()?.trim()?.toUpperCase()
        if( !(virbotTaxaMode in ['TOP', 'LCA']) ) {
            error "Invalid --virbot_taxa '${params.virbot_taxa}'. Use TOP or LCA."
        }
        def virbotIctvCsv = file(params.ictv_csv)
        if( !virbotIctvCsv.exists() ) {
            error "ICTV taxonomy file required by VirBot was not found: ${params.ictv_csv}"
        }

        def userSuppliedInstallation = params.virbot_dir != null
        def virbotInstallationPath = userSuppliedInstallation
            ? file(params.virbot_dir).toString()
            : file("${params.dbdir}/virbot/VirBot").toString()
        def virbotInstallationSource = userSuppliedInstallation
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_virbot_samples = NORMALIZE_FASTA.out.normalized_records.filter {
            prefix, type, fasta, headerMap -> type == 'rna'
        }

        // No RNA sample means no VirBot installation check or analysis task.
        ch_virbot_installation_request = ch_virbot_samples
            .take(1)
            .map { prefix, type, fasta, headerMap ->
                tuple(
                    virbotInstallationPath,
                    virbotInstallationSource,
                    params.virbot_auto_download,
                    params.virbot_repository,
                    params.virbot_revision
                )
            }

        PREPARE_VIRBOT_DATABASE(ch_virbot_installation_request)

        PREPARE_VIRBOT_DATABASE.out.installation.view {
            installation, database, metadata ->
                "VIRBOT_INSTALL installation=${installation} database=${database} metadata=${metadata.name}"
        }

        ch_virbot_bundle = PREPARE_VIRBOT_DATABASE.out.installation
            .map { installation, database, metadata ->
                tuple(
                    installation,
                    database,
                    metadata,
                    virbotSensitive,
                    virbotTaxaMode
                )
            }
            .first()

        RUN_VIRBOT(ch_virbot_samples, ch_virbot_bundle)

        RUN_VIRBOT.out.results.view {
            prefix, type, scores, virusFasta, logFile, metadata ->
                "VIRBOT sample=${prefix} type=${type} scores=${scores.name} virus_fasta=${virusFasta.name}"
        }

        ch_virbot_for_standardizer = RUN_VIRBOT.out.results.map {
            prefix, type, scores, virusFasta, logFile, metadata ->
                tuple(prefix, type, scores, virusFasta, metadata)
        }

        ch_virbot_standardizer_input = ch_virbot_for_standardizer
            .join(
                NORMALIZE_FASTA.out.normalized_records.map {
                    prefix, type, fasta, headerMap -> tuple(prefix, headerMap)
                }
            )

        STANDARDIZE_VIRBOT(
            ch_virbot_standardizer_input,
            virbotIctvCsv
        )

        STANDARDIZE_VIRBOT.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_VIRBOT.out.evidence
        )
    }

    if( runGianthunter ) {
        def gianthunterIctvCsv = file(params.ictv_csv)
        if( !gianthunterIctvCsv.exists() ) {
            error "GiantHunter ICTV taxonomy file not found: ${gianthunterIctvCsv}"
        }

        def gianthunterMinimumLength
        try {
            gianthunterMinimumLength = params.gianthunter_min_length as Integer
        }
        catch( Exception ignored ) {
            error "Invalid --gianthunter_min_length '${params.gianthunter_min_length}'. Use a positive integer."
        }
        if( gianthunterMinimumLength < 1 ) {
            error "Invalid --gianthunter_min_length '${params.gianthunter_min_length}'. Use a positive integer."
        }

        def gianthunterReject
        try {
            gianthunterReject = params.gianthunter_reject as Double
        }
        catch( Exception ignored ) {
            error "Invalid --gianthunter_reject '${params.gianthunter_reject}'. Use a number from 0 to 1."
        }
        if( !Double.isFinite(gianthunterReject) || gianthunterReject < 0.0 || gianthunterReject > 1.0 ) {
            error "Invalid --gianthunter_reject '${params.gianthunter_reject}'. Use a number from 0 to 1."
        }

        def gianthunterQueryCover
        try {
            gianthunterQueryCover = params.gianthunter_query_cover as Integer
        }
        catch( Exception ignored ) {
            error "Invalid --gianthunter_query_cover '${params.gianthunter_query_cover}'. Use an integer from 0 to 100."
        }
        if( gianthunterQueryCover < 0 || gianthunterQueryCover > 100 ) {
            error "Invalid --gianthunter_query_cover '${params.gianthunter_query_cover}'. Use an integer from 0 to 100."
        }

        def userSuppliedDatabase = params.gianthunter_db != null
        def gianthunterDatabasePath = userSuppliedDatabase
            ? file(params.gianthunter_db).toString()
            : file("${params.dbdir}/gianthunter/gianthunter_db_v1").toString()
        def gianthunterDatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_gianthunter_samples = NORMALIZE_FASTA.out.normalized_records.filter {
            prefix, type, fasta, headerMap -> type == 'dna'
        }

        // No DNA sample means no GiantHunter database check or analysis task.
        ch_gianthunter_database_request = ch_gianthunter_samples
            .take(1)
            .map { prefix, type, fasta, headerMap ->
                tuple(
                    gianthunterDatabasePath,
                    gianthunterDatabaseSource,
                    params.gianthunter_auto_download,
                    params.gianthunter_database_release,
                    params.gianthunter_database_url,
                    params.gianthunter_database_sha256
                )
            }

        PREPARE_GIANTHUNTER_DATABASE(ch_gianthunter_database_request)

        PREPARE_GIANTHUNTER_DATABASE.out.database.view { database, metadata ->
            "GIANTHUNTER_DB database=${database} metadata=${metadata.name}"
        }

        ch_gianthunter_bundle = PREPARE_GIANTHUNTER_DATABASE.out.database
            .map { database, metadata -> tuple(database, metadata) }
            .first()

        RUN_GIANTHUNTER(
            ch_gianthunter_samples,
            ch_gianthunter_bundle
        )

        RUN_GIANTHUNTER.out.results.view {
            prefix, type, prediction, annotations, logFile, metadata ->
                "GIANTHUNTER sample=${prefix} type=${type} prediction=${prediction.name}"
        }

        ch_gianthunter_for_standardizer = RUN_GIANTHUNTER.out.results.map {
            prefix, type, prediction, annotations, logFile, metadata ->
                tuple(prefix, type, prediction, annotations, metadata)
        }

        ch_gianthunter_standardizer_input = ch_gianthunter_for_standardizer
            .join(
                NORMALIZE_FASTA.out.normalized_records.map {
                    prefix, type, fasta, headerMap -> tuple(prefix, headerMap)
                }
            )

        STANDARDIZE_GIANTHUNTER(
            ch_gianthunter_standardizer_input,
            gianthunterIctvCsv
        )

        STANDARDIZE_GIANTHUNTER.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_GIANTHUNTER.out.evidence
        )
    }

    if( runVicat ) {
        def allowedVicatSensitivity = [
            'faster', 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive',
            'very-sensitive', 'ultra-sensitive'
        ]
        if( !(params.vicat_diamond_sensitivity in allowedVicatSensitivity) ) {
            error "Invalid --vicat_diamond_sensitivity '${params.vicat_diamond_sensitivity}'. " +
                "Use one of: ${allowedVicatSensitivity.join(', ')}."
        }
        def vicatOrfTaxonomySupport = params.vicat_orf_taxonomy_support as Double
        if( !Double.isFinite(vicatOrfTaxonomySupport) ||
            vicatOrfTaxonomySupport < 0.5 || vicatOrfTaxonomySupport > 1.0 ) {
            error '--vicat_orf_taxonomy_support must be between 0.5 and 1.0.'
        }
        def vicatContigTaxonomySupport = params.vicat_contig_taxonomy_support as Double
        if( !Double.isFinite(vicatContigTaxonomySupport) ||
            vicatContigTaxonomySupport < 0.5 || vicatContigTaxonomySupport > 1.0 ) {
            error '--vicat_contig_taxonomy_support must be between 0.5 and 1.0.'
        }
        def vicatLocusOverlap = params.vicat_locus_overlap as Double
        if( !Double.isFinite(vicatLocusOverlap) ||
            vicatLocusOverlap <= 0.0 || vicatLocusOverlap > 1.0 ) {
            error '--vicat_locus_overlap must be greater than 0 and at most 1.0.'
        }
        def vicatMinimumQueryCover = params.vicat_min_query_cover as Integer
        if( vicatMinimumQueryCover < 0 || vicatMinimumQueryCover > 100 ) {
            error '--vicat_min_query_cover must be between 0 and 100.'
        }
        def vicatMinimumBitscore = params.vicat_min_bitscore as Double
        if( !Double.isFinite(vicatMinimumBitscore) || vicatMinimumBitscore < 0.0 ) {
            error '--vicat_min_bitscore must be zero or greater.'
        }
        def vicatTopPercent = params.vicat_top_percent as Double
        if( vicatTopPercent < 0.0 || vicatTopPercent > 100.0 ) {
            error '--vicat_top_percent must be between 0 and 100.'
        }
        def vicatBlockSize = params.vicat_block_size as Double
        if( !Double.isFinite(vicatBlockSize) || vicatBlockSize <= 0.0 ) {
            error '--vicat_block_size must be greater than zero.'
        }
        def vicatIndexChunks = params.vicat_index_chunks as Integer
        if( vicatIndexChunks < 1 ) {
            error '--vicat_index_chunks must be at least 1.'
        }

        def userSuppliedDatabase = params.vicat_db != null
        def vicatDatabasePath = userSuppliedDatabase
            ? file(params.vicat_db).toString()
            : file(params.vicat_managed_db).toString()
        def vicatDatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        def vicatSourceProteins = params.vicat_metavr_proteins == null
            ? ''
            : file(params.vicat_metavr_proteins).toString()
        def vicatSourceMetadata = params.vicat_metavr_metadata == null
            ? ''
            : file(params.vicat_metavr_metadata).toString()

        if( (vicatSourceProteins && !vicatSourceMetadata) ||
            (!vicatSourceProteins && vicatSourceMetadata) ) {
            error 'A local viCAT build requires both --vicat_metavr_proteins and --vicat_metavr_metadata.'
        }

        ch_vicat_database_request = NORMALIZE_FASTA.out.normalized_records
            .take(1)
            .map { prefix, type, fasta, headerMap ->
                tuple(
                    vicatDatabasePath,
                    vicatDatabaseSource,
                    vicatSourceProteins,
                    vicatSourceMetadata,
                    params.vicat_metavr_proteins_sha256,
                    params.vicat_metavr_metadata_sha256,
                    params.vicat_expected_uvigs,
                    params.vicat_expected_proteins,
                    params.vicat_expected_representatives
                )
            }

        PREPARE_VICAT_DATABASE(ch_vicat_database_request)

        PREPARE_VICAT_DATABASE.out.database.view { database, metadata ->
            "VICAT_DB database=${database} metadata=${metadata.name}"
        }

        ch_vicat_database = PREPARE_VICAT_DATABASE.out.database
            .map { database, metadata -> tuple(database, metadata) }
            .first()

        PREDICT_VICAT_ORFS(NORMALIZE_FASTA.out.normalized_records)

        PREDICT_VICAT_ORFS.out.orfs.view {
            prefix, type, proteins, orfMap, headerMap ->
                "VICAT_ORFS sample=${prefix} type=${type} proteins=${proteins.name}"
        }

        RUN_VICAT_DIAMOND(
            PREDICT_VICAT_ORFS.out.orfs,
            ch_vicat_database
        )

        RUN_VICAT_DIAMOND.out.results.view {
            prefix, type, orfMap, headerMap, diamond ->
                "VICAT_DIAMOND sample=${prefix} type=${type} alignments=${diamond.name}"
        }

        STANDARDIZE_VICAT(
            RUN_VICAT_DIAMOND.out.results,
            ch_vicat_database
        )

        STANDARDIZE_VICAT.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_VICAT.out.evidence
        )
    }

    // Prepare or validate the persistent VITAP database now. Taxonomic
    // assignment runs later against the post-gate, provirus-refined FASTA.
    if( runVitap ) {
        def userSuppliedDatabase = params.vitap_db != null
        if( userSuppliedDatabase && params.vitap_vmr != null ) {
            error 'Use either --vitap_db PATH or --vitap_vmr PATH, not both.'
        }

        def vitapDatabasePath = userSuppliedDatabase
            ? file(params.vitap_db).toString()
            : file(params.vitap_dir).toString()
        def vitapDatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'
        def vitapVmrSource = params.vitap_vmr == null
            ? ''
            : file(params.vitap_vmr).toString()
        def vitapLabel = params.vitap_db_label == null
            ? ''
            : params.vitap_db_label.toString()
        def vitapAutoDownload = parseBooleanParameter(
            params.vitap_auto_download,
            '--vitap_auto_download'
        )
        def vitapUpdateDatabase = parseBooleanParameter(
            params.vitap_update_database,
            '--vitap_update_database'
        )
        def vitapCleanupSource = parseBooleanParameter(
            params.vitap_cleanup_source,
            '--vitap_cleanup_source'
        )

        ch_vitap_database_request = Channel.of(
            tuple(
                vitapDatabasePath,
                vitapDatabaseSource,
                vitapVmrSource,
                vitapLabel,
                vitapAutoDownload,
                vitapUpdateDatabase,
                vitapCleanupSource,
                file("${projectDir}/bin/prepare_vitap_vmr.py"),
                file("${projectDir}/bin/run_vitap_update.py")
            )
        )

        PREPARE_VITAP_DATABASE(ch_vitap_database_request)

        PREPARE_VITAP_DATABASE.out.database.view { database, metadata ->
            "VITAP_DB database=${database} metadata=${metadata.name}"
        }

        ch_vitap_database = PREPARE_VITAP_DATABASE.out.database
    }

    // Prepare or validate the persistent vConTACT3 reference database. The
    // future analysis module will consume the post-gate, provirus-refined
    // FASTA alongside this version-resolved database bundle.
    if( runVcontact3 ) {
        def userSuppliedDatabase = params.vcontact3_db != null
        def vcontact3DatabasePath = userSuppliedDatabase
            ? file(params.vcontact3_db).toString()
            : file(params.vcontact3_dir).toString()
        def vcontact3DatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'
        def vcontact3AutoDownload = parseBooleanParameter(
            params.vcontact3_auto_download,
            '--vcontact3_auto_download'
        )
        def vcontact3UpdateDatabase = parseBooleanParameter(
            params.vcontact3_update_database,
            '--vcontact3_update_database'
        )
        def vcontact3CleanupArchive = parseBooleanParameter(
            params.vcontact3_cleanup_archive,
            '--vcontact3_cleanup_archive'
        )

        ch_vcontact3_database_request = Channel.of(
            tuple(
                vcontact3DatabasePath,
                vcontact3DatabaseSource,
                vcontact3AutoDownload,
                vcontact3UpdateDatabase,
                vcontact3CleanupArchive,
                file("${projectDir}/bin/validate_vcontact3_database.py")
            )
        )

        PREPARE_VCONTACT3_DATABASE(ch_vcontact3_database_request)

        PREPARE_VCONTACT3_DATABASE.out.database.view { database, metadata ->
            "VCONTACT3_DB database=${database} metadata=${metadata.name}"
        }

        ch_vcontact3_database = PREPARE_VCONTACT3_DATABASE.out.database
    }

    ch_discovery_evidence_by_sample = ch_discovery_evidence
        .map { prefix, tool, evidence -> tuple(prefix, evidence) }
        .groupTuple()

    ch_discovery_gate_inputs = NORMALIZE_FASTA.out.normalized_records
        .map { prefix, type, fasta, headerMap ->
            tuple(prefix, type, fasta, headerMap)
        }
        .join(ch_discovery_evidence_by_sample, remainder: true)
        .map { joined ->
            if( joined.size() < 5 || joined[4] == null ) {
                // A path input cannot stage an empty collection in Nextflow
                // 26.04.x. Stage a valid header-only placeholder, while the
                // explicit count tells the module not to pass it to Python.
                return tuple(
                    joined[0],
                    joined[1],
                    joined[2],
                    joined[3],
                    0,
                    file("${projectDir}/assets/empty_discovery_evidence.tsv")
                )
            }
            tuple(
                joined[0],
                joined[1],
                joined[2],
                joined[3],
                joined[4].size(),
                joined[4]
            )
        }

    DISCOVERY_GATE(ch_discovery_gate_inputs)

    DISCOVERY_GATE.out.candidates.view {
        prefix, type, candidates, audit, summary ->
            "DISCOVERY_GATE sample=${prefix} type=${type} candidates=${candidates.name} audit=${audit.name}"
    }

    if( runCheckv ) {
        def userSuppliedDatabase = params.checkv_db != null
        def checkvDatabasePath = userSuppliedDatabase
            ? file(params.checkv_db).toString()
            : file("${params.dbdir}/checkv/checkv_db").toString()
        def checkvDatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'
        def checkvAutoDownload = parseBooleanParameter(
            params.checkv_auto_download,
            '--checkv_auto_download'
        )

        ch_checkv_database_request = Channel.of(
            tuple(
                checkvDatabasePath,
                checkvDatabaseSource,
                checkvAutoDownload
            )
        )

        PREPARE_CHECKV_DATABASE(ch_checkv_database_request)

        PREPARE_CHECKV_DATABASE.out.database.view { database, metadata ->
            "CHECKV_DB database=${database} metadata=${metadata.name}"
        }

        ch_checkv_database = PREPARE_CHECKV_DATABASE.out.database

        RUN_CHECKV(
            DISCOVERY_GATE.out.candidates,
            ch_checkv_database
        )

        RUN_CHECKV.out.results.view {
            prefix, type, quality, completeness, contamination, completeGenomes,
            proviruses, log, metadata ->
                "CHECKV sample=${prefix} type=${type} quality=${quality.name} metadata=${metadata.name}"
        }

        ch_checkv_for_standardizer = RUN_CHECKV.out.results.map {
            prefix, type, quality, completeness, contamination, completeGenomes,
            proviruses, log, metadata ->
                tuple(
                    prefix,
                    type,
                    quality,
                    completeness,
                    contamination,
                    completeGenomes,
                    metadata
                )
        }

        ch_checkv_standardizer_input = ch_checkv_for_standardizer.join(
            DISCOVERY_GATE.out.candidates.map {
                prefix, type, candidates, audit, summary -> tuple(prefix, candidates)
            }
        )

        STANDARDIZE_CHECKV(ch_checkv_standardizer_input)

        STANDARDIZE_CHECKV.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_provirus_evidence = ch_provirus_evidence.mix(
            STANDARDIZE_CHECKV.out.evidence
        )
    }

    ch_provirus_evidence_by_sample = ch_provirus_evidence
        .map { prefix, tool, evidence -> tuple(prefix, evidence) }
        .groupTuple()

    ch_provirus_refinement_inputs = DISCOVERY_GATE.out.candidates
        .join(ch_provirus_evidence_by_sample, remainder: true)
        .map { joined ->
            if( joined.size() < 6 || joined[5] == null ) {
                return tuple(
                    joined[0],
                    joined[1],
                    joined[2],
                    joined[3],
                    joined[4],
                    0,
                    file("${projectDir}/assets/empty_discovery_evidence.tsv"),
                    allowCt3OnlyRefinement
                )
            }
            tuple(
                joined[0],
                joined[1],
                joined[2],
                joined[3],
                joined[4],
                joined[5].size(),
                joined[5],
                allowCt3OnlyRefinement
            )
        }

    REFINE_PROVIRAL_REGIONS(ch_provirus_refinement_inputs)

    REFINE_PROVIRAL_REGIONS.out.refined.view {
        prefix, type, refinedFasta, regionMap, boundaryAudit, summary ->
            "REFINED sample=${prefix} type=${type} fasta=${refinedFasta.name} map=${regionMap.name}"
    }

    if( runVitap ) {
        RUN_VITAP(
            REFINE_PROVIRAL_REGIONS.out.refined,
            ch_vitap_database,
            vitapIncludeLowConfidence
        )

        RUN_VITAP.out.results.view {
            prefix, type, regionMap, best, allLineages, fallback, log, metadata ->
                "VITAP sample=${prefix} type=${type} best=${best.name} all=${allLineages.name}"
        }

        STANDARDIZE_VITAP(
            RUN_VITAP.out.results,
            ch_vitap_database
        )

        STANDARDIZE_VITAP.out.evidence.view { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
    }
}
