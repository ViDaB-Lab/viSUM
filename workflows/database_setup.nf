// Preparation-only workflow. Uses the same process modules and parameter contracts as analysis.
include { PREPARE_GENOMAD_DATABASE } from '../modules/local/genomad_database'
include { PREPARE_VIRSORTER2_DATABASE } from '../modules/local/virsorter2_database'
include { PREPARE_CENOTETAKER3_DATABASE } from '../modules/local/cenotetaker3_database'
include { PREPARE_DEEP6_DATABASE } from '../modules/local/deep6_database'
include { PREPARE_DEEPMICROCLASS2 } from '../modules/local/deepmicroclass2_installation'
include { PREPARE_VIRBOT_DATABASE } from '../modules/local/virbot_installation'
include { PREPARE_GIANTHUNTER_DATABASE } from '../modules/local/gianthunter_database'
include { PREPARE_VICAT_DATABASE } from '../modules/local/vicat_database'
include { PREPARE_CHECKV_DATABASE } from '../modules/local/checkv_database'
include { PREPARE_VITAP_DATABASE } from '../modules/local/vitap_database'
include { PREPARE_VCONTACT3_DATABASE } from '../modules/local/vcontact3_database'
include { PREPARE_VICAT_NONVIRAL_DATABASE } from '../modules/local/vicat_nonviral_database'

def parseBooleanParameter(value, name) {
    def normalized = value?.toString()?.trim()?.toLowerCase()
    if( !(normalized in ['true', 'false']) ) { error "${name} must be true or false." }
    return normalized == 'true'
}
def parsePositiveIntegerParameter(value, name) {
    if( !(value?.toString() ==~ /[1-9][0-9]*/) ) { error "${name} must be a positive integer." }
    return value.toString().toInteger()
}

workflow DATABASE_SETUP {
    def flags = [:]
    ['genomad','virsorter2','cenotetaker3','deep6','deepmicroclass2','virbot','gianthunter','vicat','checkv','vitap','vcontact3'].each { tool ->
        flags[tool] = parseBooleanParameter(params["run_${tool}"], "--run_${tool}")
    }
    if( !flags.values().any { it } ) { error 'Setup requires at least one enabled database-bearing tool.' }
    if( parseBooleanParameter(params.run_tesorter, '--run_tesorter') ) {
        log.warn('TEsorter has no standalone database preparation process; its environment/reference runtime must be checked in the analysis smoke test.')
    }
    println "viSUM setup-only selection: ${flags.findAll { key, enabled -> enabled }.keySet().join(', ')}"
    // A completion token serializes installations across different process types.
    ready = Channel.value(true)
    if( flags.genomad ) {
        def userSuppliedDatabase = params.genomad_db != null
        def genomadDatabasePath = userSuppliedDatabase
            ? file(params.genomad_db).toString()
            : file("${params.dbdir}/genomad/genomad_db").toString()
        def genomadDatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_genomad_database_request = ready.map { ignored ->
            tuple(
                genomadDatabasePath,
                genomadDatabaseSource,
                params.genomad_auto_download
            )
        }

        PREPARE_GENOMAD_DATABASE(ch_genomad_database_request)

        ready = PREPARE_GENOMAD_DATABASE.out.database.map { result -> true }.first()
    }
    if( flags.virsorter2 ) {
        def userSuppliedDatabase = params.virsorter2_db != null
        def virsorter2DatabasePath = userSuppliedDatabase
            ? file(params.virsorter2_db).toString()
            : file("${params.dbdir}/virsorter2/db").toString()
        def virsorter2DatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_virsorter2_database_request = ready.map { ignored ->
            tuple(
                virsorter2DatabasePath,
                virsorter2DatabaseSource,
                params.virsorter2_auto_download
            )
        }

        PREPARE_VIRSORTER2_DATABASE(ch_virsorter2_database_request)

        ready = PREPARE_VIRSORTER2_DATABASE.out.database.map { result -> true }.first()
    }
    if( flags.cenotetaker3 ) {
        def userSuppliedDatabase = params.ct3_db != null
        def cenotetaker3DatabasePath = userSuppliedDatabase
            ? file(params.ct3_db).toString()
            : file("${params.dbdir}/cenotetaker3/ct3_DBs").toString()
        def cenotetaker3DatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        ch_cenotetaker3_database_request = ready.map { ignored ->
            tuple(
                cenotetaker3DatabasePath,
                cenotetaker3DatabaseSource,
                params.ct3_auto_download,
                params.ct3_hmm_db_version
            )
        }

        PREPARE_CENOTETAKER3_DATABASE(ch_cenotetaker3_database_request)

        ready = PREPARE_CENOTETAKER3_DATABASE.out.database.map { result -> true }.first()
    }
    if( flags.deep6 ) {
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


        ch_deep6_database_request = ready.map { ignored ->
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

        ready = PREPARE_DEEP6_DATABASE.out.database.map { result -> true }.first()
    }
    if( flags.deepmicroclass2 ) {
        def userSuppliedInstallation = params.deepmicroclass2_dir != null
        def deepmicroclass2InstallationPath = userSuppliedInstallation
            ? file(params.deepmicroclass2_dir).toString()
            : file("${params.tooldir}/deepmicroclass2/DeepMicroClass2").toString()
        def deepmicroclass2InstallationSource = userSuppliedInstallation
            ? 'user-supplied'
            : 'viSUM-managed'


        ch_deepmicroclass2_installation_request = ready.map { ignored ->
                tuple(
                    deepmicroclass2InstallationPath,
                    deepmicroclass2InstallationSource,
                    params.deepmicroclass2_auto_download,
                    params.deepmicroclass2_repository,
                    params.deepmicroclass2_revision
                )
            }

        PREPARE_DEEPMICROCLASS2(ch_deepmicroclass2_installation_request)

        ready = PREPARE_DEEPMICROCLASS2.out.installation.map { result -> true }.first()
    }
    if( flags.virbot ) {
        def userSuppliedInstallation = params.virbot_dir != null
        def virbotInstallationPath = userSuppliedInstallation
            ? file(params.virbot_dir).toString()
            : file("${params.dbdir}/virbot/VirBot").toString()
        def virbotInstallationSource = userSuppliedInstallation
            ? 'user-supplied'
            : 'viSUM-managed'


        ch_virbot_installation_request = ready.map { ignored ->
                tuple(
                    virbotInstallationPath,
                    virbotInstallationSource,
                    params.virbot_auto_download,
                    params.virbot_repository,
                    params.virbot_revision
                )
            }

        PREPARE_VIRBOT_DATABASE(ch_virbot_installation_request)

        ready = PREPARE_VIRBOT_DATABASE.out.installation.map { result -> true }.first()
    }
    if( flags.gianthunter ) {
        def userSuppliedDatabase = params.gianthunter_db != null
        def gianthunterDatabasePath = userSuppliedDatabase
            ? file(params.gianthunter_db).toString()
            : file("${params.dbdir}/gianthunter/gianthunter_db_v1").toString()
        def gianthunterDatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'


        ch_gianthunter_database_request = ready.map { ignored ->
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

        ready = PREPARE_GIANTHUNTER_DATABASE.out.database.map { result -> true }.first()
    }
    if( flags.vicat ) {
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

        def userSuppliedDatabase = params.vicat_db != null
        if( userSuppliedDatabase && vicatSourceProteins ) {
            error 'Choose --vicat_db OR the two MetaVR source files, not both.'
        }
        def vicatBaseDatabase = userSuppliedDatabase
            ? file(params.vicat_db).toString()
            : file(params.vicat_managed_db).toString()
        def vicatDatabasePath = vicatBaseDatabase
        def vicatDatabaseSource = userSuppliedDatabase
            ? 'user-supplied'
            : 'viSUM-managed'

        // Resolve both sides before starting a potentially expensive viral build.
        def vicatNonviralDatabasePath = file(params.vicat_nonviral_db).toString()
        def nonviralRaw = [params.vicat_nonviral_manifest, params.vicat_nonviral_package_root, params.vicat_nonviral_metadata_root]
        if( nonviralRaw.any { it != null } && !nonviralRaw.every { it != null } ) {
            error 'A raw nonviral build requires manifest, package_root, and metadata_root together.'
        }
        if( params.vicat_nonviral_classified_dir != null && nonviralRaw.any { it != null } ) {
            error 'Choose --vicat_nonviral_classified_dir OR the three raw NCBI source inputs, not both.'
        }
        def nonviralSources = ([params.vicat_nonviral_classified_dir] + nonviralRaw).collect {
            it == null ? '' : file(it).toString()
        }
        nonviralSources.findAll { it }.each { source ->
            if( !file(source).exists() ) { error "Nonviral source input does not exist: ${source}" }
        }
        if( !file(vicatNonviralDatabasePath).exists() && !nonviralSources.any { it } ) {
            error 'No viCAT nonviral database or build inputs supplied. Set --vicat_nonviral_db to a completed database, or supply classified/raw NCBI inputs.'
        }
        parsePositiveIntegerParameter(params.vicat_nonviral_build_cpus, '--vicat_nonviral_build_cpus')
        def nonviralHelpers = ['prepare_vicat_nonviral_build.py', 'classify_vicat_nonviral_references.py',
            'build_vicat_nonviral_database.sh', 'prepare_vicat_nonviral_cluster_metadata.py',
            'parse_diamond_dbinfo.py'].collect { file("${projectDir}/bin/${it}") }

        ch_vicat_database_request = ready.map { ignored ->
                tuple(
                    vicatDatabasePath,
                    vicatDatabaseSource,
                    vicatBaseDatabase,
                    '',
                    vicatSourceProteins,
                    vicatSourceMetadata,
                    params.vicat_metavr_proteins_sha256,
                    params.vicat_metavr_metadata_sha256,
                    params.vicat_expected_uvigs,
                    params.vicat_expected_proteins,
                    params.vicat_expected_representatives,
                    '',
                    '',
                    0.90,
                    0.80,
                    50
                )
            }

        PREPARE_VICAT_DATABASE(ch_vicat_database_request)



        // Serialize the two heavyweight builds; their memory requests add up.
        ch_vicat_nonviral_database_request = PREPARE_VICAT_DATABASE.out.database
            .take(1)
            .map { viralDatabase, viralSetupMetadata ->
                tuple(vicatNonviralDatabasePath, 'class-aware-nonviral',
                    nonviralSources[0], nonviralSources[1], nonviralSources[2], nonviralSources[3], nonviralHelpers)
            }

        PREPARE_VICAT_NONVIRAL_DATABASE(ch_vicat_nonviral_database_request)




        ready = PREPARE_VICAT_NONVIRAL_DATABASE.out.database.map { result -> true }.first()
    }
    if( flags.checkv ) {
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

        ch_checkv_database_request = ready.map { ignored ->
            tuple(
                checkvDatabasePath,
                checkvDatabaseSource,
                checkvAutoDownload
            )
        }

        PREPARE_CHECKV_DATABASE(ch_checkv_database_request)

        ready = PREPARE_CHECKV_DATABASE.out.database.map { result -> true }.first()
    }
    if( flags.vitap ) {
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

        ch_vitap_database_request = ready.map { ignored ->
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
        }

        PREPARE_VITAP_DATABASE(ch_vitap_database_request)

        ready = PREPARE_VITAP_DATABASE.out.database.map { result -> true }.first()
    }
    if( flags.vcontact3 ) {
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

        ch_vcontact3_database_request = ready.map { ignored ->
            tuple(
                vcontact3DatabasePath,
                vcontact3DatabaseSource,
                vcontact3AutoDownload,
                vcontact3UpdateDatabase,
                vcontact3CleanupArchive
            )
        }

        PREPARE_VCONTACT3_DATABASE(ch_vcontact3_database_request)

        ready = PREPARE_VCONTACT3_DATABASE.out.database.map { result -> true }.first()
    }
    ready.view { 'viSUM setup complete: all selected preparation processes passed. See outdir/database_setup for per-tool provenance.' }
}
