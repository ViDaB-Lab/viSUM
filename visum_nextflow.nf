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
include { PREPARE_VICAT_NONVIRAL_DATABASE } from './modules/local/vicat_nonviral_database'
include { PREDICT_VICAT_ORFS } from './modules/local/predict_vicat_orfs'
include { RUN_VICAT_DIAMOND } from './modules/local/run_vicat_diamond'
include { RUN_VICAT_NONVIRAL_DIAMOND } from './modules/local/run_vicat_nonviral_diamond'
include { PREPARE_VICAT_REFERENCE_SUBSET } from './modules/local/prepare_vicat_reference_subset'
include { PREPARE_VICAT_HITS } from './modules/local/prepare_vicat_hits'
include { STANDARDIZE_VICAT } from './modules/local/standardize_vicat'
include { PROJECT_VICAT_REFINED } from './modules/local/project_vicat_refined'
include { DISCOVERY_GATE } from './modules/local/discovery_gate'
include { PREPARE_CHECKV_DATABASE } from './modules/local/checkv_database'
include { RUN_CHECKV } from './modules/local/run_checkv'
include { STANDARDIZE_CHECKV } from './modules/local/standardize_checkv'
include { REFINE_PROVIRAL_REGIONS } from './modules/local/refine_proviral_regions'
include { RUN_TESORTER } from './modules/local/run_tesorter'
include { STANDARDIZE_TESORTER } from './modules/local/standardize_tesorter'
include { PREPARE_VITAP_DATABASE } from './modules/local/vitap_database'
include { RUN_VITAP } from './modules/local/run_vitap'
include { STANDARDIZE_VITAP } from './modules/local/standardize_vitap'
include { PREPARE_VCONTACT3_DATABASE } from './modules/local/vcontact3_database'
include { RUN_VCONTACT3 } from './modules/local/run_vcontact3'
include { STANDARDIZE_VCONTACT3 } from './modules/local/standardize_vcontact3'
include { VIHARMONY } from './modules/local/viharmony'


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


def parseNonnegativeIntegerParameter(rawValue, parameterName) {
    Integer parsedValue
    try {
        parsedValue = rawValue as Integer
    }
    catch( Exception ignored ) {
        throw new IllegalArgumentException(
            "Invalid value '${rawValue}' for ${parameterName}. Use a nonnegative integer."
        )
    }

    if( parsedValue < 0 ) {
        throw new IllegalArgumentException(
            "Invalid value '${rawValue}' for ${parameterName}. Use a nonnegative integer."
        )
    }
    return parsedValue
}


def viewChannel(enabled, channel, formatter) {
    if( enabled ) {
        channel.view(formatter)
    }
}


def validateEvidenceArtifacts(prefix, stage, expectedTools, actualTools, evidenceFiles) {
    def expected = expectedTools.collect { it.toString() }.toSet()
    def actual = actualTools.collect { it.toString() }.toSet()
    def missing = expected.findAll { !actual.contains(it) }.sort()
    def unexpected = actual.findAll { !expected.contains(it) }.sort()

    if( actualTools.size() != evidenceFiles.size() ) {
        throw new IllegalStateException(
            "Evidence channel contract failed for sample '${prefix}' at ${stage}: " +
            "received ${actualTools.size()} tool labels but ${evidenceFiles.size()} files."
        )
    }
    if( missing ) {
        throw new IllegalStateException(
            "Missing required evidence for sample '${prefix}' at ${stage}: " +
            missing.join(', ') + '. A valid zero-call run must still emit a header-only evidence file.'
        )
    }
    if( unexpected ) {
        throw new IllegalStateException(
            "Unexpected evidence for sample '${prefix}' at ${stage}: " +
            unexpected.join(', ') + '.'
        )
    }
}


def visumHelp() {
    return '''
viSUM 0.1.0 — viral sequence discovery, refinement, and taxonomy harmonization

USAGE
  Single FASTA (recommended local launcher):
    ./visum -c visum.config \\
      --input INPUT.fasta --prefix SAMPLE --type dna|rna [options]

  Multiple FASTAs:
    ./visum -c visum.config \\
      --prefix_many samples.csv --indir INPUT_DIRECTORY [options]

  Help only:
    ./visum -c visum.config --help

  Database preparation only (enabled --run_* tools, no samples):
    ./visum -c visum.config --setup

INPUT
  --setup                     Prepare enabled tools only, without analysis inputs.
                               Selection uses the same --run_* flags as analysis.
  --input PATH                 Input FASTA for a single sample.
  --prefix TEXT                Sample identifier for single-sample mode.
                               Allowed: letters, numbers, ., _, and -.
  --type dna|rna               Input type. RNA mode supports mixed DNA/RNA
                               viruses in metatranscriptomic assemblies.
  --prefix_many PATH           CSV with header: prefix,type,fasta.
  --indir PATH                 Base directory for FASTAs in --prefix_many.

OUTPUT AND RESOURCES
  --outdir PATH                Results directory [results].
  --dbdir PATH                 Persistent managed databases [databases].
  --tooldir PATH               Persistent managed tool bundles [tools].
  --max_cpus INT               Aggregate local CPU scheduler cap through
                               ./visum; per-task ceiling otherwise [8].
  --max_memory MEMORY          Per-task memory ceiling [48 GB].
  --harmonizer_audit MODE      none, compact, or full [compact].
  --rna_pair_homology_floor BOOL  Experimental RNA Deep6/viCAT-only homology floor [false].
  --show_channel_messages BOOL Print detailed channel emissions [false].
  --ictv_csv PATH              Canonical ICTV rank table
                               [assets/ICTV_VMR_MSL41.csv].

PROGRAM SELECTION
  --run_genomad BOOL           geNomad discovery [true].
  --run_virsorter2 BOOL        VirSorter2 discovery [true].
  --run_cenotetaker3 BOOL      Cenote-Taker 3 discovery [true].
  --run_deep6 BOOL             Deep6; runs only for RNA inputs [true].
  --run_deepmicroclass2 BOOL   DeepMicroClass2; DNA inputs only [true].
  --run_virbot BOOL            VirBot; RNA inputs only [true].
  --run_gianthunter BOOL       GiantHunter; DNA inputs only [true].
  --run_vicat BOOL             viCAT discovery/taxonomy [false].
  --run_checkv BOOL            CheckV candidate quality/refinement [true].
  --run_tesorter BOOL          TEsorter retroelement evidence [true].
  --run_vitap BOOL             VITAP taxonomy refinement [true].
  --run_vcontact3 BOOL         vConTACT3 taxonomy refinement [true].

COMMON ANALYSIS OPTIONS
  --genomad_score_calibration BOOL  Calibrate geNomad scores when supported
                                    by sample size [true].
  --genomad_splits INT              0 is fastest; increase to lower memory [0].
  --vs2_min_length INT              VirSorter2 minimum length [1500].
  --vs2_min_score FLOAT             VirSorter2 minimum score [0.5].
  --ct3_minlen_circ INT             CT3 circular minimum length [1000].
  --ct3_circ_minhall INT            CT3 circular hallmark minimum [1].
  --ct3_minlen_linear INT           CT3 linear minimum length [1000].
  --ct3_linear_minhall INT          CT3 linear hallmark minimum [1].
  --deep6_minlen INT                Deep6 minimum sequence length [250].
  --deep6_min_score FLOAT           Deep6 minimum winning score [0.7].
  --deep6_median_multiplier FLOAT   Required winner/median ratio [1.25].
  --deepmicroclass2_model MODEL     8class, high_precision, or 300bp [8class].
  --virbot_sensitive BOOL           Add VirBot DIAMOND search [false].
  --virbot_taxa TOP|LCA             VirBot taxonomy method [TOP].
  --gianthunter_min_length INT      GiantHunter minimum length [3000].
  --gianthunter_reject FLOAT        Minimum aligned-protein fraction [0.1].
  --gianthunter_query_cover INT     Minimum query coverage percent [40].
  --allow_ct3_only_refinement BOOL  Permit CT3-only provirus trimming [false; opt-in].
  --vicat_provirus_min_overlap FLOAT  viCAT cluster overlap needed to corroborate CheckV [0.5].
  --vitap_include_low_confidence BOOL
                                    Retain VITAP low-confidence calls [false].
  --vcontact3_db_domain MODE        both, prokaryotes, or eukaryotes [both].
  --vcontact3_min_taxonomy_length INT
                                    Shorter calls remain audit-only [1000].

viCAT OPTIONS
  --vicat_diamond_sensitivity MODE  faster through ultra-sensitive [sensitive].
  --vicat_min_bitscore FLOAT        Minimum DIAMOND bit score [50].
  --vicat_min_query_cover INT       Minimum query coverage percent [30].
  --vicat_top_percent FLOAT         Top-hit score window [5].
  --vicat_block_size FLOAT          DIAMOND block size [2.0].
  --vicat_index_chunks INT          DIAMOND index chunks [4].
  --vicat_locus_overlap FLOAT       GV/RV ORF overlap threshold [0.80].
  --vicat_competitive_min_margin FLOAT
                                    Minimum viral/cellular bitscore margin [0.05].
  --vicat_cluster_min_viral_loci auto|INT
                                    Clean-sequence threshold [auto: DNA=2, RNA=1].
                                    Cellular context always requires >=2 loci.
  --vicat_dna_single_locus_rescue off|strict
                                    Guarded clean single-locus DNA discovery [strict].
  --vicat_cluster_max_neutral_gap INT
                                    Ambiguous/uninformative loci bridged in a cluster [1].
  --vicat_orf_taxonomy_support FLOAT
                                    Within-ORF taxonomy support [0.60].
  --vicat_contig_taxonomy_support FLOAT
                                    Across-ORF taxonomy support [0.60].

DATABASE AND INSTALLATION OVERRIDES
  --genomad_db PATH            Existing geNomad database.
  --virsorter2_db PATH         Existing VirSorter2 database.
  --ct3_db PATH                Existing Cenote-Taker 3 database.
  --deep6_dir PATH             Existing Deep6 installation.
  --deep6_model PATH           Existing Deep6 model directory.
  --deepmicroclass2_dir PATH   Existing DeepMicroClass2 installation.
  --virbot_dir PATH            Existing VirBot installation.
  --gianthunter_db PATH        Existing GiantHunter database.
  --checkv_db PATH             Existing CheckV database.
  --vicat_db PATH              Existing viCAT database.
  --vicat_nonviral_db PATH     Existing nonviral database or local build destination.
  --vicat_nonviral_classified_dir PATH  Classified inputs for a local nonviral build.
  --vicat_nonviral_manifest PATH       NCBI source assembly manifest (alternative to classified inputs).
  --vicat_nonviral_package_root PATH   Extracted NCBI protein package root.
  --vicat_nonviral_metadata_root PATH  NCBI feature-table root.
  --vicat_metavr_proteins PATH MetaVR proteins for a local viCAT build.
  --vicat_metavr_metadata PATH MetaVR metadata for a local viCAT build.
  --vitap_db PATH              Existing VITAP database.
  --vitap_vmr PATH             VMR workbook/CSV for a VITAP build.
  --vitap_update_database BOOL Check for a newer VMR [false].
  --vcontact3_db PATH          Existing vConTACT3 database or manifest.
  --vcontact3_update_database BOOL
                               Resolve and install the latest release [false].

PER-TOOL RESOURCE OVERRIDES
  --genomad_cpus INT           --virsorter2_cpus INT
  --ct3_cpus INT               --deep6_cpus INT
  --deepmicroclass2_cpus INT   --virbot_cpus INT
  --gianthunter_cpus INT       --vicat_orf_cpus INT
  --vicat_cpus INT             --vicat_standardizer_cpus INT
  --checkv_cpus INT
  --tesorter_cpus INT
  --vitap_cpus INT             --vcontact3_cpus INT

  Each tool also has a corresponding --*_memory and --*_time option in
  visum.config. Tool CPU requests are capped by --max_cpus.

NOTES
  BOOL values must be true or false. Use -resume to reuse completed work.
  Advanced version pins, download URLs, checksums, and database-build controls
  are documented in visum.config and normally should not be changed.
'''.stripIndent().trim()
}


include { DATABASE_SETUP } from './workflows/database_setup'

workflow {

    def helpRequested = params.help == true ||
        params.help?.toString()?.trim()?.toLowerCase() == 'true'
    if( helpRequested ) {
        println visumHelp()
        return
    }

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
    def launcherCpuCap = System.getenv('VISUM_LOCAL_CPU_CAP')
    def cpuCapMode = System.getenv('VISUM_CPU_CAP_MODE')
    if( launcherCpuCap != null ) {
        def enforcedCpus = parsePositiveIntegerParameter(launcherCpuCap, 'VISUM_LOCAL_CPU_CAP')
        if( enforcedCpus != maxCpus ) {
            error "Launcher CPU cap (${enforcedCpus}) does not match --max_cpus (${maxCpus})."
        }
    }
    else {
        log.warn("Direct Nextflow invocation detected: --max_cpus=${maxCpus} caps each task but not aggregate local concurrency. Use ./visum to enforce the aggregate CPU scheduler cap.")
    }
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
        '--vicat_standardizer_cpus': params.vicat_standardizer_cpus,
        '--checkv_cpus': params.checkv_cpus,
        '--tesorter_cpus': params.tesorter_cpus,
        '--vitap_cpus': params.vitap_cpus,
        '--vcontact3_cpus': params.vcontact3_cpus,
    ]
    analysisCpuParameters.each { parameterName, rawValue ->
        parsePositiveIntegerParameter(rawValue, parameterName)
    }
    parseNonnegativeIntegerParameter(
        params.ct3_linear_minhall,
        '--ct3_linear_minhall'
    )
    parseNonnegativeIntegerParameter(
        params.ct3_circ_minhall,
        '--ct3_circ_minhall'
    )

    println(
        "viSUM resource controls: max_cpus=${maxCpus} " +
        "(${launcherCpuCap != null ? (cpuCapMode ?: 'aggregate scheduler cap enforced') : 'per-task ceiling only'}), " +
        "max_memory=${params.max_memory} (per-task ceiling); " +
        "preferred tool CPUs (each capped at max_cpus)=" + analysisCpuParameters.collect { name, value ->
            "${name.substring(2)}=${value}"
        }.join(', ')
    )

    if( parseBooleanParameter(params.setup, '--setup') ) {
        if( params.input != null || params.prefix_many != null || params.indir != null || params.prefix != null || params.type != null ) {
            error '--setup does not accept analysis input parameters. Remove input/prefix/type/prefix_many/indir.'
        }
        DATABASE_SETUP()
        return
    }

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
    parseBooleanParameter(
        params.genomad_score_calibration,
        '--genomad_score_calibration'
    )
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
    def runTesorter = parseBooleanParameter(params.run_tesorter, '--run_tesorter')
    def runVitap = parseBooleanParameter(params.run_vitap, '--run_vitap')
    def runVcontact3 = parseBooleanParameter(
        params.run_vcontact3,
        '--run_vcontact3'
    )
    def showChannelMessages = parseBooleanParameter(
        params.show_channel_messages,
        '--show_channel_messages'
    )

    // These contracts distinguish a successful zero-call result (a
    // header-only evidence table) from a missing process result. Tool sets
    // differ by molecule type only where the discovery programs do.
    def sharedDiscoveryTools = []
    if( runGenomad ) sharedDiscoveryTools << 'genomad'
    if( runVirsorter2 ) sharedDiscoveryTools << 'virsorter2'
    if( runCenotetaker3 ) sharedDiscoveryTools << 'cenotetaker3'
    if( runVicat ) sharedDiscoveryTools << 'vicat'

    def dnaDiscoveryTools = []
    if( runDeepmicroclass2 ) dnaDiscoveryTools << 'deepmicroclass2'
    if( runGianthunter ) dnaDiscoveryTools << 'gianthunter'

    def rnaDiscoveryTools = []
    if( runDeep6 ) rnaDiscoveryTools << 'deep6'
    if( runVirbot ) rnaDiscoveryTools << 'virbot'

    def expectedDiscoveryToolsByType = [
        dna: (sharedDiscoveryTools + dnaDiscoveryTools).unique().sort(),
        rna: (sharedDiscoveryTools + rnaDiscoveryTools).unique().sort(),
    ]

    def sharedProvirusTools = []
    if( runGenomad ) sharedProvirusTools << 'genomad'
    if( runVirsorter2 ) sharedProvirusTools << 'virsorter2'
    if( runCenotetaker3 ) sharedProvirusTools << 'cenotetaker3'
    if( runVicat ) {
        sharedProvirusTools << 'vicat'
        sharedProvirusTools << 'vicat_boundary'
    }
    if( runCheckv ) sharedProvirusTools << 'checkv'

    def rnaPairHomologyFloor = parseBooleanParameter(params.rna_pair_homology_floor, '--rna_pair_homology_floor')
    def harmonyOnlyTools = []
    if( runVicat ) harmonyOnlyTools << 'vicat_context'
    if( runCheckv ) harmonyOnlyTools << 'checkv'
    if( runTesorter ) harmonyOnlyTools << 'tesorter'
    if( runVitap ) harmonyOnlyTools << 'vitap'
    if( runVcontact3 ) harmonyOnlyTools << 'vcontact3'

    def expectedHarmonyToolsByType = [
        dna: (expectedDiscoveryToolsByType.dna + harmonyOnlyTools).unique().sort(),
        rna: (expectedDiscoveryToolsByType.rna + harmonyOnlyTools +
            (runVicat ? ['vicat_loci'] : [])).unique().sort(),
    ]
    def vcontact3DbDomain = params.vcontact3_db_domain
        ?.toString()
        ?.trim()
        ?.toLowerCase()
    if( runVcontact3 && !(vcontact3DbDomain in ['both', 'prokaryotes', 'eukaryotes']) ) {
        error "Invalid --vcontact3_db_domain '${params.vcontact3_db_domain}'. Use both, prokaryotes, or eukaryotes."
    }
    def vitapIncludeLowConfidence = parseBooleanParameter(
        params.vitap_include_low_confidence,
        '--vitap_include_low_confidence'
    )
    def allowCt3OnlyRefinement = parseBooleanParameter(
        params.allow_ct3_only_refinement,
        '--allow_ct3_only_refinement'
    )
    def vicatProvirusMinOverlap = params.vicat_provirus_min_overlap as double
    if( vicatProvirusMinOverlap < 0.0 || vicatProvirusMinOverlap > 1.0 ) {
        error "Invalid --vicat_provirus_min_overlap '${params.vicat_provirus_min_overlap}'. Use a value from 0 to 1."
    }
    def harmonizerAudit = params.harmonizer_audit?.toString()?.trim()?.toLowerCase()
    if( !(harmonizerAudit in ['none', 'compact', 'full']) ) {
        error "Invalid --harmonizer_audit '${params.harmonizer_audit}'. Use none, compact, or full."
    }
    def vcontact3MinTaxonomyLength = params.vcontact3_min_taxonomy_length as int
    if( vcontact3MinTaxonomyLength < 0 ) {
        error "Invalid --vcontact3_min_taxonomy_length '${params.vcontact3_min_taxonomy_length}'. Use zero or a positive integer."
    }
    def ictvCsv = file(params.ictv_csv)
    if( !ictvCsv.exists() ) {
        error "Canonical ICTV taxonomy file was not found: ${ictvCsv}"
    }

    // Every standardizer emits sparse, threshold-qualified evidence. These
    // channels are merged and grouped by sample for the discovery gate.
    ch_discovery_evidence = Channel.empty()
    // viHARMONY consumes every standardized evidence table, including tools
    // that do not participate in the first discovery gate.
    ch_harmony_evidence = Channel.empty()
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
        "TEsorter=${runTesorter}, " +
        "VITAP=${runVitap}, " +
        "vConTACT3=${runVcontact3}"
    )

    NORMALIZE_FASTA(ch_samples)

    viewChannel(showChannelMessages, NORMALIZE_FASTA.out.normalized_records) { prefix, type, fasta, headerMap ->
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

        viewChannel(showChannelMessages, PREPARE_GENOMAD_DATABASE.out.database) { database, metadata ->
            "GENOMAD_DB database=${database} metadata=${metadata.name}"
        }

        ch_genomad_database = PREPARE_GENOMAD_DATABASE.out.database
            .map { database, metadata -> database }
            .first()

        RUN_GENOMAD(
            NORMALIZE_FASTA.out.normalized_records,
            ch_genomad_database
        )

        viewChannel(showChannelMessages, RUN_GENOMAD.out.results) { prefix, type, review, virusSummary, virusFasta, virusGenes, virusProteins, plasmidSummary, metadata ->
            "GENOMAD sample=${prefix} type=${type} virus_summary=${virusSummary.name} output=${review.name}"
        }

        ch_genomad_for_standardizer = RUN_GENOMAD.out.results.map {
            prefix, type, review, virusSummary, virusFasta, virusGenes, virusProteins, plasmidSummary, metadata ->
                tuple(prefix, type, virusSummary, virusGenes, plasmidSummary, metadata)
        }

        ch_genomad_standardizer_input = ch_genomad_for_standardizer
            .join(
                NORMALIZE_FASTA.out.normalized_records.map { prefix, type, fasta, headerMap ->
                    tuple(prefix, headerMap)
                }
            )

        STANDARDIZE_GENOMAD(ch_genomad_standardizer_input)

        viewChannel(showChannelMessages, STANDARDIZE_GENOMAD.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_GENOMAD.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_GENOMAD.out.evidence)
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

        viewChannel(showChannelMessages, PREPARE_VIRSORTER2_DATABASE.out.database) { database, metadata ->
            "VIRSORTER2_DB database=${database} metadata=${metadata.name}"
        }

        ch_virsorter2_database = PREPARE_VIRSORTER2_DATABASE.out.database
            .map { database, metadata -> database }
            .first()

        RUN_VIRSORTER2(
            NORMALIZE_FASTA.out.normalized_records,
            ch_virsorter2_database
        )

        viewChannel(showChannelMessages, RUN_VIRSORTER2.out.results) {
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

        viewChannel(showChannelMessages, STANDARDIZE_VIRSORTER2.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_VIRSORTER2.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_VIRSORTER2.out.evidence)
        ch_provirus_evidence = ch_provirus_evidence.mix(
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

        viewChannel(showChannelMessages, PREPARE_CENOTETAKER3_DATABASE.out.database) { database, metadata ->
            "CENOTETAKER3_DB database=${database} metadata=${metadata.name}"
        }

        ch_cenotetaker3_database = PREPARE_CENOTETAKER3_DATABASE.out.database
            .map { database, metadata -> database }
            .first()

        RUN_CENOTETAKER3(
            NORMALIZE_FASTA.out.normalized_records,
            ch_cenotetaker3_database
        )

        viewChannel(showChannelMessages, RUN_CENOTETAKER3.out.results) {
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

        viewChannel(showChannelMessages, STANDARDIZE_CENOTETAKER3.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_CENOTETAKER3.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_CENOTETAKER3.out.evidence)
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

        viewChannel(showChannelMessages, PREPARE_DEEP6_DATABASE.out.database) {
            installation, models, metadata ->
                "DEEP6_DB installation=${installation} models=${models} metadata=${metadata.name}"
        }

        ch_deep6_bundle = PREPARE_DEEP6_DATABASE.out.database
            .map { installation, models, metadata ->
                tuple(installation, models, metadata)
            }
            .first()

        RUN_DEEP6(ch_deep6_samples, ch_deep6_bundle)

        viewChannel(showChannelMessages, RUN_DEEP6.out.results) { prefix, type, scores, logFile, metadata ->
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

        viewChannel(showChannelMessages, STANDARDIZE_DEEP6.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_DEEP6.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_DEEP6.out.evidence)
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

        viewChannel(showChannelMessages, PREPARE_DEEPMICROCLASS2.out.installation) {
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

        viewChannel(showChannelMessages, RUN_DEEPMICROCLASS2.out.results) {
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

        viewChannel(showChannelMessages, STANDARDIZE_DEEPMICROCLASS2.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_DEEPMICROCLASS2.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_DEEPMICROCLASS2.out.evidence)
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

        viewChannel(showChannelMessages, PREPARE_VIRBOT_DATABASE.out.installation) {
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

        viewChannel(showChannelMessages, RUN_VIRBOT.out.results) {
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
            ictvCsv
        )

        viewChannel(showChannelMessages, STANDARDIZE_VIRBOT.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_VIRBOT.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_VIRBOT.out.evidence)
    }

    if( runGianthunter ) {
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

        viewChannel(showChannelMessages, PREPARE_GIANTHUNTER_DATABASE.out.database) { database, metadata ->
            "GIANTHUNTER_DB database=${database} metadata=${metadata.name}"
        }

        ch_gianthunter_bundle = PREPARE_GIANTHUNTER_DATABASE.out.database
            .map { database, metadata -> tuple(database, metadata) }
            .first()

        RUN_GIANTHUNTER(
            ch_gianthunter_samples,
            ch_gianthunter_bundle
        )

        viewChannel(showChannelMessages, RUN_GIANTHUNTER.out.results) {
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
            ictvCsv
        )

        viewChannel(showChannelMessages, STANDARDIZE_GIANTHUNTER.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_GIANTHUNTER.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_GIANTHUNTER.out.evidence)
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
        def vicatCompetitiveMinMargin = params.vicat_competitive_min_margin as Double
        if( !Double.isFinite(vicatCompetitiveMinMargin) ||
            vicatCompetitiveMinMargin < 0.0 || vicatCompetitiveMinMargin >= 1.0 ) {
            error '--vicat_competitive_min_margin must be at least 0 and less than 1.0.'
        }
        def vicatClusterMinViralLoci = params.vicat_cluster_min_viral_loci.toString().trim()
        if( !vicatClusterMinViralLoci.equalsIgnoreCase('auto') &&
            (!vicatClusterMinViralLoci.isInteger() || vicatClusterMinViralLoci.toInteger() < 1) ) {
            error '--vicat_cluster_min_viral_loci must be auto or a positive integer.'
        }
        def vicatClusterMaxNeutralGap = params.vicat_cluster_max_neutral_gap as Integer
        if( vicatClusterMaxNeutralGap < 0 ) {
            error '--vicat_cluster_max_neutral_gap must be zero or greater.'
        }
        def vicatDnaSingleLocusRescue = params.vicat_dna_single_locus_rescue.toString().trim()
        if( !(vicatDnaSingleLocusRescue in ['off', 'strict']) ) {
            error "--vicat_dna_single_locus_rescue must be 'off' or 'strict'."
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

        ch_vicat_database_request = NORMALIZE_FASTA.out.normalized_records
            .take(1)
            .map { prefix, type, fasta, headerMap ->
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

        viewChannel(showChannelMessages, PREPARE_VICAT_DATABASE.out.database) { database, metadata ->
            "VICAT_DB database=${database} metadata=${metadata.name}"
        }

        ch_vicat_database = PREPARE_VICAT_DATABASE.out.database
            .map { database, metadata -> tuple(database, metadata) }
            .first()

        // Serialize the two heavyweight builds; their memory requests add up.
        ch_vicat_nonviral_database_request = PREPARE_VICAT_DATABASE.out.database
            .take(1)
            .map { viralDatabase, viralSetupMetadata ->
                tuple(vicatNonviralDatabasePath, 'class-aware-nonviral',
                    nonviralSources[0], nonviralSources[1], nonviralSources[2], nonviralSources[3], nonviralHelpers)
            }

        PREPARE_VICAT_NONVIRAL_DATABASE(ch_vicat_nonviral_database_request)

        viewChannel(showChannelMessages, PREPARE_VICAT_NONVIRAL_DATABASE.out.database) { database, metadata ->
            "VICAT_NONVIRAL_DB database=${database} metadata=${metadata.name}"
        }

        ch_vicat_nonviral_database = PREPARE_VICAT_NONVIRAL_DATABASE.out.database
            .map { database, metadata -> tuple(database, metadata) }
            .first()

        PREDICT_VICAT_ORFS(NORMALIZE_FASTA.out.normalized_records)

        viewChannel(showChannelMessages, PREDICT_VICAT_ORFS.out.orfs) {
            prefix, type, proteins, orfMap, headerMap ->
                "VICAT_ORFS sample=${prefix} type=${type} proteins=${proteins.name}"
        }

        RUN_VICAT_DIAMOND(
            PREDICT_VICAT_ORFS.out.orfs,
            ch_vicat_database
        )

        RUN_VICAT_NONVIRAL_DIAMOND(
            PREDICT_VICAT_ORFS.out.orfs,
            ch_vicat_nonviral_database
        )

        viewChannel(showChannelMessages, RUN_VICAT_DIAMOND.out.results) {
            prefix, type, orfMap, headerMap, diamond ->
                "VICAT_DIAMOND sample=${prefix} type=${type} alignments=${diamond.name}"
        }

        ch_vicat_diamonds = RUN_VICAT_DIAMOND.out.results
            .map { prefix, type, orfMap, headerMap, diamond -> diamond }
            .collect()

        PREPARE_VICAT_REFERENCE_SUBSET(
            ch_vicat_diamonds,
            ch_vicat_database
        )

        viewChannel(showChannelMessages, PREPARE_VICAT_REFERENCE_SUBSET.out.subset) {
            referenceSubset, referenceSubsetMetadata ->
                "VICAT_REFERENCE_SUBSET references=${referenceSubset.name}"
        }

        // Both process inputs are value channels (`collect()` plus the
        // prepared database value), so this output is already reusable.
        ch_vicat_reference_subset = PREPARE_VICAT_REFERENCE_SUBSET.out.subset

        ch_vicat_dual_alignments = RUN_VICAT_DIAMOND.out.results
            .map { prefix, type, orfMap, headerMap, viralDiamond ->
                tuple(prefix, type, orfMap, headerMap, viralDiamond)
            }
            .join(
                RUN_VICAT_NONVIRAL_DIAMOND.out.results.map {
                    prefix, type, orfMap, headerMap, nonviralDiamond ->
                        tuple(prefix, nonviralDiamond)
                }
            )
            .map { prefix, type, orfMap, headerMap, viralDiamond, nonviralDiamond ->
                tuple(prefix, type, orfMap, headerMap, viralDiamond, nonviralDiamond)
            }

        PREPARE_VICAT_HITS(
            ch_vicat_dual_alignments,
            ch_vicat_reference_subset,
            ch_vicat_nonviral_database
        )

        viewChannel(showChannelMessages, PREPARE_VICAT_HITS.out.results) {
            prefix, type, orfMap, headerMap, preparedHits, preparedMetadata ->
                "VICAT_HITS_PREPARED sample=${prefix} hits=${preparedHits.name}"
        }

        // Treat the policy implementation as a task input. This makes a
        // decision-rule edit invalidate STANDARDIZE_VICAT under `-resume`
        // while preserving the expensive DIAMOND and hit-preparation cache.
        ch_vicat_standardizer_script = Channel.value(
            file("${projectDir}/bin/standardize_vicat.py")
        )
        STANDARDIZE_VICAT(
            PREPARE_VICAT_HITS.out.results,
            ch_vicat_standardizer_script
        )

        viewChannel(showChannelMessages, STANDARDIZE_VICAT.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_discovery_evidence = ch_discovery_evidence.mix(
            STANDARDIZE_VICAT.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_VICAT.out.evidence)
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_VICAT.out.context)
        // Both the RNA homology floor and audited CheckV exception consume loci.
        if( runVicat ) {
            ch_harmony_evidence = ch_harmony_evidence.mix(
                STANDARDIZE_VICAT.out.loci.combine(
                    NORMALIZE_FASTA.out.normalized_records
                        .filter { prefix, type, fasta, headerMap -> type == 'rna' }
                        .map { prefix, type, fasta, headerMap -> tuple(prefix, type) }, by: 0
                ).map { prefix, loci, type -> tuple(prefix, 'vicat_loci', loci) }
            )
        }
        ch_provirus_evidence = ch_provirus_evidence.mix(STANDARDIZE_VICAT.out.evidence)
        ch_provirus_evidence = ch_provirus_evidence.mix(
            STANDARDIZE_VICAT.out.provirus.map { prefix, tool, evidence ->
                tuple(prefix, 'vicat_boundary', evidence)
            }
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

        viewChannel(showChannelMessages, PREPARE_VITAP_DATABASE.out.database) { database, metadata ->
            "VITAP_DB database=${database} metadata=${metadata.name}"
        }

        // Database preparation emits once; convert that emission to a value
        // so every sample can reuse it rather than consuming a one-shot item.
        ch_vitap_database = PREPARE_VITAP_DATABASE.out.database
            .map { database, metadata -> tuple(database, metadata) }
            .first()
    }

    // Prepare or validate the persistent vConTACT3 reference database.
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
                vcontact3CleanupArchive
            )
        )

        PREPARE_VCONTACT3_DATABASE(ch_vcontact3_database_request)

        viewChannel(showChannelMessages, PREPARE_VCONTACT3_DATABASE.out.database) { database, metadata ->
            "VCONTACT3_DB database=${database} metadata=${metadata.name}"
        }

        ch_vcontact3_database = PREPARE_VCONTACT3_DATABASE.out.database
            .map { database, metadata -> tuple(database, metadata) }
            .first()
    }

    // A sized group key releases each sample as soon as all evidence expected
    // for that molecule type arrives. An unbounded groupTuple() would wait for
    // the entire mixed channel to close and let one slow sample delay all
    // otherwise-complete samples.
    ch_discovery_evidence_by_sample = ch_discovery_evidence
        .combine(
            NORMALIZE_FASTA.out.normalized_records.map {
                prefix, type, fasta, headerMap -> tuple(prefix, type)
            },
            by: 0
        )
        .map { prefix, tool, evidence, type ->
            tuple(
                groupKey(prefix, expectedDiscoveryToolsByType[type].size()),
                tool,
                evidence
            )
        }
        .groupTuple()
        .map { sampleKey, tools, evidenceFiles ->
            tuple(sampleKey.getGroupTarget(), tools, evidenceFiles)
        }

    ch_discovery_gate_inputs = NORMALIZE_FASTA.out.normalized_records
        .map { prefix, type, fasta, headerMap ->
            tuple(prefix, type, fasta, headerMap)
        }
        .join(ch_discovery_evidence_by_sample, remainder: true)
        .map { joined ->
            def expectedTools = expectedDiscoveryToolsByType[joined[1]]
            def actualTools = joined.size() >= 5 && joined[4] != null ? joined[4] : []
            def evidenceFiles = joined.size() >= 6 && joined[5] != null ? joined[5] : []
            validateEvidenceArtifacts(
                joined[0], 'discovery_gate', expectedTools, actualTools, evidenceFiles
            )

            if( evidenceFiles.isEmpty() ) {
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
                evidenceFiles.size(),
                evidenceFiles
            )
        }

    DISCOVERY_GATE(ch_discovery_gate_inputs)

    viewChannel(showChannelMessages, DISCOVERY_GATE.out.candidates) {
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

        viewChannel(showChannelMessages, PREPARE_CHECKV_DATABASE.out.database) { database, metadata ->
            "CHECKV_DB database=${database} metadata=${metadata.name}"
        }

        ch_checkv_database = PREPARE_CHECKV_DATABASE.out.database
            .map { database, metadata -> tuple(database, metadata) }
            .first()

        RUN_CHECKV(
            DISCOVERY_GATE.out.candidates,
            ch_checkv_database
        )

        viewChannel(showChannelMessages, RUN_CHECKV.out.results) {
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

        viewChannel(showChannelMessages, STANDARDIZE_CHECKV.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_provirus_evidence = ch_provirus_evidence.mix(
            STANDARDIZE_CHECKV.out.evidence
        )
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_CHECKV.out.evidence)
    }

    ch_provirus_evidence_by_sample = ch_provirus_evidence
        .map { prefix, tool, evidence ->
            tuple(groupKey(prefix, sharedProvirusTools.size()), tool, evidence)
        }
        .groupTuple()
        .map { sampleKey, tools, evidenceFiles ->
            tuple(sampleKey.getGroupTarget(), tools, evidenceFiles)
        }

    ch_provirus_refinement_inputs = DISCOVERY_GATE.out.candidates
        .join(ch_provirus_evidence_by_sample, remainder: true)
        .map { joined ->
            def actualTools = joined.size() >= 6 && joined[5] != null ? joined[5] : []
            def evidenceFiles = joined.size() >= 7 && joined[6] != null ? joined[6] : []
            validateEvidenceArtifacts(
                joined[0], 'provirus_refinement', sharedProvirusTools,
                actualTools, evidenceFiles
            )

            if( evidenceFiles.isEmpty() ) {
                return tuple(
                    joined[0],
                    joined[1],
                    joined[2],
                    joined[3],
                    joined[4],
                    0,
                    file("${projectDir}/assets/empty_discovery_evidence.tsv"),
                    allowCt3OnlyRefinement,
                    vicatProvirusMinOverlap
                )
            }
            tuple(
                joined[0],
                joined[1],
                joined[2],
                joined[3],
                joined[4],
                evidenceFiles.size(),
                evidenceFiles,
                allowCt3OnlyRefinement,
                vicatProvirusMinOverlap
            )
        }

    ch_provirus_refiner_script = Channel.value(
        file("${projectDir}/bin/refine_proviral_regions.py")
    )
    REFINE_PROVIRAL_REGIONS(
        ch_provirus_refinement_inputs,
        ch_provirus_refiner_script
    )

    viewChannel(showChannelMessages, REFINE_PROVIRAL_REGIONS.out.refined) {
        prefix, type, refinedFasta, regionMap, boundaryAudit, summary ->
            "REFINED sample=${prefix} type=${type} fasta=${refinedFasta.name} map=${regionMap.name}"
    }

    if( runVicat ) {
        ch_vicat_refinement_projection = REFINE_PROVIRAL_REGIONS.out.refined
            .map { prefix, type, refinedFasta, regionMap, boundaryAudit, summary ->
                tuple(prefix, type, regionMap)
            }
            .join(
                STANDARDIZE_VICAT.out.loci.map { prefix, loci -> tuple(prefix, loci) }
            )
            .map { prefix, type, regionMap, loci ->
                tuple(prefix, type, regionMap, loci)
            }

        PROJECT_VICAT_REFINED(ch_vicat_refinement_projection)

        viewChannel(showChannelMessages, PROJECT_VICAT_REFINED.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} scope=refined_region evidence=${evidence.name}"
        }
        ch_harmony_evidence = ch_harmony_evidence.mix(
            PROJECT_VICAT_REFINED.out.evidence
        )
    }

    if( runTesorter ) {
        RUN_TESORTER(REFINE_PROVIRAL_REGIONS.out.refined)

        viewChannel(showChannelMessages, RUN_TESORTER.out.results) {
            prefix, type, regionMap, sequenceMap, classifications, domains, domainGff, log, metadata ->
                "TESORTER sample=${prefix} type=${type} classifications=${classifications.name}"
        }

        STANDARDIZE_TESORTER(RUN_TESORTER.out.results)

        viewChannel(showChannelMessages, STANDARDIZE_TESORTER.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }

        ch_harmony_evidence = ch_harmony_evidence.mix(
            STANDARDIZE_TESORTER.out.evidence
        )

        // Keep the refined candidate tuple unchanged, but release it only
        // after TEsorter completes for the same sample.
        ch_refined_for_taxonomy = REFINE_PROVIRAL_REGIONS.out.refined
            .join(RUN_TESORTER.out.completed)
            .map { prefix, type, refinedFasta, regionMap, boundaryAudit,
                   refinementSummary, completedType ->
                if( type != completedType ) {
                    error "TEsorter completion type mismatch for sample '${prefix}'"
                }
                tuple(
                    prefix, type, refinedFasta, regionMap,
                    boundaryAudit, refinementSummary
                )
            }
    } else {
        ch_refined_for_taxonomy = REFINE_PROVIRAL_REGIONS.out.refined
    }

    if( runVicat ) {
        // This is a completion barrier, not another biological vote. It
        // guarantees that every viCAT-enabled sample has completed the
        // coordinate projection before viHARMONY is allowed to run.
        ch_refined_for_harmony = ch_refined_for_taxonomy
            .join(
                PROJECT_VICAT_REFINED.out.projection.map {
                    prefix, projection -> tuple(prefix, projection)
                }
            )
            .map { prefix, type, refinedFasta, regionMap, boundaryAudit,
                   refinementSummary, projection ->
                tuple(
                    prefix, type, refinedFasta, regionMap,
                    boundaryAudit, refinementSummary
                )
            }
    } else {
        ch_refined_for_harmony = ch_refined_for_taxonomy
    }

    if( runVitap ) {
        RUN_VITAP(
            ch_refined_for_taxonomy,
            ch_vitap_database,
            vitapIncludeLowConfidence
        )

        viewChannel(showChannelMessages, RUN_VITAP.out.results) {
            prefix, type, regionMap, best, allLineages, fallback, log, metadata ->
                "VITAP sample=${prefix} type=${type} best=${best.name} all=${allLineages.name}"
        }

        STANDARDIZE_VITAP(
            RUN_VITAP.out.results,
            ch_vitap_database
        )

        viewChannel(showChannelMessages, STANDARDIZE_VITAP.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_VITAP.out.evidence)
    }

    if( runVcontact3 ) {
        RUN_VCONTACT3(
            ch_refined_for_taxonomy,
            ch_vcontact3_database,
            vcontact3DbDomain
        )

        viewChannel(showChannelMessages, RUN_VCONTACT3.out.results) {
            prefix, type, regionMap, assignments, metrics, logs, metadata ->
                "VCONTACT3 sample=${prefix} type=${type} domain_mode=${vcontact3DbDomain} assignments=${assignments}"
        }

        STANDARDIZE_VCONTACT3(RUN_VCONTACT3.out.results)

        viewChannel(showChannelMessages, STANDARDIZE_VCONTACT3.out.evidence) { prefix, tool, evidence ->
            "STANDARDIZED sample=${prefix} tool=${tool} evidence=${evidence.name}"
        }
        ch_harmony_evidence = ch_harmony_evidence.mix(STANDARDIZE_VCONTACT3.out.evidence)
    }

    // viCAT contributes both parent-contig and post-refinement evidence, so
    // its second artifact is included in the expected emission count even
    // though both artifacts deliberately retain the canonical `vicat` label.
    def extraHarmonyEvidenceCount = runVicat ? 1 : 0
    ch_harmony_evidence_by_sample = ch_harmony_evidence
        .combine(
            NORMALIZE_FASTA.out.normalized_records.map {
                prefix, type, fasta, headerMap -> tuple(prefix, type)
            },
            by: 0
        )
        .map { prefix, tool, evidence, type ->
            tuple(
                groupKey(
                    prefix,
                    expectedHarmonyToolsByType[type].size() +
                        extraHarmonyEvidenceCount
                ),
                tool,
                evidence
            )
        }
        .groupTuple()
        .map { sampleKey, tools, evidenceFiles ->
            tuple(sampleKey.getGroupTarget(), tools, evidenceFiles)
        }

    ch_harmony_evidence_for_sample = NORMALIZE_FASTA.out.normalized_records
        .map { prefix, type, fasta, headerMap ->
            tuple(prefix, expectedHarmonyToolsByType[type])
        }
        .join(ch_harmony_evidence_by_sample, remainder: true)
        .map { joined ->
            def expectedTools = joined[1]
            def actualTools = joined.size() >= 3 && joined[2] != null ? joined[2] : []
            def evidenceFiles = joined.size() >= 4 && joined[3] != null ? joined[3] : []
            validateEvidenceArtifacts(
                joined[0], 'viharmony', expectedTools, actualTools, evidenceFiles
            )

            if( evidenceFiles.isEmpty() ) {
                return tuple(
                    joined[0],
                    0,
                    file("${projectDir}/assets/empty_discovery_evidence.tsv")
                )
            }
            tuple(joined[0], evidenceFiles.size(), evidenceFiles)
        }

    if( runVcontact3 ) {
        ch_vcontact3_groups_for_harmony = STANDARDIZE_VCONTACT3.out.groups
            .map { prefix, groups -> tuple(prefix, 1, groups) }
    } else {
        ch_vcontact3_groups_for_harmony = ch_samples.map { prefix, type, fasta ->
            tuple(prefix, 0, file("${projectDir}/assets/empty_vcontact3_groups.tsv"))
        }
    }

    ch_harmony_inputs = ch_refined_for_harmony
        .map { prefix, type, refined, regionMap, boundaryAudit, summary ->
            tuple(prefix, type, refined, regionMap)
        }
        .join(
            NORMALIZE_FASTA.out.normalized_records.map {
                prefix, type, fasta, headerMap -> tuple(prefix, fasta, headerMap)
            }
        )
        .join(
            DISCOVERY_GATE.out.candidates.map {
                prefix, type, candidates, audit, summary -> tuple(prefix, audit)
            }
        )
        .join(ch_harmony_evidence_for_sample)
        .join(ch_vcontact3_groups_for_harmony)
        .map { prefix, type, refined, regionMap, normalized, headerMap, discoveryAudit,
               evidenceFileCount, evidenceFiles, groupCount, groupFiles ->
            tuple(
                prefix, type, normalized, headerMap, discoveryAudit, refined,
                regionMap, evidenceFileCount, evidenceFiles, groupCount,
                groupFiles, ictvCsv, harmonizerAudit,
                vcontact3MinTaxonomyLength
            )
        }

    VIHARMONY(ch_harmony_inputs)

    viewChannel(showChannelMessages, VIHARMONY.out.results) { prefix, normalizedFasta, originalFasta, metadata,
                                reviewQueue, disposition, sequenceMap, manifest,
                                databaseFasta, databaseMetadata, allCandidatesFasta,
                                reviewCandidatesFasta, provisionalNormalizedFasta,
                                provisionalOriginalFasta, provisionalMetadata ->
        "VIHARMONY sample=${prefix} fasta=${normalizedFasta.name} metadata=${metadata.name}"
    }
}
