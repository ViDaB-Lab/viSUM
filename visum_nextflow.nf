#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NORMALIZE_FASTA } from './modules/local/normalize_fasta'


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
}
