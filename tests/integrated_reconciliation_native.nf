// Packaged at project root as native_validation.nf. No search modules imported.
nextflow.enable.dsl = 2
include { VIHARMONY } from './modules/local/viharmony'

process STAGE_SAVED_INPUTS {
    tag "$sample"
    cpus 1
    memory '256 MB'
    input:
    tuple val(sample), path(saved_files)
    output:
    tuple val(sample), path('staged/*')
    script:
    def sources = saved_files.collect { "'${it}'" }.join(' ')
    """
    mkdir staged
    cp -- ${sources} staged/
    """
}

workflow {
    ch_saved = Channel.fromList(params.saved_specs).map { spec ->
        tuple(spec.sample, spec.files.collect { file(it, checkIfExists: true) })
    }
    STAGE_SAVED_INPUTS(ch_saved)
    ch_harmony = STAGE_SAVED_INPUTS.out.map { sample, staged ->
        def spec = params.saved_specs.find { it.sample == sample }
        def files = staged.collectEntries { [(it.name): it] }
        tuple(sample, spec.type, files[spec.normalized], files[spec.header],
            files[spec.gate], files[spec.refined], files[spec.regions],
            spec.evidence.size(), spec.evidence.collect { files[it] },
            1, [files[spec.groups]], file(params.ictv_csv, checkIfExists: true),
            file(params.ictv_reconciliation, checkIfExists: true),
            file(params.ictv_reconciliation + '.manifest.json', checkIfExists: true),
            'full', spec.minimum_length)
    }
    VIHARMONY(ch_harmony)
}
