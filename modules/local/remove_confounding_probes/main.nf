process REMOVE_CONFOUNDING_PROBES {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ 'docker.io/ajandriaa/methylarray:0.0.4dev' }" // TODO: Update container for ChAMP

    input:
    tuple val(samplesheet_name), path(RData_mVals)
    tuple val(samplesheet_name), path(RData_bVals)
    tuple val(samplesheet_name), path(RData_REMOVESNP)
    path(extensive_metadata)

    output:
    tuple val(samplesheet_name), path("*")  , emit: out

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'remove_confounding_probes_o2.R'
}
