process REMOVE_CONFOUNDING_PROBES {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ params.methylarray_deps_container }"

    input:
    tuple val(samplesheet_name), path(mVals)
    tuple val(samplesheet_name), path(bVals)
    tuple val(samplesheet_name), path(RData_REMOVESNP)
    path(extensive_metadata)

    output:
    tuple val(samplesheet_name), path("cbVals.filtered_probes.csv")     , emit: bVals
    tuple val(samplesheet_name), path("cmVals.filtered_probes.csv")     , emit: mVals
    tuple val(samplesheet_name), path("mSetSqFlt.filtered_probes.RData"), emit: rdata

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'remove_confounding_probes_o2.R'
}
