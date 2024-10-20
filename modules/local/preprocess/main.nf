process PREPROCESS {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ params.methylarray_deps_container }"

    input:
    tuple path(idat_folders), val(samplesheet_name)

    output:
    tuple val(samplesheet_name), path("*.csv")  , emit: csv
    tuple val(samplesheet_name), path("mSetSqFlt.RData"), emit: rdata
    tuple val(samplesheet_name), path("rgSet.RData"), emit: rdata_rgSet

    when:
    task.ext.when == null || task.ext.when

    script:
    template "preprocess_1.R"
}
