process REMOVE_SEX_CHROMOSOMES {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ params.methylarray_deps_container }"

    input:
    tuple val(samplesheet_name), path(RData_SNPPROBES)
    tuple val(samplesheet_name), path(RData_rgSet)

    output:
    tuple val(samplesheet_name), path("mVals.RData")    , emit: mVals
    tuple val(samplesheet_name), path("bVals.RData")    , emit: bVals
    tuple val(samplesheet_name), path("mVals.csv")      , emit: mVals_csv
    tuple val(samplesheet_name), path("bVals.csv")      , emit: bVals_csv
    tuple val(samplesheet_name), path("mSetSqFlt.RData"), emit: mSetSqFlt
    tuple val(samplesheet_name), path("rgSet.RData")    , emit: rgSet

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'remove_sex_chromosomes_o1.R'
}
