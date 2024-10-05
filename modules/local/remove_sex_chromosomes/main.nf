process REMOVE_SEX_CHROMOSOMES {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ 'docker.io/ajandriaa/methylarray:0.0.3dev' }"

    input:
    tuple val(samplesheet_name), path(RData_SNPPROBES)
    tuple val(samplesheet_name), path(RData_rgSet)

    output:
    tuple val(samplesheet_name), path("*.csv")  , emit: csv
    tuple val(samplesheet_name), path("mSetSqFlt.RData"), emit: rdata_mSetSqFlt
    tuple val(samplesheet_name), path("rgSet.RData"), emit: rdata_rgSet

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'remove_sex_chromosomes_o1.R'
}
