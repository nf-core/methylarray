process REMOVE_SNP_PROBES {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ 'docker.io/ajandriaa/methylarray:0.0.2dev' }"

    input:
    tuple val(samplesheet_name), path(RData_XREACTIVE)

    output:
    tuple val(samplesheet_name), path("*.csv")  , emit: csv
    tuple val(samplesheet_name), path("*.RData"), emit: rdata

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'remove_snp_probes_3.R'
}
