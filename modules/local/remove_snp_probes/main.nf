process REMOVE_SNP_PROBES {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ 'docker.io/ajandriaa/methylarray:0.0.2dev' }"

    input:
    tuple val(samplesheet_name), path(RData_XREACTIVE)

    output:
    tuple val(samplesheet_name), path("mVals_noXprob_noSNP.csv")  , emit: csv_mVals
    tuple val(samplesheet_name), path("bVals_noXprob_noSNP.csv")  , emit: csv_bVals
    tuple val(samplesheet_name), path("mSetSqFlt_noXprob_noSNP.RData"), emit: rdata

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'remove_snp_probes_3.R'
}
