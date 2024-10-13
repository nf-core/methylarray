process FIND_DMR {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ 'docker.io/ajandriaa/methylarray:0.0.4dev' }"

    input:
    tuple val(samplesheet_name), path(bVALS_SNPPROBES)
    path(extensive_metadata)

    output:
    tuple val(samplesheet_name), path("dmp_champ.%s.csv"), emit: mVals
    tuple val(samplesheet_name), path("dmp_minfi.csv"), emit: bVals

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'find_dmr_5.R'
}
