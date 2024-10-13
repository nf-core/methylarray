process ADJUST_BATCH_EFFECT {
    tag "${samplesheet_name}"
    label 'process_single'

    container "${ 'docker.io/ajandriaa/methylarray:0.0.4dev' }"

    input:
    tuple val(samplesheet_name), path(bVALS_SNPPROBES)
    path(extensive_metadata)

    output:
    tuple val(samplesheet_name), path("mVals.cell_comp.cor_bmi_age.csv"), emit: mVals
    tuple val(samplesheet_name), path("bVals.cell_comp.cor_bmi_age.csv"), emit: bVals

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'adjust_batch_effect_o4.R'
}
