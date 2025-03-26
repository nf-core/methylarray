process XREACTIVE_PROBES_FIND_REMOVE {
    tag "${RData_PREPROCESSING.baseName}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ params.methylarray_deps_container }"

    input:
    tuple val(samplesheet_name), path(RData_PREPROCESSING)
    path(genome_path)

    output:
    tuple val(samplesheet_name), path("*.csv")  , emit: csv
    tuple val(samplesheet_name), path("mSetSqFlt.RData"), emit: rdata

    when:
    task.ext.when == null || task.ext.when

    script:
    conda = workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    if (conda) {
        log.info "\u001B[32mINFO: XREACTIVE_PROBES_FIND_REMOVE will attempt to fetch DNAmCrosshyb R package from forked GitHub release as it is not hosted on Conda.\u001B[0m"
    }
    chrom_number = params.xreactive_chr_targets ? params.xreactive_chr_targets : 'all'
    template "xreactive_probes_find_remove_2.R"
}
