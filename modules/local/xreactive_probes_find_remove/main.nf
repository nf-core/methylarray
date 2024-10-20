process XREACTIVE_PROBES_FIND_REMOVE {
    tag "${RData_PREPROCESSING.baseName}"
    label 'process_medium'

    container "${ params.methylarray_deps_container }"

    input:
    tuple val(samplesheet_name), path(RData_PREPROCESSING)
    path(genome_path)

    output:
    tuple val(samplesheet_name), path("*.csv")  , emit: csv
    tuple val(samplesheet_name), path("mSetSqFlt_noXprob.RData"), emit: rdata

    when:
    task.ext.when == null || task.ext.when
    
    script:
    chrom_number = params.xreactive_chr_targets ? params.xreactive_chr_targets : 'all'
    template "xreactive_probes_find_remove_2.R"
}
