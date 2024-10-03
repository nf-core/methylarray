process XREACTIVE_PROBES_FIND_REMOVE {
    tag "${RData_PREPROCESSING.baseName}"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    conda "/Users/adrian/anaconda3/envs/methylarray" // For local development only until this is resolved: https://github.com/Bioconductor/bioconductor_docker/issues/22
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'oras://community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicmanifest_bioconductor-minfi_pruned:8e9311656c17654c' :
    //    'community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicmanifest_bioconductor-minfi_pruned:2cf4967a5ce4c623' }"

    input:
    tuple val(samplesheet_name), val(RData_PREPROCESSING)

    output:
    tuple val(samplesheet_name), path("*.csv")  , emit: csv
    tuple val(samplesheet_name), path("mSetSqFlt_noXprob.RData"), emit: rdata

    when:
    task.ext.when == null || task.ext.when
    
    script:
    // Download from: https://github.com/pjhop/DNAmCrosshyb at https://doi.org/10.5281/zenodo.4088019
    def genome_path = "${projectDir}/results/genome_bs/hg19"
    """
    Rscript ${moduleDir}/templates/xreactive_probes_find_remove_2.R ${RData_PREPROCESSING} ${genome_path}
    """
}
