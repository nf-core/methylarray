process REMOVE_SNP_PROBES {
    tag "${samplesheet_name}"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    conda "/Users/adrian/anaconda3/envs/methylarray" // For local development only until this is resolved: https://github.com/Bioconductor/bioconductor_docker/issues/22
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'oras://community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicmanifest_bioconductor-minfi_pruned:8e9311656c17654c' :
    //    'community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicmanifest_bioconductor-minfi_pruned:2cf4967a5ce4c623' }"

    input:
    tuple val(samplesheet_name), val(RData_XREACTIVE)

    output:
    tuple val(samplesheet_name), path("*.csv")  , emit: csv
    tuple val(samplesheet_name), path("*.RData"), emit: rdata

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${moduleDir}/templates/remove_snp_probes_3.R ${RData_XREACTIVE}
    """
}
