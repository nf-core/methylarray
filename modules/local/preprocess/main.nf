process PREPROCESS {
    tag "${samplesheet_name}"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    conda "/Users/adrian/anaconda3/envs/methylarray" // For local development only until this is resolved: https://github.com/Bioconductor/bioconductor_docker/issues/22
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'oras://community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicmanifest_bioconductor-minfi_pruned:8e9311656c17654c' :
    //    'community.wave.seqera.io/library/bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19_bioconductor-illuminahumanmethylation450kmanifest_bioconductor-illuminahumanmethylationepicmanifest_bioconductor-minfi_pruned:2cf4967a5ce4c623' }"

    input:
    tuple val(idat_folders), val(samplesheet_name)

    output:
    tuple val(samplesheet_name), path("*.csv")  , emit: csv
    tuple val(samplesheet_name), path("mSetSqFlt.RData"), emit: rdata

    when:
    task.ext.when == null || task.ext.when
    
    // For -profile test
    idat_folders = params.test_data ? "${projectDir}/${idat_folders}" : idat_folders

    script:
    """
    Rscript ${moduleDir}/templates/preprocess_1.R ${idat_folders} ${samplesheet_name}
    """
}
