process COVARIATE_PCA {

    tag "covariate_pca"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/covariate_pca/environment.yml"

    input:
    path bvals_aligned_rds
    path meta_aligned_csv

    output:
    path "covariates_PC_association_minuslog10p.csv", emit: pc_association_csv
    path "covariate_checks.txt"                     , emit: covariate_checks_txt
    path "ChAMP_SVD"                                , optional: true, emit: champ_svd_dir
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    covariate_pca.R \\
        --bvals ${bvals_aligned_rds} \\
        --meta  ${meta_aligned_csv} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        bioconductor-champ: \$(Rscript -e "cat(as.character(packageVersion('ChAMP')))")
        bioconductor-sva: \$(Rscript -e "cat(as.character(packageVersion('sva')))")
    END_VERSIONS
    """

    stub:
    """
    touch covariates_PC_association_minuslog10p.csv
    touch covariate_checks.txt
    mkdir ChAMP_SVD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        bioconductor-champ: 2.32.0
        bioconductor-sva: 3.50.0
    END_VERSIONS
    """
}
