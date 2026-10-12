process SPLIT_COLLAPSE {

    tag "split_collapse"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/split_collapse/environment.yml"

    input:
    path beta3_rds

    output:
    path "beta_for_DMR_uncollapsed.rds", emit: beta_dmr_rds
    path "beta_clean.rds"              , emit: beta_clean_rds
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    split_collapse.R \\
        --beta3     ${beta3_rds} \\
        --out-dmr   beta_for_DMR_uncollapsed.rds \\
        --out-clean beta_clean.rds \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-dmrcate: \$(Rscript -e "cat(as.character(packageVersion('DMRcate')))")
        r-sesame: \$(Rscript -e "cat(as.character(packageVersion('sesame')))")
    END_VERSIONS
    """

    stub:
    """
    touch beta_for_DMR_uncollapsed.rds
    touch beta_clean.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-dmrcate: 2.16.0
        r-sesame: 1.18.0
    END_VERSIONS
    """
}
