process FILTER_XY_NONCG_BEADS {

    tag "filter_xy_noncg_beads"
    label 'process_high'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/filter_xy_noncg_beads/environment.yml"

    input:
    path rg_filtered_rds
    path beta_norm_rds

    output:
    path "beta3.rds"             , emit: beta3_rds
    path "filtering_summary.tsv" , emit: filtering_summary_tsv
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    filter_xy_noncg_beads.R \\
        --rg          ${rg_filtered_rds} \\
        --beta-norm   ${beta_norm_rds} \\
        --out-beta    beta3.rds \\
        --out-summary filtering_summary.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch beta3.rds
    touch filtering_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-minfi: 1.48.0
    END_VERSIONS
    """
}
