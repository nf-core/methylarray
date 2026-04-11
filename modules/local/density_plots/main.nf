process DENSITY_PLOTS {

    tag "density_plots"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/density_plots/environment.yml"

    input:
    path rg_rds
    path beta_norm_rds

    output:
    path "06_density_beta_raw.png"   , emit: density_raw_plot
    path "07_density_beta_sesame.png", emit: density_sesame_plot
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    density_plots.R \\
        --rg        ${rg_rds} \\
        --beta-norm ${beta_norm_rds} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch 06_density_beta_raw.png
    touch 07_density_beta_sesame.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-minfi: 1.48.0
    END_VERSIONS
    """
}
