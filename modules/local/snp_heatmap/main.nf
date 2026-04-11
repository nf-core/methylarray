process SNP_HEATMAP {

    tag "snp_heatmap"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/snp_heatmap/environment.yml"

    input:
    path rg_rds

    output:
    path "05_heatmap_snp_correlation.png", optional: true, emit: snp_heatmap_plot
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    snp_heatmap.R --rg ${rg_rds} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch 05_heatmap_snp_correlation.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-minfi: 1.48.0
    END_VERSIONS
    """
}
