process CELL_COMP {

    tag "cell_comp"
    label 'process_high'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/cell_comp/environment.yml"

    input:
    path rg_rds

    output:
    path "cell_counts_blood.csv", emit: cell_counts_csv
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cell_comp.R \\
        --rg  ${rg_rds} \\
        --out cell_counts_blood.csv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
        bioconductor-flowsorted.blood.epic: \$(Rscript -e "cat(as.character(packageVersion('FlowSorted.Blood.EPIC')))")
    END_VERSIONS
    """

    stub:
    """
    touch cell_counts_blood.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-minfi: 1.48.0
        bioconductor-flowsorted.blood.epic: 2.4.0
    END_VERSIONS
    """
}
