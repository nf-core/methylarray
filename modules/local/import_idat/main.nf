process IMPORT_IDAT {

    tag { sample_sheet.baseName }
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/import_idat/environment.yml"

    input:
    path sample_sheet
    path idat_dir

    output:
    path "rg0.rds"      , emit: rg_rds
    path "targets.tsv"  , emit: targets_tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    import_idat.R \\
        --sample-sheet ${sample_sheet} \\
        --idat-dir     ${idat_dir} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch rg0.rds
    touch targets.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-minfi: 1.48.0
    END_VERSIONS
    """
}
