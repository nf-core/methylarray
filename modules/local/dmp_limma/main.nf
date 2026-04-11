process DMP_LIMMA {

    tag "dmp_limma"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/dmp_limma/environment.yml"

    input:
    path mvals_aligned_rds
    path bvals_aligned_rds
    path meta_aligned_csv

    output:
    path "DMP_limma__plate_in_model", optional: true, emit: dmp_plate_dir
    path "DMP_limma__combat"        , optional: true, emit: dmp_combat_dir
    path "DMP_limma__norm_only"     , emit: dmp_norm_dir
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    dmp_limma.R \\
        --mvals ${mvals_aligned_rds} \\
        --bvals ${bvals_aligned_rds} \\
        --meta  ${meta_aligned_csv} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-limma: \$(Rscript -e "cat(as.character(packageVersion('limma')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir DMP_limma__plate_in_model
    mkdir DMP_limma__combat
    mkdir DMP_limma__norm_only

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-limma: 3.58.0
    END_VERSIONS
    """
}
