process LOG_VERSIONS {

    tag "log_versions"
    label 'process_single'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/log_versions/environment.yml"

    output:
    path "package_versions.csv"   , emit: package_versions_csv
    path "sessionInfo.txt"        , emit: session_info_txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    log_versions.R \\
        --out-versions package_versions.csv \\
        --out-session  sessionInfo.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """

    stub:
    """
    touch package_versions.csv
    touch sessionInfo.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
    END_VERSIONS
    """
}
