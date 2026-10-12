process FINAL_EXPORTS {

    tag "final_exports"
    label 'process_low'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/final_exports/environment.yml"

    input:
    path beta_clean_rds

    output:
    path "bVals_final.rds"  , emit: bvals_rds
    path "mVals_final.rds"  , emit: mvals_rds
    path "samples_final.csv", emit: samples_csv
    path "probes_final.csv" , emit: probes_csv
    path "qc_summary.txt"   , optional: true, emit: qc_summary_txt
    path "bVals_final.csv"  , optional: true, emit: bvals_csv
    path "mVals_final.csv"  , optional: true, emit: mvals_csv
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    final_exports.R \\
        --beta-clean ${beta_clean_rds} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """

    stub:
    """
    touch bVals_final.rds
    touch mVals_final.rds
    touch samples_final.csv
    touch probes_final.csv
    touch qc_summary.txt
    touch bVals_final.csv
    touch mVals_final.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
    END_VERSIONS
    """
}
