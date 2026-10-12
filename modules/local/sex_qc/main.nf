process SEX_QC {

    tag "sex_qc"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/sex_qc/environment.yml"

    input:
    path rg_rds
    path targets_tsv
    path meta_file

    output:
    path "qc_predicted_sex.csv"          , emit: qc_pred_sex_csv
    path "08_plotSex_predicted.png"      , optional: true, emit: sex_plot
    path "sex_discordant_samples.csv"    , optional: true, emit: sex_discordant_csv
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sex_qc.R \\
        --rg      ${rg_rds} \\
        --targets ${targets_tsv} \\
        --meta    '${meta_file}' \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch qc_predicted_sex.csv
    touch 08_plotSex_predicted.png
    touch sex_discordant_samples.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-minfi: 1.48.0
    END_VERSIONS
    """
}
