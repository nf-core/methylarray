process QC_REPORT_CONTROLS {

    tag "qc_report_controls"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/qc_report_controls/environment.yml"

    input:
    path rg_rds

    output:
    path "00_qcReport.pdf"            , optional: true, emit: qc_report_pdf
    path "00_controlStrip_*.png"      , optional: true, emit: control_strip_plots
    path "03_control_bisulfite_I.png" , optional: true, emit: bisulfite_i_plot
    path "04_control_bisulfite_II.png", optional: true, emit: bisulfite_ii_plot
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    qc_report_controls.R --rg ${rg_rds} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch 00_qcReport.pdf
    touch "00_controlStrip_stub.png"
    touch 03_control_bisulfite_I.png
    touch 04_control_bisulfite_II.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-minfi: 1.48.0
    END_VERSIONS
    """
}
