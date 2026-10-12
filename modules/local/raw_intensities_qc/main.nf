process RAW_INTENSITIES_QC {

    tag "raw_intensities"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/raw_intensities_qc/environment.yml"

    input:
    path rg_rds

    output:
    path "01_plotQC_raw_intensities.png" , emit: qc_plot
    path "qc_raw.tsv"                    , emit: qc_tsv
    path "low_intensity_samples.tsv"     , optional: true, emit: low_intensity_tsv
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    raw_intensities_qc.R --rgset ${rg_rds} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch 01_plotQC_raw_intensities.png
    touch qc_raw.tsv
    touch low_intensity_samples.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-minfi: 1.48.0
    END_VERSIONS
    """
}
