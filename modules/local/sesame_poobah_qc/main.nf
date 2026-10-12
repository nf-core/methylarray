process SESAME_POOBAH_QC {

    tag "sesame_poobah_qc"
    label 'process_high'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/sesame_poobah_qc/environment.yml"

    input:
    path rg0_rds
    path targets_tsv
    path idat_dir

    output:
    path "beta_norm.rds"                           , emit: beta_norm_rds
    path "rg_filtered.rds"                         , emit: rg_filt_rds
    path "targets_filtered.tsv"                    , emit: targets_filt_tsv
    path "02_hist_failrate_per_sample_masking.png" , emit: failrate_plot
    path "masking_failrate_per_sample.csv"         , emit: failrate_per_sample_csv
    path "masking_failrate_per_probe.csv"          , emit: failrate_per_probe_csv
    path "dropped_samples.csv"                     , optional: true, emit: dropped_samples_csv
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sesame_poobah_qc.R \\
        --rg0         ${rg0_rds} \\
        --targets     ${targets_tsv} \\
        --out-beta    beta_norm.rds \\
        --out-rg      rg_filtered.rds \\
        --out-targets targets_filtered.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-sesame: \$(Rscript -e "cat(as.character(packageVersion('sesame')))")
        r-minfi: \$(Rscript -e "cat(as.character(packageVersion('minfi')))")
    END_VERSIONS
    """

    stub:
    """
    touch beta_norm.rds
    touch rg_filtered.rds
    touch targets_filtered.tsv
    touch 02_hist_failrate_per_sample_masking.png
    touch masking_failrate_per_sample.csv
    touch masking_failrate_per_probe.csv
    touch dropped_samples.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-sesame: 1.18.0
        r-minfi: 1.48.0
    END_VERSIONS
    """
}
