process ALIGN_META {

    tag "align_meta"
    label 'process_low'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/align_meta/environment.yml"

    input:
    path bvals_rds
    path mvals_rds
    path targets_tsv
    path meta_csv
    path cell_counts_csv  // may be [] (empty list) when do_estimate_cellcomp = false
    path pred_sex_csv

    output:
    path "meta_aligned.csv"              , emit: meta_aligned_csv
    path "bVals_aligned.rds"             , emit: bvals_aligned_rds
    path "mVals_aligned.rds"             , emit: mvals_aligned_rds
    path "bVals_unmapped_colnames.csv"   , optional: true, emit: bvals_unmapped_csv
    path "bVals_aligned.csv"             , optional: true, emit: bvals_aligned_csv
    path "mVals_aligned.csv"             , optional: true, emit: mvals_aligned_csv
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def cell_arg = cell_counts_csv instanceof List ? '' : "--cell-counts ${cell_counts_csv}"
    """
    align_meta.R \\
        --bvals    ${bvals_rds} \\
        --mvals    ${mvals_rds} \\
        --targets  ${targets_tsv} \\
        --meta     ${meta_csv} \\
        ${cell_arg} \\
        --pred-sex ${pred_sex_csv} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """

    stub:
    """
    touch meta_aligned.csv
    touch bVals_aligned.rds
    touch mVals_aligned.rds
    touch bVals_unmapped_colnames.csv
    touch bVals_aligned.csv
    touch mVals_aligned.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
    END_VERSIONS
    """
}
