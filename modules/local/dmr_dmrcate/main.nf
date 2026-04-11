process DMR_DMRCATE {

    tag "dmr_dmrcate"
    label 'process_medium'
    container 'quay.io/nf-core/methylarray:1.0.0dev'
    conda "${projectDir}/modules/local/dmr_dmrcate/environment.yml"

    input:
    path beta_dmr_uncollapsed_rds
    path targets_tsv
    path meta_aligned_csv

    output:
    path "DMR_DMRcate_EPICv2__plate_in_model", optional: true, emit: dmr_plate_dir
    path "DMR_DMRcate_EPICv2__combat"        , optional: true, emit: dmr_combat_dir
    path "DMR_DMRcate_EPICv2__norm_only"     , emit: dmr_norm_dir
    path "DMR_coMethDMR__plate_in_model"     , optional: true, emit: cometh_plate_dir
    path "DMR_coMethDMR__combat"             , optional: true, emit: cometh_combat_dir
    path "DMR_coMethDMR__norm_only"          , optional: true, emit: cometh_norm_dir
    path "DMR_consensus__plate_in_model"     , optional: true, emit: consensus_plate_dir
    path "DMR_consensus__combat"             , optional: true, emit: consensus_combat_dir
    path "DMR_consensus__norm_only"          , optional: true, emit: consensus_norm_dir
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def run_cometh = task.ext.run_cometh ?: 'false'
    """
    dmr_dmrcate.R \\
        --beta-uc    ${beta_dmr_uncollapsed_rds} \\
        --targets    ${targets_tsv} \\
        --meta       ${meta_aligned_csv} \\
        --run-cometh ${run_cometh} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | sed 's/R version //' | sed 's/ (.*//')
        r-dmrcate: \$(Rscript -e "cat(as.character(packageVersion('DMRcate')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir DMR_DMRcate_EPICv2__plate_in_model
    mkdir DMR_DMRcate_EPICv2__combat
    mkdir DMR_DMRcate_EPICv2__norm_only
    mkdir DMR_coMethDMR__plate_in_model
    mkdir DMR_coMethDMR__combat
    mkdir DMR_coMethDMR__norm_only
    mkdir DMR_consensus__plate_in_model
    mkdir DMR_consensus__combat
    mkdir DMR_consensus__norm_only

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        r-dmrcate: 2.16.0
    END_VERSIONS
    """
}
