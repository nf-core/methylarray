/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { LOG_VERSIONS           } from '../modules/local/log_versions'
include { IMPORT_IDAT            } from '../modules/local/import_idat'
include { RAW_INTENSITIES_QC     } from '../modules/local/raw_intensities_qc'
include { SESAME_POOBAH_QC       } from '../modules/local/sesame_poobah_qc'
include { QC_REPORT_CONTROLS     } from '../modules/local/qc_report_controls'
include { SNP_HEATMAP            } from '../modules/local/snp_heatmap'
include { DENSITY_PLOTS          } from '../modules/local/density_plots'
include { SEX_QC                 } from '../modules/local/sex_qc'
include { FILTER_XY_NONCG_BEADS  } from '../modules/local/filter_xy_noncg_beads'
include { CELL_COMP              } from '../modules/local/cell_comp'
include { SPLIT_COLLAPSE         } from '../modules/local/split_collapse'
include { FINAL_EXPORTS          } from '../modules/local/final_exports'
include { ALIGN_META             } from '../modules/local/align_meta'
include { COVARIATE_PCA          } from '../modules/local/covariate_pca'
include { DMP_LIMMA              } from '../modules/local/dmp_limma'
include { DMR_DMRCATE            } from '../modules/local/dmr_dmrcate'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLARRAY {

    take:
    ch_samplesheet // channel: path to sample sheet (--input / params.input)
    ch_idatdir     // channel: path to IDAT directory (--idat_dir)
    ch_meta        // channel: path to metadata CSV (--meta_file)

    main:

    ch_versions = channel.empty()

    //
    // Step 1: Log software versions for reproducibility
    //
    LOG_VERSIONS()
    ch_versions = ch_versions.mix(LOG_VERSIONS.out.versions)

    //
    // Step 2: Import IDAT files and build RGset / targets
    //
    IMPORT_IDAT(ch_samplesheet, ch_idatdir)
    ch_versions = ch_versions.mix(IMPORT_IDAT.out.versions)

    //
    // Step 3: Raw intensity QC (on imported RGset)
    //
    RAW_INTENSITIES_QC(IMPORT_IDAT.out.rg_rds)
    ch_versions = ch_versions.mix(RAW_INTENSITIES_QC.out.versions)

    //
    // Step 4: SeSAMe preprocessing + pOOBAH QC
    // Produces filtered RGset and normalized beta values
    //
    SESAME_POOBAH_QC(
        IMPORT_IDAT.out.rg_rds,
        IMPORT_IDAT.out.targets_tsv,
        ch_idatdir
    )
    ch_versions = ch_versions.mix(SESAME_POOBAH_QC.out.versions)

    //
    // Step 5: Control probes QC report
    //
    QC_REPORT_CONTROLS(SESAME_POOBAH_QC.out.rg_filt_rds)
    ch_versions = ch_versions.mix(QC_REPORT_CONTROLS.out.versions)

    //
    // Step 6: SNP-based sample relationship / identity heatmap
    //
    SNP_HEATMAP(SESAME_POOBAH_QC.out.rg_filt_rds)
    ch_versions = ch_versions.mix(SNP_HEATMAP.out.versions)

    //
    // Step 7: Beta value density plots (raw + SeSAMe normalized)
    //
    DENSITY_PLOTS(
        SESAME_POOBAH_QC.out.rg_filt_rds,
        SESAME_POOBAH_QC.out.beta_norm_rds
    )
    ch_versions = ch_versions.mix(DENSITY_PLOTS.out.versions)

    //
    // Step 8: Sex prediction QC
    //
    SEX_QC(
        SESAME_POOBAH_QC.out.rg_filt_rds,
        SESAME_POOBAH_QC.out.targets_filt_tsv,
        ch_meta
    )
    ch_versions = ch_versions.mix(SEX_QC.out.versions)

    //
    // Step 9: Apply independent filters (bead count, chrX/Y, non-CG).
    // rg_filt_rds is the post-sample-QC RGChannelSet from SESAME_POOBAH_QC,
    // used by minfi::getNBeads() for the bead count filter.
    //
    FILTER_XY_NONCG_BEADS(
        SESAME_POOBAH_QC.out.rg_filt_rds,
        SESAME_POOBAH_QC.out.beta_norm_rds
    )
    ch_versions = ch_versions.mix(FILTER_XY_NONCG_BEADS.out.versions)

    //
    // Step 10: Optional cell composition estimation
    // When disabled, an empty channel is passed so ALIGN_META omits cell
    // composition columns from the metadata and downstream design matrices.
    //
    def ch_cell_counts
    if (params.do_estimate_cellcomp) {
        CELL_COMP(SESAME_POOBAH_QC.out.rg_filt_rds)
        ch_cell_counts = CELL_COMP.out.cell_counts_csv
        ch_versions = ch_versions.mix(CELL_COMP.out.versions)
    } else {
        ch_cell_counts = Channel.value([])
    }

    //
    // Step 11: Split / collapse probes (DMP vs DMR beta matrices)
    //
    SPLIT_COLLAPSE(FILTER_XY_NONCG_BEADS.out.beta3_rds)
    ch_versions = ch_versions.mix(SPLIT_COLLAPSE.out.versions)

    //
    // Step 12: Final matrix exports (bVals / mVals)
    //
    FINAL_EXPORTS(SPLIT_COLLAPSE.out.beta_clean_rds)
    ch_versions = ch_versions.mix(FINAL_EXPORTS.out.versions)

    //
    // Step 13: Align metadata with methylation matrices
    //
    ALIGN_META(
        FINAL_EXPORTS.out.bvals_rds,
        FINAL_EXPORTS.out.mvals_rds,
        SESAME_POOBAH_QC.out.targets_filt_tsv,
        ch_meta,
        ch_cell_counts,
        SEX_QC.out.qc_pred_sex_csv
    )
    ch_versions = ch_versions.mix(ALIGN_META.out.versions)

    //
    // Step 14: Optional covariate PCA (controlled via ext.when in modules.config)
    //
    COVARIATE_PCA(
        ALIGN_META.out.bvals_aligned_rds,
        ALIGN_META.out.meta_aligned_csv
    )
    ch_versions = ch_versions.mix(COVARIATE_PCA.out.versions)

    //
    // Step 15: Optional DMP analysis — limma (controlled via ext.when in modules.config)
    //
    DMP_LIMMA(
        ALIGN_META.out.mvals_aligned_rds,
        ALIGN_META.out.bvals_aligned_rds,
        ALIGN_META.out.meta_aligned_csv
    )
    ch_versions = ch_versions.mix(DMP_LIMMA.out.versions)

    //
    // Step 16: Optional DMR analysis — DMRcate (controlled via ext.when in modules.config)
    //
    DMR_DMRCATE(
        SPLIT_COLLAPSE.out.beta_dmr_rds,
        SESAME_POOBAH_QC.out.targets_filt_tsv,
        ALIGN_META.out.meta_aligned_csv
    )
    ch_versions = ch_versions.mix(DMR_DMRCATE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_methylarray_software_versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
