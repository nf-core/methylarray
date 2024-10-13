/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESS             } from '../modules/local/preprocess/main'
include { XREACTIVE_PROBES_FIND_REMOVE } from '../modules/local/xreactive_probes_find_remove/main'
include { REMOVE_SNP_PROBES      } from '../modules/local/remove_snp_probes/main'
include { REMOVE_SEX_CHROMOSOMES      } from '../modules/local/remove_sex_chromosomes/main'
include { REMOVE_CONFOUNDING_PROBES      } from '../modules/local/remove_confounding_probes/main'
include { ADJUST_CELL_COMPOSITION      } from '../modules/local/adjust_cell_composition/main'
include { ADJUST_BATCH_EFFECT     } from '../modules/local/adjust_batch_effect/main'
include { FIND_DMP     } from '../modules/local/find_dmp/main'

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_methylarray_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLARRAY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_preprocessed_files = Channel.empty()
    extensive_metadata = params.confounding_probes_rm_sheet ? Channel.fromPath(params.confounding_probes_rm_sheet) : Channel.empty()

    //
    // MODULE: Run PREPROCESS
    //
    PREPROCESS (
        ch_samplesheet
    )
    ch_preprocessed_files = ch_preprocessed_files.mix(PREPROCESS.out.rdata)

    //
    // MODULE: Run XREACTIVE_PROBES_FIND_REMOVE
    //
    // Download from: https://github.com/pjhop/DNAmCrosshyb at https://doi.org/10.5281/zenodo.4088019
    genome_path = Channel.fromPath("${projectDir}/results/genome_bs/hg19")
    XREACTIVE_PROBES_FIND_REMOVE (
        PREPROCESS.out.rdata,
        genome_path
    )
    ch_preprocessed_files = ch_preprocessed_files.mix(XREACTIVE_PROBES_FIND_REMOVE.out.rdata)

    //
    // MODULE: Run REMOVE_SNP_PROBES
    //
    REMOVE_SNP_PROBES (
        XREACTIVE_PROBES_FIND_REMOVE.out.rdata
    )
    ch_preprocessed_files = ch_preprocessed_files.mix(REMOVE_SNP_PROBES.out.rdata)

    //
    // MODULE: Run REMOVE_SNP_PROBES
    //
    REMOVE_SEX_CHROMOSOMES (
        REMOVE_SNP_PROBES.out.rdata,
        PREPROCESS.out.rdata_rgSet
    )

    //
    // MODULE: Run REMOVE_CONFOUNDING_PROBES
    // NOTE: This is not completly integrated as additional insights are needed in relation to the extensive_metadata.csv file
    //
    if (params.confounding_probes_rm_sheet) {
        REMOVE_CONFOUNDING_PROBES (
            REMOVE_SNP_PROBES.out.csv_mVals,
            REMOVE_SNP_PROBES.out.csv_bVals,
            REMOVE_SNP_PROBES.out.rdata,
            extensive_metadata
        )
    }

    //
    // MODULE: Run REMOVE_SNP_PROBES
    // NOTE: Probably failes due to the smaller input file size than expected
    //
    if (params.adjust_cell_composition) {
        ADJUST_CELL_COMPOSITION (
            REMOVE_SNP_PROBES.out.csv_bVals,
        )
    }

    //
    // MODULE: Run REMOVE_CONFOUNDING_PROBES
    // NOTE: This is not completly integrated as additional insights are needed in relation to the extensive_metadata.csv file
    //
    if (params.confounding_probes_rm_sheet) {
        ADJUST_BATCH_EFFECT (
            REMOVE_SNP_PROBES.out.csv_bVals,
            extensive_metadata
        )
    }

    //
    // MODULE: Run REMOVE_CONFOUNDING_PROBES
    // NOTE: This is not completly integrated as additional insights are needed in relation to the extensive_metadata.csv file
    //
    if (params.confounding_probes_rm_sheet) {
        FIND_DMP (
            REMOVE_SNP_PROBES.out.csv_bVals,
            extensive_metadata
        )
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
