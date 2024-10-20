/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// methylarray
include { PREPROCESS                    } from '../modules/local/preprocess/main'
include { XREACTIVE_PROBES_FIND_REMOVE  } from '../modules/local/xreactive_probes_find_remove/main'
include { REMOVE_SNP_PROBES             } from '../modules/local/remove_snp_probes/main'
include { REMOVE_SEX_CHROMOSOMES        } from '../modules/local/remove_sex_chromosomes/main'
include { REMOVE_CONFOUNDING_PROBES     } from '../modules/local/remove_confounding_probes/main'
include { ADJUST_CELL_COMPOSITION       } from '../modules/local/adjust_cell_composition/main'
include { ADJUST_BATCH_EFFECT           } from '../modules/local/adjust_batch_effect/main'
include { FIND_DMP                      } from '../modules/local/find_dmp/main'
include { FIND_DMR                      } from '../modules/local/find_dmr/main'
include { FIND_BLOCKS                   } from '../modules/local/find_blocks/main'


// nf-core
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
    extensive_metadata = params.sample_metadata ? Channel.fromPath(params.sample_metadata) : Channel.empty()

    //
    // MODULE: Run PREPROCESS
    //
    PREPROCESS (
        ch_samplesheet
    )

    //
    // MODULE: Run XREACTIVE_PROBES_FIND_REMOVE
    //
    // Download from: https://github.com/pjhop/DNAmCrosshyb at https://doi.org/10.5281/zenodo.4088019 and point to the version
    // TODO: If bs_genome_path is not provided then the pipeline might resolve it itself
    genome_path = Channel.fromPath(params.bs_genome_path)
    XREACTIVE_PROBES_FIND_REMOVE (
        PREPROCESS.out.rdata,
        genome_path
    )

    //
    // MODULE: Run REMOVE_SNP_PROBES
    //
    REMOVE_SNP_PROBES (
        XREACTIVE_PROBES_FIND_REMOVE.out.rdata
    )

    //
    // Optional steps of methylarray
    //
    
    //
    // Output channel following optional steps
    //
    final_bVals_ch = Channel.empty()
    current_bVals_ch = Channel.empty()

    if (params.run_optional_steps) {
        current_bVals_ch = REMOVE_SNP_PROBES.out.csv_bVals
        if (params.remove_sex_chromosomes || params.remove_confounding_probes) { // If params.remove_confounding_probes then this has to be run
            current_bVals_ch = REMOVE_SNP_PROBES.out.rdata
            //
            // MODULE: Run REMOVE_SEX_CHROMOSOMES
            //
            REMOVE_SEX_CHROMOSOMES (
                current_bVals_ch,
                PREPROCESS.out.rdata_rgSet
            )
        }

        if (params.remove_confounding_probes || params.remove_sex_chromosomes) { // Currently depends on REMOVE_SEX_CHROMOSOMES
            //
            // MODULE: Run REMOVE_CONFOUNDING_PROBES
            //
            REMOVE_CONFOUNDING_PROBES (
                REMOVE_SEX_CHROMOSOMES.out.mVals_csv,
                REMOVE_SEX_CHROMOSOMES.out.bVals_csv,
                REMOVE_SEX_CHROMOSOMES.out.mSetSqFlt,
                extensive_metadata
            )
            current_bVals_ch = REMOVE_CONFOUNDING_PROBES.out.bVals
        }

        if (params.adjust_cell_composition) {
            //
            // MODULE: Run ADJUST_CELL_COMPOSITION
            //
            ADJUST_CELL_COMPOSITION (
                current_bVals_ch
            )
            current_bVals_ch = ADJUST_CELL_COMPOSITION.out.bVals
        }

        if (params.adjust_batch_effect) {
            //
            // MODULE: Run REMOVE_CONFOUNDING_PROBES
            //
            ADJUST_BATCH_EFFECT (
                current_bVals_ch,
                extensive_metadata
            )
            current_bVals_ch = ADJUST_BATCH_EFFECT.out.bVals
        }
    }

    //
    // Update final bVals channel
    //
    final_bVals_ch = current_bVals_ch

    //
    // MODULE: Run FIND_DMP
    //
    FIND_DMP (
        final_bVals_ch,
        extensive_metadata
    )

    //
    // MODULE: Run FIND_DMR
    //
    if (params.find_dmrs) { // Will not be able to find DMRs with test data
        FIND_DMR (
            final_bVals_ch,
            extensive_metadata
        )
    }
    //
    // MODULE OPTIONAL: Run FIND_BLOCKS
    //
    if (params.find_blocks) {
        //
        // MODULE: Run FIND_BLOCKS
        //
        FIND_BLOCKS (
            final_bVals_ch,
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
