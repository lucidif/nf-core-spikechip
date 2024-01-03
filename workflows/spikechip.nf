/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.fasta) { ch_fasta =  Channel.fromPath(params.fasta) } else { exit 1, 'Fasta reference genome not specified!' }

// Modify fasta channel to include meta data
ch_fasta_meta = ch_fasta.map{ it -> [[id:it[0].baseName], it] }.collect()

/*

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowSpikechip.initialise(params, log)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
//include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { BOWTIE2_ALIGN               } from '../modules/nf-core/bowtie2/align/main.nf'
//include { SAMBAMBA_MARKDUP            } from '../modules/nf-core/sambamba/markdup/main'
include { PICARD_MARKDUPLICATES       } from '../modules/nf-core/picard/markduplicates/main'
include { TRIMMOMATIC                 } from '../modules/nf-core/trimmomatic/main'
include { DEEPTOOLS_BAMCOVERAGE       } from '../modules/nf-core/deeptools/bamcoverage/main' 
include { BEDTOOLS_MERGE } from '../modules/nf-core/bedtools/merge/main'

include { SAMTOOLS_FAIDX              } from '../modules/local/samtools/faidx/main'
include { SAMTOOLS_SPLITSPECIES       } from '../modules/local/samtools/splitspecies/main.nf'
include { SAMTOOLS_FLAGSTAT           } from '../modules/local/samtools/flagstat/main.nf'
include { CALCULATEDOWNFACTOR         } from '../modules/local/calculatedownfactor/main.nf'
include { SAMTOOLS_DOWNSAMPLING       } from '../modules/local/samtools/downsampling/main.nf'
include { SAMTOOLS_MERGE              } from '../modules/local/samtools/merge/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SPIKECHIP {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    SAMTOOLS_FAIDX (
            ch_fasta_meta,
            [[], []]
        )

    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    if (!params.fromBAM) {

        INPUT_CHECK (
            file(params.input)
        )
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
        // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
        // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
        // ! There is currently no tooling to help you write a sample sheet schema

        //INPUT_CHECK.out.reads.view()

        //
        // MODULE: Run FastQC
        //

        TRIMMOMATIC (
            INPUT_CHECK.out.reads  
        )

        FASTQC (
            TRIMMOMATIC.out.trimmed_reads
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())


        //INPUT_CHECK.out.reads.view()

        ch_bowtie2_index = [ [:], file(params.bowtie2_index) ]

        BOWTIE2_ALIGN (
            TRIMMOMATIC.out.trimmed_reads,
            ch_bowtie2_index,
            false,
            "sort"
        )

        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

        ch_bamfiles=BOWTIE2_ALIGN.out.aligned

        //SAMBAMBA_MARKDUP (BOWTIE2_ALIGN.out.aligned)

        //BOWTIE2_ALIGN.out.aligned.view()
                
    } else {
     //input start channel   
     //[[id:SPT5_INPUT_REP2_T1, single_end:false], /mnt/c/wkdir/lucio/test_nf_spikeinchip/work/f3/67c45d925563552cb242de816b4fa3/SPT5_INPUT_REP2_T1.bam]

    ch_bamfiles=channel.fromPath(params.input)
        | splitCsv (header: true)
        | map { row -> 
            ssinfo = row.subMap ('id', 'single_end','condition','details','bam')
            [[id:ssinfo.id, single_end:ssinfo.single_end, condition:ssinfo.condition, details:ssinfo.details], ssinfo.bam]
        }

    //ch_bamfiles.view()    

    }

    PICARD_MARKDUPLICATES (
        ch_bamfiles,
        ch_fasta_meta,
        SAMTOOLS_FAIDX.out.fai.collect()
        )

    if (!params.onlyBAM) {


        SAMTOOLS_SPLITSPECIES (
            PICARD_MARKDUPLICATES.out.bam,
            params.reference_genome,
            params.spikein_genome 
        )

    ch_versions = ch_versions.mix(SAMTOOLS_SPLITSPECIES.out.versions.first())
            
        //SAMTOOLS_SPLITSPECIES.out.refbam.view()
        SAMTOOLS_FLAGSTAT (
            SAMTOOLS_SPLITSPECIES.out.bam
        )
        //SAMTOOLS_FLAGSTAT.out.flagstat.view()

        //ch_meta_hub=SAMTOOLS_FLAGSTAT.out.flagstat.map{meta,path -> meta}.collect()
        //ch_path_hub=SAMTOOLS_FLAGSTAT.out.flagstat.map{meta,path -> path}.collect()

        ch_reference_hub=SAMTOOLS_FLAGSTAT.out.reference.flatten().collect()
                            .map{ tuplain -> 
                            def evenTuple = tuplain.findAll { tuplain.indexOf(it) % 2 == 0 }
                            def oddTuple = tuplain.findAll { tuplain.indexOf(it) % 2 == 1 }

                            [evenTuple,oddTuple]

                             }

        ch_spikein_hub=SAMTOOLS_FLAGSTAT.out.spikein.flatten().collect()
                            .map{ tuplain -> 
                            def evenTuple = tuplain.findAll { tuplain.indexOf(it) % 2 == 0 }
                            def oddTuple = tuplain.findAll { tuplain.indexOf(it) % 2 == 1 }

                            [evenTuple,oddTuple]

                             }                     
        
        ch_meta_collect=SAMTOOLS_FLAGSTAT.out.reference.map{meta, path -> 
            newmeta=meta.collect()
        }

 

        //ch_reference_hub.view()
        //ch_spikein_hub.view()
  
        CALCULATEDOWNFACTOR (
            ch_reference_hub,
            ch_spikein_hub
        )

        downfl=CALCULATEDOWNFACTOR.out.downfile.flatten()
        //downfl.view()

        ch_downfact=downfl
            | splitCsv (header: true)
            | map { row -> 
                downss = row.subMap ('id', 'single_end','condition','details','downfactor')
                [[id:downss.id, single_end:downss.single_end, condition:downss.condition, details:downss.details, downfactor:downss.downfactor]]
            }

        //ch_downfact.flatten().view()
        //SAMTOOLS_SPLITSPECIES.out.bam.view()

        ch_downfact.flatten()
        .map{ meta -> 
            id=meta.id
            [id, meta]
        }
        .set{ch_downfc_by_id}

        SAMTOOLS_SPLITSPECIES.out.bam.map{meta, ref_bam, spike_bam ->
            id=meta.subMap('id')
            [id.id, ref_bam]
        }set{ch_bam_by_id}

        ch_downfc_by_id
            .combine(ch_bam_by_id, by:0)
            .map{id,meta,bam -> 
            [meta,bam]
            }
            .set{ch_downfc}

        //ch_downin.filter{it[0].downfactor < 1}.view()
        ch_downin=ch_downfc.filter{it[0].downfactor.toFloat() < 1}
        ch_nodown=ch_downfc.filter{it[0].downfactor.toFloat() >= 1}

        SAMTOOLS_DOWNSAMPLING(
            ch_downin
        )

        ch_tomerge=SAMTOOLS_DOWNSAMPLING.out.bam.mix(ch_nodown)

        //ch_tomerge.view()

        ch_tomerge.map{ meta, path -> 
                        condition=meta.subMap('condition')
                        meta=meta
                        path=path
                        [condition.condition, meta, path]
                      }
                .groupTuple()
                //.flatten()
                .map{condition, meta, path->
                    //meta.id=meta.condition
                    //id=meta.subMap('condition')
                    //single_end=meta.subMap('single_end')
                    [[id:meta.condition[0], single_end:meta.single_end[0]], path]
                }.set{ch_mergeBam}

        //ch_mergeBam.view()

        SAMTOOLS_MERGE (
            ch_mergeBam
        )

        SAMTOOLS_MERGE.out.bam.view()

        // DEEPTOOLS_BAMCOVERAGE (
        //     SAMTOOLS_MERGE.out.bam,
        //     ch_fasta_meta,
        //     SAMTOOLS_FAIDX.out.fai.collect()
        // )


 

    }        

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSpikechip.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSpikechip.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    
    
    if (!params.fromBAM) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
    }

    if (!params.onlyBAM) {
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.reference.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.spikein.collect{it[1]}.ifEmpty([]))
    }

    
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
