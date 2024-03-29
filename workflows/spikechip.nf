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
include { PICARD_MARKDUPLICATES       } from '../modules/nf-core/picard/markduplicates/main'
include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'
include { DEEPTOOLS_BAMCOVERAGE       } from '../modules/nf-core/deeptools/bamcoverage/main'
include { DEEPTOOLS_BAMCOVERAGE   as  DEEPTOOLS_BAMCOVNOCALIB  } from '../modules/nf-core/deeptools/bamcoverage/main'
include { DEEPTOOLS_BAMCOVERAGE   as  DEEPTOOLS_BAMCOVSCALING  } from '../modules/nf-core/deeptools/bamcoverage/main' 
include { BEDTOOLS_MERGE              } from '../modules/nf-core/bedtools/merge/main'
include { SAMTOOLS_FLAGSTAT       as  SAMTOOLS_DOWSAMPFLAGSTAT  } from '../modules/nf-core/samtools/flagstat/main.nf'
include { SAMTOOLS_BAM2FQ             } from '../modules/nf-core/samtools/bam2fq/main'

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
    SAMTOOLS_FAIDX.out.fai.map{meta,path -> [path]}.set{faidx_path}

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

        TRIMGALORE (
            INPUT_CHECK.out.reads  
        )

        FASTQC (
            TRIMGALORE.out.reads
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())


        //INPUT_CHECK.out.reads.view()

        ch_bowtie2_index = [ [:], file(params.bowtie2_index) ]

        BOWTIE2_ALIGN (
            TRIMGALORE.out.reads,
            ch_bowtie2_index,
            false,
            "sort"
        )

        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

        ch_bamfiles=BOWTIE2_ALIGN.out.aligned

        //BOWTIE2_ALIGN.out.aligned.view()
                
    } else {
    
    //input start channel   

    ch_bamfiles=channel.fromPath(params.input)
        | splitCsv (header: true)
        | map { row -> 
            ssinfo = row.subMap ('id', 'single_end','condition','details','analysis','bam')
            [[id:ssinfo.id, single_end:ssinfo.single_end, condition:ssinfo.condition, details:ssinfo.details, analysis:ssinfo.analysis], ssinfo.bam]
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

        SAMTOOLS_BAM2FQ (
            SAMTOOLS_SPLITSPECIES.out.refpath,
            "TRUE"
        )

        //ch_meta_hub=SAMTOOLS_FLAGSTAT.out.flagstat.map{meta,path -> meta}.collect()
        //ch_path_hub=SAMTOOLS_FLAGSTAT.out.flagstat.map{meta,path -> path}.collect()

        //SAMTOOLS_FLAGSTAT.out.reference.groupTuple{it.analysis}.view()

        SAMTOOLS_FLAGSTAT.out.reference.map{meta,path -> 
         an=meta.subMap('analysis')
         [an.analysis, meta,path]
        }.groupTuple()
        .map {an, meta, path -> 
         [meta,path] 
        }
        .set{ch_refstat_by_analysis}

        //ch_refstat_by_analysis.view()

        SAMTOOLS_FLAGSTAT.out.spikein
            .map{meta,path -> 
                an=meta.subMap('analysis')
                [an.analysis, meta,path]
            }.groupTuple()
                .map {an, meta, path -> 
                [meta,path] 
            }
            .set{ch_spikestat_by_analysis}

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
        
        //ch_reference_hub.view()
        //ch_spikein_hub.view{"spikein: ${it}"}

        //split by analysis

        //TODO filter analysis replicates without input

        if (params.calibStrategy != "no_calibration" ){
            CALCULATEDOWNFACTOR (
            ch_reference_hub,
            ch_spikein_hub
        )
        }

        //TODO filter analysis replicates without input

        //CALCULATEDOWNFACTOR.out.downfile.view()

        //START OF DOWN SAMPLING WITH INPUT
        if (params.calibStrategy == "with_input"){

            downfl=CALCULATEDOWNFACTOR.out.calibfile.flatten() //devi mantenere i progetti separati
            //downfl.view()

            ch_downfact=downfl
                | splitCsv (header: true)
                | map { row -> 
                    downss = row.subMap ('id', 'single_end','condition','details','analysis','downfactor')
                    [[id:downss.id, single_end:downss.single_end, condition:downss.condition, details:downss.details, analysis:downss.analysis, downfactor:downss.downfactor]]
                }

            //ch_downfact.flatten().view()
            //SAMTOOLS_SPLITSPECIES.out.bam.view()

            ch_downfact.flatten()
            .map{ meta -> 
                id=meta.id
                [id, meta]
            }
            .set{ch_downfc_by_id}

            SAMTOOLS_SPLITSPECIES.out.refpath.map{meta, ref_bam, ref_bai ->
                id=meta.subMap('id')
                [id.id, ref_bam, ref_bai]
            }set{ch_bam_by_id}

            ch_downfc_by_id
                .combine(ch_bam_by_id, by:0)
                .map{id,meta,bam, bai -> 
                [meta,bam, bai]
                }
                .set{ch_downfc}

            //ch_downin.filter{it[0].downfactor < 1}.view()
            ch_downin=ch_downfc.filter{it[0].downfactor.toFloat() < 1}.map{meta,bam,bai->[meta,bam]}
            ch_nodown=ch_downfc.filter{it[0].downfactor.toFloat() >= 1}//.map{meta,bam,bai->[meta,bam]}

            //ch_downin.view()
            //ch_nodown.view()

            SAMTOOLS_DOWNSAMPLING(
                ch_downin
            )

            SAMTOOLS_DOWSAMPFLAGSTAT(SAMTOOLS_DOWNSAMPLING.out.bam)


            //SAMTOOLS_DOWNSAMPLING.out.bam.view{"SAMTOOLS_DOWNSAMPLING.out : ${it}"}

            ch_tomerge=SAMTOOLS_DOWNSAMPLING.out.bam
                                .mix(ch_nodown)
                                .map{meta, path, path2 -> 
                                                [meta,path]
                                                }        

            ch_tomerge.map{ meta, path -> 
                            condition=meta.subMap('condition')
                            analysis=meta.subMap('analysis')
                            meta=meta
                            path=path
                            [condition.condition + "_" + analysis.analysis, meta, path]
                        }
                    .groupTuple()
                    
                    .map{condition, meta, path ->
                        [[id:meta.analysis[0] + "_" + meta.condition[0] , single_end:meta.single_end[0]], path]
                    }.set{ch_mergeBam}

            //ch_mergeBam.view{"ch_mergeBam : ${it}"}

            SAMTOOLS_MERGE (
                ch_mergeBam
            )

            //SAMTOOLS_MERGE.out.bam.view{"SAMTOOLS_MERGE.out.bam : ${it}"}
            //ch_fasta_meta.view{"ch_fasta_meta : ${it}"}

            //SAMTOOLS_FAIDX.out.fai.view{"SAMTOOLS_FAIDX.out.fai : ${it}"}

            ch_tocov=SAMTOOLS_MERGE.out.bam.mix(SAMTOOLS_DOWNSAMPLING.out.bam) 
            ch_todeepcoverage=ch_tocov.mix(ch_nodown)

            //ch_tocov.view()
            //ch_todeepcoverage.view()
            
            //SAMTOOLS_MERGE.out.bam.view()

            DEEPTOOLS_BAMCOVERAGE (
                ch_todeepcoverage,
                params.fasta,
                faidx_path
            )

        }
        //END OF DOWNSAMPLING WITH INPUT

        //START OF DOWNSAMPLING WITHOUT INPUT

        if (params.calibStrategy == "no_calibration") {

            //SAMTOOLS_SPLITSPECIES.out.refpath.view()

            DEEPTOOLS_BAMCOVNOCALIB ( 
                SAMTOOLS_SPLITSPECIES.out.refpath,
                params.fasta,
                faidx_path  
            )
        }

        if (params.calibStrategy == "without_input") {
            
            noin=CALCULATEDOWNFACTOR.out.calibfile.flatten()
            //noin.view()

            ch_noin=noin
                | splitCsv (header: true)
                | map { row -> 
                    noinss = row.subMap ('id', 'single_end','condition','details','analysis','noinputNorm')
                    [[id:noinss.id, single_end:noinss.single_end, condition:noinss.condition, details:noinss.details, analysis:noinss.analysis, noinputNorm:noinss.noinputNorm]]
                }

            //ch_noin.view()

            ch_noin.flatten()
            .map{ meta -> 
                id=meta.id
                [id, meta]
            }
            .set{ch_noin_by_id}

        //ch_downin.view()
        //SAMTOOLS_SPLITSPECIES.out.refpath.view()

        SAMTOOLS_SPLITSPECIES.out.refpath.map{meta, bampath, baipath ->
            id=meta.subMap('id')
            [id.id, bampath, baipath]
        }set{ch_id_bam_downin}

        //ch_id_bam_downin.view()

        //ch_id_bam_downin.combine(ch_noin_by_id, by:0).view()

        ch_noin_by_id
            .combine(ch_id_bam_downin, by:0)
            .map{id,meta,noinbam, noinbai -> 
            [meta, noinbam, noinbai]
            }
            .set{ch_cov_scaling}

        //ch_cov_scaling.view()
        

        DEEPTOOLS_BAMCOVSCALING ( //scaling parameter is in module config file
            ch_cov_scaling,
            params.fasta,
            faidx_path
        )
     
        }

        //END OF DOWNSAMPLING WITHOUT INPUT
 
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

    if (!params.onlyBAM && params.calibStrategy == "with_input") {
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_DOWSAMPFLAGSTAT.out.flagstat.collect{it[1]}.ifEmpty([]))
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
