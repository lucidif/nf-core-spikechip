process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(input_files, stageAs: "?/*")
    //tuple val(meta2), path(fasta)
    //tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}_mg.bam"), path("${prefix}_mg.bam.bai") ,  emit: bam
    //tuple val(meta), path("${prefix}.cram"), optional:true, emit: cram
    //tuple val(meta), path("*.csi")         , optional:true, emit: csi
    //tuple val(meta), path("*.crai")        , optional:true, emit: crai
    path  "versions.yml"                                  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    //def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    //def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        view \\
        -H \\
        $input_files > header.sam

    samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        $args \\
        ${prefix}_mg.bam \\
        $input_files

    samtools \\
        sort \\
        --threads ${task.cpus-1} \\
        -o ${prefix}_sort_mg.bam \\
        ${prefix}_mg.bam

    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        ${prefix}_mg.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args   ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
   //def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    //def index_type = file_type == "bam" ? "csi" : "crai"
    //def index = args.contains("--write-index") ? "touch ${prefix}.${index_type}" : ""
    """
    touch ${prefix}_mg.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
