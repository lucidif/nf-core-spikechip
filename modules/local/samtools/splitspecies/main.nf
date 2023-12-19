process SAMTOOLS_SPLITSPECIES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(input)
    val   ref_name
    val   spikein_name
    //tuple val(meta2), path(fasta)


    output:
    //tuple val(meta), path("*_${ref_name}.bam"), path("*_${ref_name}.bam.bai") ,  emit: refbam
    //tuple val(meta), path("*_${spikein_name}.bam"), path("*_${spikein_name}.bam.bai") ,  emit: spikebam
    tuple val(meta), path("*_${ref_name}.bam"), path("*_${spikein_name}.bam") ,  emit: bam
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    //def args2 = task.ext.args2 ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    //def reference = fasta ? "--reference ${fasta}" : ""
    //def readnames = qname ? "--qname-file ${qname}": ""
    // def file_type = args.contains("--output-fmt sam") ? "sam" :
    //                 args.contains("--output-fmt bam") ? "bam" :
    //                 args.contains("--output-fmt cram") ? "cram" :
    //                 input.getExtension()

    """
    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        -h \\
        ${input} | grep -v ${spikein_name} | sed s/${ref_name}\\_chr/chr/g | samtools view -bhS - > ${meta.id}_${ref_name}.bam

    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        -h \\
        ${input} | grep -v ${ref_name} | sed s/${spikein_name}\\_chr/chr/g | samtools view -bhS - > ${meta.id}_${spikein_name}.bam 

    samtools index -M *.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:

    """
    touch ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
