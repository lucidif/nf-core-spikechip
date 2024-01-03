process SAMTOOLS_DOWNSAMPLING {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(input)//, path(index)
    //tuple val(meta2), path(fasta)
    //path qname

    output:
    tuple val(meta), path("*.bam"),  emit: bam
    //tuple val(meta), path("*.cram"), emit: cram,    optional: true
    //tuple val(meta), path("*.sam"),  emit: sam,     optional: true
    //tuple val(meta), path("*.bai"),  emit: bai,     optional: true
    //tuple val(meta), path("*.csi"),  emit: csi,     optional: true
    //tuple val(meta), path("*.crai"), emit: crai,    optional: true
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def inname = "${meta.id}_${params.reference_genome}.bam"
    def outname = "${meta.id}_ds"
    def downFactor = meta.downfactor
    //def reference = fasta ? "--reference ${fasta}" : ""
    //def readnames = qname ? "--qname-file ${qname}": ""
    // def file_type = args.contains("--output-fmt sam") ? "sam" :
    //                 args.contains("--output-fmt bam") ? "bam" :
    //                 args.contains("--output-fmt cram") ? "cram" :
    //                 input.getExtension()
    //if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        -s ${downFactor} \\
        --subsample-seed 123 \\
        --subsample ${downFactor} \\
        -b \\
        -o pgin.bam \\
        ${input}

    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        -H \\
        pgin.bam > pg.header.sam

    awk '!/^@PG/' pg.header.sam > nopg.header.sam

    samtools \\
        reheader \\
        nopg.header.sam \\
        pgin.bam > ${outname}.bam
   

    rm pgin.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    // def file_type = args.contains("--output-fmt sam") ? "sam" :
    //                 args.contains("--output-fmt bam") ? "bam" :
    //                 args.contains("--output-fmt cram") ? "cram" :
    //                 input.getExtension()
    //if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    //def index = args.contains("--write-index") ? "touch ${prefix}.csi" : ""

    """
    touch ${outname}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
