
include { SAMTOOLS_VIEW } from '../modules/nf-core/samtools/view/main'
include { SAMBAMBA_MARKDUP } from '../modules/nf-core/sambamba/markdup/main' 

workflow SPLIT_SPIKEIN { 

    take:
        ch_mergeref_bam
        ref_name
        spikein_name

    main:

        //1)sambamba markdup --tmpdir /data/tmp/ -r -t 56 $filename\_UniqMapped_sorted.bam $filename\_UniqMapped_sorted_rmdup.bam

        SAMBAMBA_MARKDUP (ch_mergeref_bam)

        //2)samtools view -h $filename\_UniqMapped_sorted_rmdup.bam | grep -v $spikegenome | sed s/$genome\_chr/chr/g | samtools view -bhS - > $filename\_$genome.UniqMapped_sorted_rmdup.bam
        
        SAMTOOLS_VIEW ( )
        

        //3)samtools view -h $filename\_UniqMapped_sorted_rmdup.bam | grep -v $genome | sed s/$spikegenome\_chr/chr/g | samtools view -bhS - > $filename\_$spikegenome.UniqMapped_sorted_rmdup.bam



}