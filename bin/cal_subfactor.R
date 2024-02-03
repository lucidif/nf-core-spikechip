#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

in_meta <- args[1]
scale_fact<- as.numeric(as.character(args[2]))
calib_c <- as.numeric(as.character(args[3]))
mode <- args[4] #"with_input" "without_input" TODO add calibration mode selector

g1<-gsub("\\[\\[",", \\[",in_meta)
g2<-gsub(", \\[","?",g1)
g3<-gsub(", \\[","?",g2)
g4<-gsub("\\]\\?","?",g3)
g5<-gsub("\\]\\]","?",g4)

g6<-strsplit(g5, "\\?")[[1]]

meta_id<-grep("id\\:",g6)
#path_id<-meta_id+1
meta_unformat<-g6[meta_id]

#estrei in nomi delle colonne
univSep<-gsub("\\:","\\?",meta_unformat[1])
univSep<-gsub("\\, ","\\?",univSep)
col_noformat<-strsplit(univSep,"\\?")[[1]]
col_names<-col_noformat[seq(1, length(col_noformat), by = 2)]

meta.format<-matrix(nrow=length(meta_unformat), ncol=length(col_names))
colnames(meta.format)<-col_names
#formattiamo i meta
for (i in 1:length(meta_unformat)){
  meta.split<-strsplit(meta_unformat[i],"\\, ")[[1]]
  for(j in 1:length(col_names)){
    meta.format[i,j]<-gsub(paste0(col_names[j],"\\:"),"",meta.split[j])
  }
  
}

#formatta i path
#path.unformat<-g4[path_id]
#path.unformat<-gsub(" ","",path.unformat)
#path.format=matrix(ncol=2,nrow=length(path.unformat))
#colnames(path.format)<-c("ref_file","spikein_file")

# for(i in 1:length(path.unformat)){
#   
#   paths<-strsplit(path.unformat[i],",")[[1]]
#   path.format[i,1]<-paths[1]
#   path.format[i,2]<-paths[2]
#   
# }

#path.format<-gsub("\\]","",path.format)
#path.format[,2]<-gsub("\\]","",path.format[,2])

#inmatrix<-cbind(meta.format, path.format)

save(meta.format, file="inputs.RData")
save(in_meta,scale_fact, calib_c, mode,file="parameters.RData")

input.table<-as.data.frame(meta.format)


#colnames=("id", "single_end", "condition", "details", "analysis", "calib.eq", "perc_spikein_in_sample", "perc_spikein_in_input", "scaling_factor")

ref_aln_reads<-c()
spikein_aln_reads<-c()

for (i in 1:length(input.table$id)){
    ref_file=paste0(input.table$id[i],"_ref.flagstat")
    spikein_file=paste0(input.table$id[i],"_spike.flagstat")
    
    ref_stats<-read.delim(ref_file,sep="\n",header=FALSE)
    ref_aln_reads[i]<-as.numeric(strsplit(ref_stats[7,1]," \\+ ")[[1]][1])
    
    spikein_stats<-read.delim(spikein_file,sep="\n",header=FALSE)
    spikein_aln_reads[i]<-as.numeric(strsplit(spikein_stats[7,1]," \\+ ")[[1]][1])

  }

  reads.table<-data.frame(id=input.table$id, single_end=input.table$single_end, condition=input.table$condition, details=input.table$details,analysis=input.table$analysis,ref_aln_reads, spikein_aln_reads)



if(mode =="with_input"){

  for (j in 1:length(unique(reads.table$analysis))){
  
  project.name<-unique(reads.table$analysis)[j]
  project.table<- subset(reads.table, reads.table$analysis==project.name)
  
  #ref_aln_reads<-c()
  #spikein_aln_reads<-c()
  calib.eq<-c()
  
  perc_spikein_in_sample<-c()
  perc_spikein_in_input<-c()
  
  project<-c()
  sampname<-c()
  
  scaling_factor<-c()
  
  no_input_scaling<-c()
  
  for (i in 1:length(project.table$id)){
    
    project[i]<-project.name
    sampname[i]<-project.table$id[i]
        
    ref_file=paste0(project.table$id[i],"_ref.flagstat")
    spikein_file=paste0(project.table$id[i],"_spike.flagstat")
    
    #ref_stats<-read.delim(ref_file,sep="\n",header=FALSE)
    #ref_aln_reads[i]<-as.numeric(strsplit(ref_stats[7,1]," \\+ ")[[1]][1])
    
    #spikein_stats<-read.delim(spikein_file,sep="\n",header=FALSE)
    #spikein_aln_reads[i]<-as.numeric(strsplit(spikein_stats[7,1]," \\+ ")[[1]][1])
    
    condition.table<-subset(project.table,project.table$details==project.table$details[i])  
    condition.input<-condition.table[which(condition.table$condition=="INPUT"),]
    
    input_ref_file<-paste0(condition.input$id[1],"_ref.flagstat")
    input_spikein_file<-paste0(condition.input$id[1],"_spike.flagstat")
    
    input_ref<-read.delim(input_ref_file,sep="\n",header=FALSE)
    input_ref_aln<-as.numeric(strsplit(input_ref[7,1]," \\+ ")[[1]][1])
    
    perc_spikein_in_sample[i]=project.table$spikein_aln_reads[i]/(project.table$ref_aln_reads[i]+spikein_aln_reads[i])
    
    if(nrow(condition.input)!=0){ #if there are no inputs skip input based normalization
      input_spikein<-read.delim(input_spikein_file,sep="\n",header=FALSE)
      input_spikein_aln<-as.numeric(strsplit(input_spikein[7,1]," \\+ ")[[1]][1])
      
      #Fursova 2019
      calib.eq[i]<-(1/project.table$spikein_aln_reads[i])*(input_spikein_aln/input_ref_aln)
      
      #Larson 2019
      perc_spikein_in_input[i]=input_spikein_aln/(input_ref_aln+input_spikein_aln)
      scaling_factor[i]=perc_spikein_in_input[i]/perc_spikein_in_sample[i]
    }
    
    
  }
  
  res.table<-data.frame(id=project.table$id, single_end=project.table$single_end, condition=project.table$condition, details=project.table$details,analysis=project.table$analysis,ref_aln_reads, spikein_aln_reads,calib.eq, perc_spikein_in_sample, perc_spikein_in_input, scaling_factor)
  
    if(j==1){
      out.table<-res.table
    } else {
      out.table<-rbind(out.table,res.table)
    
  }
  
  
  
  
}

    out.table.noin<-subset(out.table,out.table$condition!="INPUT")
    alpha=1/max(out.table.noin$calib.eq)
    out.table.noin<-cbind(out.table.noin,alpha=rep(alpha,length(out.table.noin$calib.eq)))
    out.table.noin<-cbind(out.table.noin,downfactor=out.table.noin$calib.eq*out.table.noin$alpha)
    out.table.export<-out.table.noin[,c("id", "single_end", "condition", "details", "analysis", "downfactor")]
    save(out.table.export, file=paste0("norm_output.RData"))
    write.table(out.table.export, file=paste0("allsamples_calib.txt"),
                sep=",", col.names = TRUE, quote=FALSE, row.names = FALSE)

}

#no input normalization
#constant=max(out.table$spikein_aln_reads)
#constant=1000000

#constant/out.table$spikein_aln_reads

if(mode == "without_input"){
  
  out.table<-cbind(reads.table,noinputNorm=scale_fact*(calib_c/reads.table$spikein_aln_reads))

  out.table.noinNorm<-out.table[,c("id", "single_end", "condition", "details", "analysis", "noinputNorm")]

  write.table(out.table.noinNorm, file=paste0("allsamples_calib.txt"),
              sep=",", col.names = TRUE, quote=FALSE, row.names = FALSE)

}







