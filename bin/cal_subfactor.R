#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

in_meta <- args[1]

g1<-gsub("\\[\\[",", \\[",in_meta)
g2<-gsub(", \\[","?",g1)
g3<-gsub("\\],","?",g2)
g4<-strsplit(g3, "\\?")[[1]]
meta_id<-grep("id\\:",g4)
path_id<-meta_id+1
meta_unformat<-g4[meta_id]

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
path.unformat<-g4[path_id]
path.unformat<-gsub(" ","",path.unformat)
path.format=matrix(ncol=2,nrow=length(path.unformat))
colnames(path.format)<-c("ref_file","spikein_file")

for(i in 1:length(path.unformat)){
  
  paths<-strsplit(path.unformat[i],",")[[1]]
  path.format[i,1]<-paths[1]
  path.format[i,2]<-paths[2]
  
}

path.format<-gsub("\\]","",path.format)
#path.format[,2]<-gsub("\\]","",path.format[,2])

inmatrix<-cbind(meta.format, path.format)

save(inmatrix, file="inputs.RData")

input.table<-as.data.frame(inmatrix)
ref_aln_reads<-c()
spikein_aln_reads<-c()
calib.eq<-c()
for (i in 1:length(input.table$id)){
  ref_stats<-read.delim(input.table$ref_file[i],sep="\n",header=FALSE)
  ref_aln_reads[i]<-as.numeric(strsplit(ref_stats[7,1]," \\+ ")[[1]][1])
  
  spikein_stats<-read.delim(input.table$spikein_file[i],sep="\n",header=FALSE)
  spikein_aln_reads[i]<-as.numeric(strsplit(spikein_stats[7,1]," \\+ ")[[1]][1])

  condition.table<-subset(input.table,input.table$condition==input.table$condition[i])  
  condition.input<-condition.table[which(condition.table$details=="INPUT"),]
  
  input_ref<-read.delim(condition.input$ref_file[1],sep="\n",header=FALSE)
  input_ref_aln<-as.numeric(strsplit(input_ref[7,1]," \\+ ")[[1]][1])
  
  input_spikein<-read.delim(condition.input$spikein_file[1],sep="\n",header=FALSE)
  input_spikein_aln<-as.numeric(strsplit(input_spikein[7,1]," \\+ ")[[1]][1])
  
  calib.eq[i]<-(1/spikein_aln_reads[i])*(input_spikein_aln/input_ref_aln)
    
}

out.table<-data.frame(id=input.table$id, single_end=input.table$single_end, condition=input.table$condition, details=input.table$details,ref_aln_reads, spikein_aln_reads,calib.eq)

out.table.noin<-subset(out.table,out.table$details!="INPUT")

alpha_ref_calib<-out.table.noin[which(out.table.noin$calib.eq==max(out.table.noin$calib.eq)),"calib.eq"]

alpha=1/max(out.table.noin$calib.eq)

out.table.downfact<-data.frame(out.table.noin, down.factor=out.table.noin$calib.eq*alpha)

out.table.export<-merge(input.table,out.table.downfact[,c("id","down.factor")],by=1)

#total_ref_aln_reads<-sum(ref_aln_reads)
#total_spikein_aln_reads<-sum(spikein_aln_reads)





