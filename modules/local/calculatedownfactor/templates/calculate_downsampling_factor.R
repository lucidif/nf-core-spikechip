#!/usr/bin/env Rscript

# usage:
# Rscript /path/to/make_twosample_mle_spp.R sample_bam control_bam chromsizes_file mle_output_file"

## LIBRARIES

################################################
# FUNCTIONS
################################################

# sortbychr <- function(x, chrcol="chr", stcol="start", endcol=NULL, chrorder=paste("chr", c(seq(22), "X", "Y"), sep="")) {
#         if (!(is.null(endcol))) {
#                 x <- x[order(x[,endcol]),]
#         }
#         x <- x[order(x[,stcol]),]
#         chrs <- ordered(x[,chrcol], levels=chrorder)
#         # chrs <- ordered(tmp, levels=paste("chr", c(seq(22), "X", "Y"), sep=""))
#         x <- x[order(chrs),]
# }

################################################
## PARAMS
################################################

in_meta <- "${meta}"

g1<-gsub("\\\\[\\\\[",", \\\\[",in_meta)
g2<-gsub(", \\\\[","?",g1)
g3<-gsub("\\\\],","?",g2)
g4<-strsplit(g3, "?")[[1]]
meta_id<-grep("id:",g4)
path_id<-meta_id+1
meta_unformat<-g4[meta_id]

#estrei in nomi delle colonne
univSep<-gsub(":","?",meta_unformat[1])
univSep<-gsub(", ","?",univSep)
col_noformat<-strsplit(univSep,"?")[[1]]
col_names<-col_noformat[seq(1, length(col_noformat), by = 2)]

meta.format<-matrix(nrow=length(meta_unformat), ncol=length(col_names))
colnames(meta.format)<-col_names
#formattiamo i meta
for (i in 1:length(meta_unformat)){
  meta_split<-strsplit(meta_unformat[i],", ")[[1]]
  for(j in 1:length(col_names)){
    torem<-paste0(col_names[j],':')
    meta.format[i,j]<-gsub(torem,'?',meta_split[j])
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

inmatrix<-cbind(meta.format, path.format)

save(inmatrix, file="inputs.RData")


################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

# sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
# print(sessionInfo())
# sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

# r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
# deseq2.version <- as.character(packageVersion('DESeq2'))

# writeLines(
#     c(
#         '"${task.process}":',
#         paste('    r-base:', r.version),
#         paste('    bioconductor-deseq2:', deseq2.version)
#     ),
# 'versions.yml')

################################################
################################################
################################################
################################################
