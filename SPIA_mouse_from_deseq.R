options(stringsAsFactors=F)
library(SPIA)
#library(pathview)

#files <- commandArgs(trailingOnly = TRUE)
files <- "/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi.txt"
#files <- "/nfs/nobackup2/research/bertone/mxenoph/bh/forMM9/deseq/WT_EpivsKO_Epi.txt"

mouse.data<-read.delim("/nfs/research2/bertone/user/remco/data/ens70_mouseids.txt",as.is=TRUE,quote="")
#mouse.data<-read.delim("/nfs/research2/bertone/user/remco/data/mouse_ids.txt",as.is=TRUE,quote="")

entrez.data <- mouse.data[!is.na(mouse.data$EntrezGene.ID),]

data.list = lapply(files, read.delim,as.is=TRUE,quote="")
names(data.list) <- sub(".txt","",files)

run.spia <- function(name,deseq.results) {

  de.data <- deseq.results[deseq.results$id %in% entrez.data$Ensembl.Gene.ID,]

  de.data$log2FoldChange[de.data$log2FoldChange > 10] <-  10
  de.data$log2FoldChange[de.data$log2FoldChange < -10] <- -10
  
  all <-  entrez.data$EntrezGene.ID[match(de.data$id,entrez.data$Ensembl.Gene.ID)]

  de <- de.data$log2FoldChange[!is.na(de.data$padj) & de.data$padj < 0.20]
  names(de) <- as.vector(entrez.data$EntrezGene.ID[match(de.data$id[!is.na(de.data$padj) & de.data$padj < 0.20],entrez.data$Ensembl.Gene.ID)])
  de <- de[!duplicated(names(de))]
  
  res=spia(de=de,all=all,organism="mmu",nB=1000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)

  sign.res <- res[res$pGFdr < 0.1,]

  write.table(sign.res, file=paste(name,"_enriched_pathways0.2_loose.txt", sep=""), quote=FALSE, row.names=FALSE, sep="\t")

}

for (i in 1:length(data.list)) {
  run.spia(names(data.list)[i],data.list[[i]])
}


