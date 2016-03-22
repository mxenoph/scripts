source("~/source/Rscripts/deseq-functions.R")

#Why not using newCountDataSet-FromHTSeqCount function?
#This is for differentially expressed introns
design<-read.table("/nfs/research2/bertone/user/mxenoph/hendrich/hendrich_design_intron.txt", header=TRUE, sep='\t', as.is=TRUE, quote="")
#Directory is "" if in the design file the full path is given
cds<-newCountDataSetFromHTSeqCount(design, directory="")

#Estimate size factors
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cds <- estimateDispersions(cds)

#experimental conditions
conditions<-unique(gsub(".*_", "", design[[3]]))
exonDE<-read.delim("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi-de.txt", as.is=T, quote="")
names<-"/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol"

#Call function de_pair to compare to conditions WT Vs KO under same comdition
de_pair(cds, paste("WT_",conditions[[1]],sep=""), paste("KO_",conditions[[1]],sep=""), names, exonDE)

