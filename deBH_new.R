# Correct one to use

source("~/source/Rscripts/deseq-functions.R")


#Why not using newCountDataSet-FromHTSeqCount function?
design<-read.table("/nfs/research2/bertone/user/mxenoph/hendrich/hendrich_design_4nbinom.txt", header=TRUE, sep='\t', as.is=TRUE, quote="")

#Directory is "" if in the design file the full path is given
cds<-newCountDataSetFromHTSeqCount(design, directory="")

#Estimate size factors
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
conditions<-design[[3]]
cds <- estimateDispersions(cds)

#experimental conditions
conditions<-unique(gsub(".*_", "", design[[3]]))


#REMEMBER need to fill in arguments for when calling de_pair (i.e. cds (1st arg) and file name where gene Id->gene name) -DONE
for (c in c(1:5)) {
 #Call function de_pair to compare to conditions WT Vs KO under same comdition
 #de_pair(cds, paste("WT_",conditions[[c]],sep=""), paste("KO_",conditions[[c]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
 NbinomTest(cds, paste("WT_",conditions[[c]],sep=""), paste("KO_",conditions[[c]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
 #Second iteration for comparing same cells under treatment and no treatment conditions
 #why is 2i Vs N2B27 or 2i Vs NoLif informative though?
 for (cprime in c(2:5)) {
  if (c < cprime) {
 #Second iteration for comparing same cells under treatment and no treatment conditions
   #de_pair(cds, paste("WT_",conditions[[c]],sep=""),paste("WT_",conditions[[cprime]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
   NbinomTest(cds, paste("WT_",conditions[[c]],sep=""),paste("WT_",conditions[[cprime]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
 #Second iteration for comparing same cells under treatment and no treatment conditions
   #de_pair(cds, paste("KO_",conditions[[c]],sep=""),paste("KO_",conditions[[cprime]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
   NbinomTest(cds, paste("KO_",conditions[[c]],sep=""),paste("KO_",conditions[[cprime]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
  }
 }
}


# General comparison WT - KO

#conds <- c("WT","WT","WT","WT","WT","WT","WT","WT","KO","KO","KO","KO","KO","KO","KO","KO")
#cds<-newCountDataSet( my_matrix, conds )
#cds <- estimateSizeFactors( cds )
#sizeFactors( cds )
#cds <- estimateDispersions(cds)
#REMEMBER: fill in the missing args
#de_pair(cds, "WT","KO", "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
