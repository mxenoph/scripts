source("~/source/Rscripts/deseq-functions.R")

#Run on the hendrich data; argument should be the directory of the files: REMEMBER to correct it if you set more args
my_matrix<-make_matrix("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/")

#Set a vector for the conditions; REMEMBER to check if your files are retrieved in the same order-SAME ORDER
conds <- c("WT_2i","WT_Lif","WT_N2B27","WT_NoLif","WT_2i","WT_Lif","WT_N2B27","WT_NoLif","KO_2i","KO_Lif","KO_N2B27","KO_NoLif","KO_2i","KO_Lif","KO_N2B27","KO_NoLif")

cds<-newCountDataSet( my_matrix, conds )

#Why not using newCountDataSet-FromHTSeqCount function?
#design<-data.frame(sample.names= , count.files= , condition=conds)
#cds<-newCountDataSet-FromHTSeqCount(design, directory="")

#Estimate size factors
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cds <- estimateDispersions(cds)
#cds <- estimateDispersions(cds,method="per-condition")

#experimental conditions
conditions<-c("2i","Lif","N2B27","NoLif")


#REMEMBER need to fill in arguments for when calling de_pair (i.e. cds (1st arg) and file name where gene Id->gene name) -DONE
for (c in c(1:4)) {
 #Call function de_pair to compare to conditions WT Vs KO under same comdition
 de_pair(cds, paste("WT_",conditions[[c]],sep=""), paste("KO_",conditions[[c]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
 #Second iteration for comparing same cells under treatment and no treatment conditions
 #why is 2i Vs N2B27 or 2i Vs NoLif informative though?
 for (cprime in c(2:4)) {
  if (c < cprime) {
 #Second iteration for comparing same cells under treatment and no treatment conditions
   de_pair(cds, paste("WT_",conditions[[c]],sep=""),paste("WT_",conditions[[cprime]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
 #Second iteration for comparing same cells under treatment and no treatment conditions
   de_pair(cds, paste("KO_",conditions[[c]],sep=""),paste("KO_",conditions[[cprime]],sep=""), "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
  }
 }
}


# General comparison WT - KO

conds <- c("WT","WT","WT","WT","WT","WT","WT","WT","KO","KO","KO","KO","KO","KO","KO","KO")
cds<-newCountDataSet( my_matrix, conds )
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cds <- estimateDispersions(cds)
#REMEMBER: fill in the missing args
de_pair(cds, "WT","KO", "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol")
