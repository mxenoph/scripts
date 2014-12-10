library(GenomicRanges)
library(GenomicFeatures)
library(biovizBase)
library(biomaRt)
library(ggplot2)
library(ggbio)
library(rtracklayer)


txNexons<-function(args){
#passing the host to be used as an argument e.g.:args<-"jan2013.archive.ensembl.org"
txdb<-makeTranscriptDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filters="", id_prefix="ensembl_", host=args, miRBaseBuild=NA)

allTx <- transcriptsBy(txdb, by='gene')
}