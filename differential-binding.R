#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")

parser <-  ArgumentParser(description="Perform differential binding analysis")
parser$add_argument('-s', '--sheet', metavar= "file", required='True', type= "character", help= "Sample sheet must have SampleID,Condition,Treatment,Replicate,bamReads,bamControl,Peaks (last 3 are the path to the respective bed files)")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args <- parser$parse_args()

parent_path <- file.path(args$out, 'differential_binding')
plot_path <- file.path(parent_path, 'plots')
#}}}

#setwd("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/macs/")

#SampleSheet must have SampleID,Condition,Treatment,Replicate,bamReads,bamControl,Peaks (last 3 are the path to the respective bed files)
DB<-function(sheet, fdr=0.05){
   library(DiffBind)

   output <- gsub(".csv", 'pdf', gsub(".*/", '', sheet))

   #Cause it will otherwise complain and terminate when you run it on a cluster
   pdf(output)
   #working object is of dba class

   #Check how many Conditions are there
   line <- unlist(strsplit(tolower(readLines(sheet[1], n=1)), ','))
   if(!"condition" %in% line) stop("Sample Sheet does not contain a Condition column")

   protein <- dba(sampleSheet= sheet)
   #Will print the correlation heatmap using occupancy(peak caller score) data
   plot(protein)

   sheet_table <-read.table(sheet, sep=',')
   conditions <- sort(levels(sheet_table$condition), decreasing=T)

   nCond <- match("condition", line)
   tmp <- rep('NULL', length(line))
   tmp[nCond] <- NA
   nCond <- sort(levels(read.table(sheet, colClasses=tmp, sep=",", head=T)[[1]]), decreasing=T)

   #calculates a binding matrix, scores based on read counts (affinity scores) rather than the previous one that plotted confidence scores.
   #DBA_SCORE_READS_MINUS => read count for interval-read count for interval in control..might be a better way of doing this
   #bScaleControl=TRUE => scale the control reads based on relative library size
   wobj <- dba.count(wobj, score=DBA_SCORE_READS_MINUS, bScaleControl=TRUE, bParallel=TRUE)
   btag <- TRUE

   if(sum(wobj$masks[[nCond[1]]]) < 2 | sum(wobj$masks[[nCond[2]]]) < 2){
      wobj <- dba.contrast(wobj, wobj$masks[[nCond[1]]], wobj$masks[[nCond[2]]], nCond[1], nCond[2])
      btag <- FALSE
   }
   else{
      #step is actually optional cause it's set up
      wobj <- dba.contrast(wobj, categories=DBA_CONDITION)
   }


   wobj <- dba.analyze(
                       wobj,
                       method= DBA_DESEQ,
                       bSubControl= TRUE,
                       #bFullLibrarySize= TRUE,
                       #SOS if no replicates this option should be set to FALSE
                       bTagwise= btag,
                       bParallel= TRUE)

   wobj.DB <- dba.report(wobj,
                         method=DBA_DESEQ,
                         #only sites with an absolute Fold value greater than equal to this will be included in the report
                         #fold= 2,
                         th= fdr,
                         initString= fdr,
                         file= paste(gsub(".csv", '', gsub(".*/", '', sheet)), '.diffbind', sep=''))
   wobj.df <- dba.report(wobj, method=DBA_DESEQ, th=fdr, DataType = DBA_DATA_FRAME)


   write.table(wobj.df, file=gsub(".csv", 'pdf', gsub(".*/", '', sheet)), sep="\t", quote=FALSE, row.names=FALSE)
   #plot(wobj)
   #dba.plotMA(wobj)
   #PCA based on the afinity data for all sites
   #dba.plotPCA(wobj, DBA_CONDITION)
   #PCA based on the affiity data for the DB sites
   #dba.plotPCA(wobj, contrast=1, th=.05)
   dev.off()
   return(list("db"=wobj.DB))
}

#myPeakFiles <- list.files(path=getwd(), full.name=T )
#samplesheets must point to modified macs peaks.bed with 4 columns instead of five
mySampleSheets <- list.files("..", full.name=T, pattern="csv")

#system.time(DB(mySampleSheets[1]), .01)


#Comment in and out depending what you are running
#Mi2b_2i <- DB(mySampleSheets[1], 0.05)
#Mi2b_Epi <- DB(mySampleSheets[2], 0.05)
#Mi2b_serum <- DB(mySampleSheets[3], 0.05)



#source("filewith get.annotation")
#mm9.annot<-get.annotation('may2012.archive.ensembl.org', 'mmusculus_gene_ensembl')





