source("~/source/Rscripts/deseq-functions.R")

#Function to create the cds--calls RunNbinom for further analysis--flag should be "uniform" if design only has samples sequenced with one library preparation type
MakeCDSnRunTest<-function(design, flag, geneNameFile ){
 GeneNames<-read.delim(geneNameFile, as.is=TRUE, quote="")
 
 cds<-newCountDataSetFromHTSeqCount(design, directory="")

 #Estimate size factors
 cds <- estimateSizeFactors( cds )
 #sizeFactors( cds )
 cds <- estimateDispersions(cds)

 if(flag=="uniform") {
  #Experimental conditions
  conditions<-unique(gsub(".*_", "", design[[3]]))

  RunNbinom(cds, conditions, GeneNames)
 }

 if(flag=="glm"){
  conditions<-unique(design[[3]])
  ret<-RunGLM(cds, conditions, GeneNames)
 }
 return(list(cds=cds, ret=ret))
}

RunNbinom<-function(cds, conditions, geneNames){
 #REMEMBER need to fill in arguments for when calling NbinomTest (i.e. cds (1st arg) and file name where gene Id->gene name) -DONE
 for (c in c(1:length(conditions))) {
  #Call function NbinomTest to compare to conditions WT Vs KO under same comdition
  NbinomTest(cds, paste("WT_",conditions[[c]],sep=""), paste("KO_",conditions[[c]],sep=""), geneNames)

  #If only one condition in the design there is no need for the second loop
  if (c==length(conditions)) break()

  #Second iteration for comparing same cells under treatment and no treatment conditions
  for (cprime in c(2:length(conditions))) {
   if (c < cprime) {
    #Second iteration for comparing same cells under treatment and no treatment conditions
    NbinomTest(cds, paste("WT_",conditions[[c]],sep=""),paste("WT_",conditions[[cprime]],sep=""), geneNames)
    #Second iteration for comparing same cells under treatment and no treatment conditions
    NbinomTest(cds, paste("KO_",conditions[[c]],sep=""),paste("KO_",conditions[[cprime]],sep=""), geneNames)
   }
  }
 }
}


RunGLM<-function(cds, conditions, geneNames){
 fit1<-fitNbinomGLMs(cds, count~libType + condition)
 fit0<-fitNbinomGLMs(cds, count~libType)
 
 pvalsGLM<-nbinomGLMTest(fit1, fit0)
 padjGLM<-p.adjust(pvalsGLM, method="BH")

 #DE genes at FDR<0.05
 de<-fit1[which(!is.na(padjGLM) & padjGLM<0.05),]
 matchGeneNames<-geneNames$mgi_symbol[match(rownames(de), geneNames$ensembl_gene_id)]
 deWithGName<-cbind(de, Names=matchGeneNames)

 #If statement to get the name for the output file according to the comparison made by the GLM
 if(conditions[1]==gsub("condition", "", colnames(fit1)[3])){ conditions<-rev(conditions)}
 write.table(deWithGName, file=paste(conditions[1], "vs", conditions[2], "-de.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE) 

 return(list(fit1=fit1, padjGLM=padjGLM, de=de))
}


###End-of-Functions###


geneNameFile<-"/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol"

design<-read.table("/nfs/research2/bertone/user/mxenoph/hendrich/hendrich_design.txt", header=TRUE, sep='\t', as.is=TRUE, quote="")

#For the SE data
designSE<-design[which(design$libType=="single-end"), !colnames(design) %in% "libType"]
designSE<-designSE[order(designSE$condition, decreasing=TRUE),]
#return<-MakeCDSnRunTest(designSE, "uniform", geneNameFile)

#For the PE data
designPE<-design[which(design$libType=="paired-end"), !colnames(design) %in% "libType"]
designPE<-designPE[order(designPE$condition, decreasing=TRUE),]
#return<-MakeCDSnRunTest(designPE, "uniform", geneNameFile)

for(pe in unique(gsub(".*_", "", designPE[[3]]))){
 for(se in unique(gsub(".*_", "", designSE[[3]]))){
  designMix<-subset(design, condition %in% c(paste("WT_", pe, sep=""), paste("KO_", se, sep="")))
  designMix<-designMix[order(designMix$condition, decreasing=TRUE),]
  myreturn<-MakeCDSnRunTest(designMix, "glm", geneNameFile)

  designMix<-subset(design, condition %in% c(paste("WT_", pe, sep=""), paste("WT_", se, sep="")))
  myreturn1<-MakeCDSnRunTest(designMix, "glm", geneNameFile)

  designMix<-subset(design, condition %in% c(paste("KO_", pe, sep=""), paste("KO_", se, sep=""))) 
  myreturn2<-MakeCDSnRunTest(designMix, "glm", geneNameFile)
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
