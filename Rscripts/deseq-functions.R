library("DESeq")
library("plyr")
library("ggplot2")

### FUNCTIONS ###

#Get htseq-count files and create a data matrix
make_matrix<-function(dir){
 
 #Get the files-REMEMBER to fill in the paste arguments with the <samout> option set for htseq once you run it-DONE
 files<-list.files(path=dir, pattern= glob2rx(paste("*",".counts", sep="")), full.names=TRUE)
 #Save the names (cells and conditions) for the files-REMEMBER to fill it in-DONE
 #sub("pattern", "replacement", x)
 names<-gsub("\\..*", "", files)
 names<-gsub(".*/", "", names)

 #Read the files and create a list of tables
 counts_list<-lapply(files, read.table, row.names=1)
 #Initialize a matrix
 counts_matrix<- NULL
 #Iterate through the elements of the list (i.e. the tables)

 for(i in 1:length(counts_list)){
 #for(i in 1:length(counts_list)){
  #Stop and produce an error message indicating the first of the elements of which were not true if not all tables have the same genes (rownames) 
  stopifnot(all(rownames(counts_list[[i]])== rownames(counts_list[[1]])))
  #make a matrix. Our table only has two columns, with the first one (V1) being NULL because contains the genes and the V2 containing the counts for each gene
  counts_matrix<-cbind(counts_matrix, counts_list[[i]]$V2)
 }
 #Get the gene names
 rownames(counts_matrix)<-rownames(counts_list[[1]])
 colnames(counts_matrix)<-names
 
 #Removing the last 5 lines of the HTSeq-Count output
 counts_matrix<-subset(counts_matrix, grepl("ENSMUSG", rownames(counts_matrix)))
 return(counts_matrix)
}

#Make a plotting function
plotDE<-function(res, fdr){
 df<-data.frame(res, FDR=ifelse(res$padj<fdr, paste("<", fdr, sep=" "), paste(">", fdr, sep=" ")))

 p<-qplot(log(baseMean), log2FoldChange, data=df, aes(x="baseMean", y="log2FoldChange"))
 p<-p + geom_point(aes(colour=FDR))
 print (p)
}

#Make a plotting function for introns--take into account exon de analysis
plotDEnIntron<-function(res, fdr){
 df<-data.frame(res, FDR=ifelse(res$padj<fdr, paste("<", fdr, sep=" "), paste(">", fdr, sep=" ")))

 p<-ggplot(df, aes(x=log(baseMean), y=log2FoldChange, colour=FDR, shape=DEexon, size=2))
 p<-p + geom_point(aes(colour=FDR, size=FDR, shape=DEexon), size=2)
 print (p)
}

#Run DESeq and compare 2 conditions
NbinomTest <- function(cds, cond1, cond2, geneNames, exonDE=NULL) {
 #Run the negative binomial test on the CountDataSet(cds with estimatedDispersions) for the 2 conditions
 res <- nbinomTest( cds, cond1, cond2 )

 if(is.null(exonDE)==FALSE) res<-data.frame(res, DEexon=ifelse(res$id %in% exonDE$id, "yes", "no"))

 #will need to give it a more descriptive name and check if the conditions are printed in the plot
 pdf(paste(cond1, "vs", cond2, ".pdf", sep=""))
 ifelse(is.null(exonDE)==FALSE, plotDEnIntron(res, .05), plotDE(res, .05))
 dev.off()

 #Filter for significant genes for FDR threshold of 5% (also don't take into account genes for which padj was not defined (i.e. NA) because fold change is NaN (i.e. division of 0 by 0 (counts/size/factor-??)))
 resSig <- res[!is.na(res$padj) & res$padj < .05, ]
 #Order the genes based on p-value
 resSig<<-resSig[ order(resSig$pval), ]

 #Get the Gene Name for the e! Gene Ids of the de significant genes-REMEMBER: cannot give remco's file cause he got the genes names probably based on gene ids for MM9
 #mouse_data<-read.delim(id_to_name_file, as.is=TRUE, quote="")
 geneNames<-read.delim(geneNames, as.is=TRUE, quote="")

 #Get the strongly downregulated of the significant genes--this will just order it so the downreg come first--the list still contains the upreg
 downRegSig<-resSig[order(resSig$foldChange, -resSig$baseMean), ] 
 down_mouse_names<-geneNames$mgi_symbol[ match(downRegSig$id,geneNames$ensembl_gene_id) ]
 down_genes<-cbind(downRegSig,Name=down_mouse_names)

 write.table(down_genes, file=paste(cond1,"vs",cond2,"-downRegSig.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

 #Get the strongly upregulated of the significant genes
 upRegSig<-resSig[order(-resSig$foldChange, -resSig$baseMean), ]
 up_mouse_names<-geneNames$mgi_symbol[ match(upRegSig$id,geneNames$ensembl_gene_id) ]
 up_genes<-cbind(upRegSig,Name=up_mouse_names)

 write.table(up_genes, file=paste(cond1,"vs",cond2,"-upRegSig.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

 #Change the match arguments if you have different header--what happens for ids that don't have an mgi symbol?
 de_mouse_names<-geneNames$mgi_symbol[ match(resSig$id,geneNames$ensembl_gene_id) ]
 #DE significant genes (id to name)
 degenes<-cbind(resSig,Name=de_mouse_names)
 #Save an object of the de genes
 save(degenes, file=paste(cond1,"vs",cond2,"-de.Rda",sep=""))  

 #Write output-all de significant genes with gene name instead of gene id
 write.table(degenes, file=paste(cond1,"vs",cond2,"-de.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

 #Write output for all genes ordered by p-val
 resAll<-res[ order(res$pval), ]

 all_mouse_names<-geneNames$mgi_symbol[ match(resAll$id,geneNames$ensembl_gene_id) ]
 allgenes<-cbind(resAll,Name=all_mouse_names)
 #Save an object for all genes 
 save(allgenes, file=paste(cond1,"vs",cond2,"-de.Rda",sep=""))
 
 #Write output for all genes with gene name instead of id
 write.table(allgenes, file=paste(cond1,"vs",cond2,".txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
}
