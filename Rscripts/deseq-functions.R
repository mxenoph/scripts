library("DESeq2")
#library("DESeq")
library("plyr")
library("ggplot2")
library("gplots")
library("gridBase")
library("gridExtra")

# purple, lime, dark red, pinkish, brown, mustard, hotpink, blue, green, grey,
default_colors = c("#B25FCD", "#87CE4F", "#C05038",
                   "#C6AEA9", "#494235", "#B99F47",
                   "#A94777", "#6F7FB5", "#76C2A0",
                   "#A2AFB9")

# More elegant box plot.- should be moved to plotting functions and source that file
boxplot = function (..., col) {
    # Plot very narrow, minimal boxes.
    pars = list(boxwex = 0.4, staplewex = 0,
                medlwd = 1, whisklty = 1,
                frame.plot = F, outpch = 20, outcex = 0.3)
    if (! missing(col))
        pars = c(pars, list(boxcol = col, whiskcol = col, outcol = col))
    graphics::boxplot(..., pars = pars)
}

# Color scheme for divergent colors.
divergent_colors = colorRampPalette(c('#603D71', 'white', '#A4B962'))(30)

# Less intrusive defaults for heatmap.2.
heatmap.2 = function (...)
        gplots::heatmap.2(..., trace = 'none', density.info = 'none',
                          col = divergent_colors, margins = c(8,8))

### FUNCTIONS ###
# create a master file with all counts from the counts directory# {{{
export_counts = function(counts_path) {
    if (file.exists(file.path(counts_path, 'experiment-counts.tsv'))) {
        #check.names=FALSE ensures that no X is put in front of colnames starting with a number
        counts = read.table(file.path(counts_path, 'experiment-counts.tsv'), row.names = 1, header = T, sep = "\t", check.names = FALSE)
    } else {
        counts_files = list.files(counts_path, pattern = "\\.counts$", full.name = T)
        counts_per_sample <- mclapply(counts_files, read.delim, row.names = 1, header = F)
        names(counts_per_sample) = basename(file_path_sans_ext(counts_files))
        
        # Ensure that all files have the same identifiers before merging
        check = unlist(lapply(2:length(counts_per_sample)-1, function(i){
                              s=i+1
                              unlist(lapply(s:length(counts_per_sample), function(j) {
                                            check = identical(rownames(counts_per_sample[[i]]), rownames(counts_per_sample[[j]]))
                                            names(check) = paste0(names(counts_per_sample[i]), ':', names(counts_per_sample[j]))
                                            return(check) })) }))
        if (any(check == FALSE)) {print(check); stop('Count files do not contain the same order and number feature IDs. Please inspect the files and rerun.')}
        
        counts = do.call(cbind.data.frame, counts_per_sample)
        names(counts) = basename(file_path_sans_ext(counts_files))
        write.table(counts, file.path(counts_path, 'experiment-counts.tsv'), row.names = T, col.names = T, quote = F, sep = "\t")
    }
    return(counts)
}# }}}

#Make a plotting function# {{{
plotDE<-function(res, fdr){
 df<-data.frame(res, FDR=ifelse(res$padj<fdr, paste("<", fdr, sep=" "), paste(">", fdr, sep=" ")))

 p<-qplot(log(baseMean), log2FoldChange, data=df, aes(x="baseMean", y="log2FoldChange"))
 p<-p + geom_point(aes(colour=FDR))
 print (p)
}
# }}}

#Make a plotting function for introns--take into account exon de analysis# {{{
plotDEnIntron<-function(res, fdr){
 df<-data.frame(res, FDR=ifelse(res$padj<fdr, paste("<", fdr, sep=" "), paste(">", fdr, sep=" ")))

 p<-ggplot(df, aes(x=log(baseMean), y=log2FoldChange, colour=FDR, shape=DEexon, size=2))
 p<-p + geom_point(aes(colour=FDR, size=FDR, shape=DEexon), size=2)
 print (p)
}
# }}}

#Run DESeq and compare 2 conditions# {{{
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
}# }}}



# Plotting function specific to plotting counts# {{{
# contrast needs to be a character vector of the contrasts used in cds
plot_counts = function(counts, contrast, color = default_colors){
    if(missing(color)){
        color = default_colors
        # If color not passed in the call will use the default colors
        unique_contrast = unique(contrast)
        color = color[1:length(unique_contrast)]
        names(color) = unique_contrast
    }

    layout(matrix(c(1, 2), nrow = 1), widths = c(0.5, 0.5))
    op = par(mar = c(5.1, 5.1, 4.1, 1), oma = c(0, 0, 0, 0))
    
    # Use global minimum value > 0 as pseudocount.
    eps = min(counts[counts > 0])
    boxplot(counts + eps, las = 2, frame.plot = F, log = 'y', col = color[contrast])

    correlated = cor(counts, method = 'spearman')
    pcs = prcomp(correlated)
    explained_variance = summary(pcs)$importance['Proportion of Variance', ]
    plot(PC2 ~ PC1, pcs$x,
         col = color[contrast], pch = 16,
         xlab = sprintf('PC1 (%.0f%% variance)', explained_variance[1] * 100),
         ylab = sprintf('PC2 (%.0f%% variance)', explained_variance[2] * 100))
    legend('topright', bty = 'n', legend = names(color[unique_contrast]), fill = color)
    text(pcs$x[,1], pcs$x[,2], labels = gsub(".*\\.", '', names(pcs$x[,1])), cex = 0.7, pos = 3)

    par(op)

    margin = c(5,5)
    heatmap.2(correlated, ColSideColors = color[contrast])
    
    # Reset par to defaults
    par(mfrow=c(1,1))
    par(mar = c(5, 4, 4, 2) + 0.1)
}
# }}}

