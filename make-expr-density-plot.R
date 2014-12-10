#!/usr/bin/env Rscript

source("~/source/Rscripts/functions.R")# {{{
source("~/source/Rscripts/granges-functions.R")
x <- c('plyr', 'ggplot2', 'gridExtra', 'reshape')
lapply(x, suppressMessages(library), character.only=T)# }}}

load("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/mm9.e67.ann.Rdata")
gene_assoc_IDs <- read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/MM9.maps/mm9.e67.associated.geneIDs.txt", header=TRUE, as.is=TRUE, sep="\t")

args <- commandArgs(trailingOnly=TRUE)
peaks <- macs2GRanges(args[1])
expr <- deseq2vect(args[2])
# User defined subset list i.e. ChEA/M2-Chd4.Epi.SUZ12_CommonTargets.txt
subsUL <- read.table(args[3], as.is=T, sep="\t")[[1]]
out <- args[4]
desc <- args[5]

# Check if 'plots' dir exists in the wdr and create if not
plots <- file.path(out, 'plots', '')
dir.create(plots)

genes <- my.annotation$genes
extended_genes <- resize(promoters(genes, upstream=2000), width=width(genes)+2500)

ov <- findOverlaps(extended_genes, peaks)
subs <- names(extended_genes[unique(queryHits(ov))])

if(length(unique(grepl('ENS',subsUL))) > 1){
    stop("List provided for subsetting does not include a single type of identifiers.")
} else if(unique(grepl('ENS',subsUL)) == 'FALSE'){
    subsUL <- gene_assoc_IDs[toupper(gene_assoc_IDs$Associated.Gene.Name) %in% subsUL, 'Ensembl.Gene.ID']
}

maxXLim <- round(max(abs(min(expr$df$log2FoldChange)), abs(max(expr$df$log2FoldChange))) + 0.5)

pdf(paste0(plots, desc, 'ExpressionDensity.pdf'), paper='a4')
x <- ggplot(data=expr$df, aes(x=log2FoldChange))
x <- x + geom_density(aes(color=padj < 0.05), alpha=.3)
#x <- x + coord_cartesian(xlim= c(-maxXLim, maxXLim))
# Manually setting the xlim to -4 and +4
x <- x + coord_cartesian(xlim= c(-4.5, 4.5))
x <- x + theme(plot.title=element_text(hjust= 150),
               legend.position= "top")

# Get plots legend
legend <- g_legend(x)
x <- x + theme(legend.position= "none")

grid.arrange(legend, textGrob(' '),
             ggplotGrob(x), textGrob("ALL", rot=270),
             ggplotGrob(x %+% expr$df[rownames(expr$df) %in% subs,]), textGrob('Bound', rot=270),
             ggplotGrob(x %+% expr$df[rownames(expr$df) %in% subsUL,]), textGrob('Common Targets', rot=270), nrow=4, ncol=2, widths=unit(c(12,2), "cm"),
             # Removes the new page printed by ggplot before the plots
             newpage=FALSE)
dev.off()

