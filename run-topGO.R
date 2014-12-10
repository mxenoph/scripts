#!/usr/bin/env Rscript

x <- c('topGO', 'ggplot2', 'org.Mm.eg.db')# {{{
lapply(x, suppressMessages(library), character.only=T)
source("~/source/Rscripts/functions.R")# }}}

gene_assoc_IDs <- read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/MM9.maps/mm9.e67.associated.geneIDs.txt", header=TRUE, as.is=TRUE, sep="\t")

args <- commandArgs(trailingOnly=TRUE)
f <- read.table(args[1], sep="\t")
descr <- args[2]
out_dir <- args[3]

binary <- gene_assoc_IDs[['Ensembl.Gene.ID']]
# Check if 'plots' dir exists in the wdr and create if not
plot_dir <- file.path(out_dir, 'plots', '')
dir.create(plot_dir)
setwd(plot_dir)

#de <- deseq2vect(f)
eIDs <- gene_assoc_IDs[toupper(gene_assoc_IDs[['Associated.Gene.Name']]) %in% f[,1, drop=T], 'Ensembl.Gene.ID']
binary <- factor(as.integer(gene_assoc_IDs[['Ensembl.Gene.ID']] %in% eIDs)); names(binary) <- gene_assoc_IDs[['Ensembl.Gene.ID']]
print(head(binary))

#GOtermEnr(de$binary, descr, "BP", "org.Mm.eg.db")
t <- GOtermEnr(binary, descr, "BP", "org.Mm.eg.db")
pdf(paste0(info,"GO.pdf"), paper='a4')
ggplot(t$table,aes(y=Significant/Annotated,x=Term,fill=weight01Fisher)) + geom_bar(stat="identity", aes(order=desc(weight01Fisher))) + scale_fill_gradient(low="red",high="yellow")   + coord_flip()
dev.off()


