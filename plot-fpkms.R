#!/usr/bin/env Rscript

source("~/source/Rscripts/functions.R")# {{{
source("~/source/Rscripts/granges-functions.R")
x <- c('topGO', 'ggplot2')
lapply(x, suppressMessages(library), character.only=T)# }}}

args <- commandArgs(trailingOnly=TRUE)
fpkms <- args[1]
# Either as file or string separated by ':'
genesIDs2plot <- args[2]


ensembl<- useMart(host='jan2013.archive.ensembl.org', dataset='mmusculus_gene_ensembl', biomart='ENSEMBL_MART_ENSEMBL')

gene_data <- getBM(attributes=c('ensembl_gene_id', 'go_id'), filters='go_id', values=paste0('GO:', c('0003365','0003379','0042074','0006932','0061450','0016477','0048870','0070358','2000145','0051674','0051234','0051641','0007163','0030952','0016332','0040011')), mart=ensembl)
e70AssocName <- getBM(attributes=c('ensembl_gene_id', 'external_gene_id'), mart=ensembl)
gene_data$name <- e70AssocName[match(gene_data$ensembl_gene_id, e70AssocName$ensembl_gene_id), 'external_gene_id']
gene_data$de_epi <- gene_data$ensembl_gene_id %in% de[['WT_EpivsKO_Epi-de.txt']][, 'id']
go.fpkms <- fpkms[gene_data[gene_data$de_epi == TRUE, 'ensembl_gene_id'], grep('Epi', colnames(fpkms))]
plot.data <- melt(go.fpkms, id=rownames(go.fpkms))
plot.data[grep("2lox_Epi|3Flox_Epi", plot.data[['X2']]), 'Group'] <-'WT'
plot.data[grep("spl2_Epi|3KO_Epi", plot.data[['X2']]), 'Group'] <-'KO'
plot.data[, 'Condition'] <- 'EpiSCs'
plot.data$GO.ID <- gene_data[match(plot.data$X1, gene_data$ensembl_gene_id), 'go_id']
levels(plot.data$X1) <- e70AssocName[match(levels(plot.data$X1), e70AssocName$ensembl_gene_id), 'external_gene_id']

pdf("go.maria.pdf")
lapply(unique(plot.data[, 'GO.ID']), function(x){
       df <- plot.data[plot.data$GO.ID== x,]
       p <- ggplot(df, aes(x=Condition, y=value, color=Group)) + geom_point()
       p <- p + geom_errorbar(stat = "hline", yintercept = "mean", width=0.8,aes(ymax=..y..,ymin=..y..))
       p <- p + facet_grid(. ~ X1) + theme_bw() + xlab('') + ylab("FPKM") +theme(panel.grid.major.x = element_blank())
       p <- p + theme(strip.text.x = element_text(size = 10, colour = "black", angle = 90)) + ggtitle(x)
       print(p)
       rm(p)
})
dev.off()

write.table(gene_data, file="go4maria_table.txt", rownames=F, sep="\t")

