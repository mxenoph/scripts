source("~/source/Rscripts/functions.R")# {{{
#for selecting substrings in perl-like way
x <- c('topGO', 'org.Mm.eg.db', 'VennDiagram')
lapply(x, suppressMessages(library), character.only=T)# }}}

path <- "/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal"
wt <- file.path(path, "WT_2ivsWT_Epi.txt")
ko <- file.path(path, "KO_2ivsKO_Epi.txt")

plot_path <- file.path(path, 'plots')
dir.create(plot_path, recursive= TRUE)

wt <- deseq2vect(wt)
ko <- deseq2vect(ko)

sets_union <- table(rownames(wt$de) %in% rownames(ko$de))[['TRUE']]
wt_total <- nrow(wt$de)
ko_total <- nrow(ko$de)

pdf(file.path(plot_path, 'transition.pdf'), paper= 'a4')
draw.pairwise.venn(wt_total, ko_total, sets_union,
                   category= c('wt', 'ko'),
                   scaled = TRUE)
dev.off()

common_genes <- rownames(wt$de)[rownames(wt$de) %in% rownames(ko$de)]
wt_specific <- wt$de[!rownames(wt$de) %in% rownames(ko$de)]
ko_specific <- ko$de[!rownames(ko$de) %in% rownames(wt$de)]

gene_assoc_IDs <- read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/MM9.maps/mm9.e67.associated.geneIDs.txt", header=TRUE, as.is=TRUE, sep="\t")

f <- read.table(args[1], sep="\t")
descr <- 'transition'

binary <- gene_assoc_IDs[['Ensembl.Gene.ID']]
# Check if 'plots' dir exists in the wdr and create if not
plot_dir <- plot_path

#de <- deseq2vect(f)
eIDs <- gene_assoc_IDs[toupper(gene_assoc_IDs[['Associated.Gene.Name']]) %in% , 'Ensembl.Gene.ID']
binary <- factor(as.integer(gene_assoc_IDs[['Ensembl.Gene.ID']] %in% eIDs)); names(binary) <- gene_assoc_IDs[['Ensembl.Gene.ID']]
print(head(binary))

#GOtermEnr(de$binary, descr, "BP", "org.Mm.eg.db")
t <- GOtermEnr(binary, descr, "BP", "org.Mm.eg.db")
pdf(paste0(info,"GO.pdf"), paper='a4')
ggplot(t$table,aes(y=Significant/Annotated,x=Term,fill=weight01Fisher)) + geom_bar(stat="identity", aes(order=desc(weight01Fisher))) + scale_fill_gradient(low="red",high="yellow")   + coord_flip()
dev.off()


