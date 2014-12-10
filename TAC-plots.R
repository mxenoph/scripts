load("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/TestMe/FPKMs.Rdata")
write.table(fpkms, file="/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/TestMe/fpkms.txt", row.names=T, col.names=T, sep="\t")


condition <- mclapply(c('2i', '_Lif', 'Epi'), function(i)grep(i,colnames(fpkms)))
names(condition) <- c('2i','Lif', 'Epi')
average_fpkms <- mclapply(condition, function(v){ 
                          wt <- grep("2lox|3Flox", colnames(fpkms)); 
                          WT <- rowMeans(fpkms[,v[v %in% wt]])
                          KO <- rowMeans(fpkms[,v[!v %in% wt]])
                          data.frame(WT,KO)
})
# combine to a single dataframe
average_fpkms <- do.call(cbind, average_fpkms)
write.table(average_fpkms, file="/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/TestMe/average_fpkms.txt", row.names=T, col.names=T, sep="\t", quote=F)
deseq_report <- list('2i'= deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_2ivsKO_2i.txt")$raw,
           'Lif' = deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_LifvsKO_Lif.txt")$raw,
           'Epi' = deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi.txt")$raw)

deseq_report_fpkm <- mclapply(c('2i', 'Lif', 'Epi'), function(i){
         v <- grep(i, colnames(average_fpkms))
         df <- merge(deseq_report[[i]], average_fpkms[,v], by='row.names')
         colnames(df) <- gsub('Row.names', 'id', colnames(df))
         colnames(df) <- gsub(".*WT", 'wt_fpkm', colnames(df))
         colnames(df) <- gsub(".*KO", 'ko_fpkm', colnames(df))
         return(df)
           })
write.table(deseq_report_fpkm[[3]], file="/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/bowtie/test/WT_EpivsKO_Epi.fpkms.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(deseq_report_fpkm[[1]], file="/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/bowtie/test/WT_2ivsKO_2i.fpkms.txt", row.names=F, col.names=T, sep="\t", quote=F)



library(rtracklayer)
source("~/source/Rscripts/annotation-functions.R")
annotationInfo('mm9')
genes <- import(genes_gtf, asRangedData=F)

markers <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/markers.txt", header=T, stringsAsFactors=F)
markers <- lapply(as.list(markers), function(x){
                  x <- x[!is.na(x)]
                  id <- mclapply(x, function(n) values(genes[values(genes)[['gene_name']] == n ])[['gene_id']] )
                  names(id) <- x
                  unlist(id)
})

require(reshape)# {{{
require(ggplot)
library(gridExtra)
plot_fpkm <- function(fpkms, markers, de, name, what='p2'){
    initial_order <- rownames(fpkms)
    source("~/source/Rscripts/ggplot-functions.R")
    fpkms <- melt(fpkms, id=rownames(fpkms))
    fpkms[grep("2lox|3Flox", fpkms[['X2']]), 'Group'] <- 'WT'
    fpkms[grep("spl2|3KO", fpkms[['X2']]), 'Group'] <- 'KO'
    fpkms[, 'Condition'] <- gsub(".*_", '', fpkms[['X2']])
    fpkms$X1 <- factor(fpkms$X1, levels=initial_order)

    for (x in levels(factor(fpkms$Condition))){
        fpkms[fpkms$Condition == x & fpkms$X1 %in% names(de[[x]][de[[x]]==TRUE]), 'Group'] <- paste0(fpkms[fpkms$Condition == x & fpkms$X1 %in% names(de[[x]][de[[x]]==TRUE]), 'Group'], '_de')
    }
    fpkms[['X1']] <- unlist(mclapply(as.character(fpkms[['X1']]), function(x) names(markers[markers == x])))

    fpkms[, 'Condition'] <- gsub("2i", 'ESCs 2i', fpkms[['Condition']])
    fpkms[, 'Condition'] <- gsub("Lif", 'ESCs SL', fpkms[['Condition']])
    fpkms[, 'Condition'] <- gsub("Epi", 'EpiSCs', fpkms[['Condition']])

    #gg <- c(gg_color_hue(2), c("#A1A1A1","#C9C9C9"))
    gg <- c(c("#AC2A21","#007378"), c("#FF9087","#1AD9DE"))
    names(gg)<- c('KO_de', 'WT_de', 'KO', 'WT')
    group_color <- fpkms$Group
    group_color <- as.character(gg[match(group_color,names(gg))])
    fpkms$Group <- factor(fpkms$Group, levels=c('WT','KO', 'WT_de', 'KO_de'))

    if(what == 'p1'){
        p <- ggplot(fpkms, aes(x=Condition, y=value, color=Group)) + geom_point() + scale_color_manual(values=gg)
        p <- p + geom_errorbar(stat = "hline", yintercept = "mean", width=0.8,aes(ymax=..y..,ymin=..y..))
        p <- p + facet_grid(~ X1,) + theme_bw() + xlab('') + ylab("FPKM")
        p <- p + theme(panel.grid.major.x= element_blank(),
                       strip.text.x = element_text(size = 10, colour = "black", angle = 90),
                       axis.text.x=element_text(angle=90, vjust=1)) + ggtitle(name)
        print(p)
   } else if (what == 'p2') {
        fpkms[as.character(fpkms[['X1']]) %in% names(markers)[1:3], 'State'] <- 'Pluripotency'
        fpkms[as.character(fpkms[['X1']]) %in% names(markers)[4:7], 'State'] <- 'Naive Pl.'
        fpkms[as.character(fpkms[['X1']]) %in% names(markers)[8:10], 'State'] <- 'Primed Pl.'
        fpkms[as.character(fpkms[['X1']]) %in% names(markers)[11:15], 'State'] <- 'Endoderm'
        fpkms[as.character(fpkms[['X1']]) %in% names(markers)[c(16:17,19)], 'State'] <- 'Ectoderm'
        fpkms[as.character(fpkms[['X1']]) %in% names(markers)[18], 'State'] <- 'Mesoderm'
        fpkms[as.character(fpkms[['X1']]) %in% names(markers)[c(20,21)], 'State'] <- 'JAK/STAT'
        
        p <- ggplot(fpkms, aes(x=X1, y=value, color=Group)) + geom_point(size=2) + scale_color_manual(values=gg)
        p <- p + geom_errorbar(stat = "hline", yintercept = "mean", width=0.8,aes(ymax=..y..,ymin=..y..))
        p <- p + facet_grid(Condition ~ State, scales="free")
        p <- p + scale_y_continuous(breaks= c(seq(0, 100, by=20), seq(100, 800, by=200)))
        p <- p + theme_bw() + xlab('') + ylab("FPKM")
        p <- p + theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0),
                       axis.text.y = element_text(size=7, face="plain"),
                       panel.grid.major.x = element_blank(),
                       axis.text.x=element_text(angle=90, vjust=1, face="plain")) + ggtitle(name)
        print(p)
    }

}# }}}


source("~/source/Rscripts/functions.R")
de <- list('2i'= deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_2ivsKO_2i.txt")$binary,
           'Lif' = deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_LifvsKO_Lif.txt")$binary,
           'Epi' = deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi.txt")$binary)


# using remcos de genes apart from epis
de <- list('2i'= deseq2vect("/nfs/research2/bertone/user/remco/hendrich/deseq/WT_2ivsKO_2i.txt", padj=0.1)$binary,
           'Lif' = deseq2vect("/nfs/research2/bertone/user/remco/hendrich/deseq/WT_LifvsKO_Lif.txt", padj=0.1)$binary,
           'Epi' = deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi.txt", padj=0.1)$binary)

lineage_marker <- c(markers$PL_general, markers$PL_naive, markers$PL_primed, markers$Endoderm, markers$Mesoderm, markers$Ectoderm, markers$JAK_STAT)
lineage_marker <- lineage_marker[c('Nanog', 'Pou5f1', 'Sox2', 'Klf2', 'Klf4',
             'Klf5', 'Zfp42', 'Lefty1', 'Lefty2', 'T',
             'Foxa2', 'Gata4', 'Gata6', 'Pdgfra', 'Eomes',
             'Otx2', 'Bmp4', 'Fgf5', 'Bmp4', 'Pax6', 'Lif', 'Lifr')]
lineage_marker <- lineage_marker[!duplicated(lineage_marker)]

pdf('fpkms.pdf')# {{{
condition <- unlist(mclapply(c('Epi'), function(i)grep(i,colnames(fpkms))))
for (x in names(markers)[6:14]) {
       f <- fpkms[as.character(markers[[x]]), condition]
       plot_fpkm(f, markers[[x]], de, x, what='p1')
}

condition <- unlist(mclapply(c('2i', '_Lif', 'Epi'), function(i)grep(i,colnames(fpkms))))
f <- fpkms[as.character(lineage_marker), condition]
plot_fpkm(f, lineage_marker, de, '', what='p2')

grid.newpage()
condition <- unlist(mclapply(c('2i', '_Lif'), function(i)grep(i,colnames(fpkms))))
f <- fpkms[as.character(lineage_marker[c('Sox2', 'Pou5f1', 'Nanog', 'Zfp42', 'Klf4', 'Klf5')]), condition]
plot_fpkm(f, lineage_marker[c('Sox2', 'Pou5f1', 'Nanog', 'Zfp42', 'Klf4', 'Klf5')], de, '', what='p1')
dev.off()# }}}

de <- deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi.txt")# {{{
plot_ma <- function(df, fdr=0.05){
    #df <- data.frame(df, FDR=ifelse(df$padj < fdr, paste("<", fdr, sep=" "), paste(">", fdr, sep=" ")))
    df <- subset(df, baseMean !=0)
    df <- data.frame(df, FDR=ifelse(df$padj < fdr, paste("<", fdr, sep=" "),
                                    ifelse(df$padj < fdr+0.05, paste("<", fdr+0.05, sep=" "), paste(">", fdr+0.05, sep=" "))))

    gg <- c("#F8766D", "#00BFC4", "#D4D4D4")
    names(gg) <- c(paste("<", fdr, sep=" "), paste("<", fdr+0.05, sep=" "), paste(">", fdr+0.05, sep=" "))
    # Define axis limits
    fc <- df$log2FoldChange
    ylim <- c(-1,1) * quantile(abs(fc[is.finite(fc)]), probs=0.99) * 1.1
    
    pch <- ifelse(df$log2FoldChange < ylim[1], 2, ifelse(df$log2FoldChange > ylim[2], 16, 20))
    df$pch <- as.factor(pch)

    p <- ggplot(df, aes(x=baseMean, y=log2FoldChange, colour=FDR, shape=pch)) + geom_point(size=1)
    p <- p + scale_x_continuous(trans=log_trans(), breaks=c(1e0,1e02,1e04,1e06)) + scale_y_continuous(limits = ylim)
    p <- p + scale_color_manual(values=gg) + scale_shape_manual(values=c(2,6,20), guide=FALSE)
    p <- p + geom_abline(intercept=0, slope=0, colour="#ff000080", size=2)
    p <- p + theme(legend.position= c(1,1),
                  legend.justification=c(1,1))
    p <- p + ylab(expression(log[2] ~ fold ~ change)) + xlab("mean of normalised counts")
    print(p)
}
pdf('plot_ma.pdf')
plot_ma(de$raw)
dev.off()# }}}

load("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/GOres.Rdata")
library(topGO)
plot_go_fpkm <- function(fpkms, markers, de, name){# {{{
    source("~/source/Rscripts/ggplot-functions.R")
    fpkms <- melt(fpkms, id=rownames(fpkms))
    fpkms[grep("2lox|3Flox", fpkms[['X2']]), 'Group'] <- 'WT'
    fpkms[grep("spl2|3KO", fpkms[['X2']]), 'Group'] <- 'KO'
    fpkms[, 'Condition'] <- gsub(".*_", '', fpkms[['X2']])

    for (x in levels(factor(fpkms$Condition))){
        fpkms[fpkms$Condition == x & fpkms$X1 %in% names(de[[x]][de[[x]]==TRUE]), 'Group'] <- paste0(fpkms[fpkms$Condition == x & fpkms$X1 %in% names(de[[x]][de[[x]]==TRUE]), 'Group'], '_de')
    }

    gg <- c(gg_color_hue(2), c("#964944","#047578"))
    names(gg)<- c('KO', 'WT', 'KO_de', 'WT_de')
    group_color <- fpkms$Group
    group_color <- as.character(gg[match(group_color,names(gg))])
    if(length(unique(fpkms$X1)) > 30) fpkms <- fpkms[fpkms$Group == 'WT_de' | fpkms$Group == 'KO_de' ,]
    if(length(unique(fpkms$X1)) > 30){
        sorted_de <- deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi.txt")$de
        sorted_de <- rownames(sorted_de[order(abs(-sorted_de$log2FoldChange)),])
        sorted_de <- sorted_de[1:30]

        print(head(fpkms))
        print(head(sorted_de))
        fpkms <- fpkms[as.character(fpkms$X1) %in% sorted_de,]
        fpkms <- fpkms[as.character(fpkms$X1) %in% c("ENSMUSG00000001288"),]
        print(head(fpkms))
    }
    levels(fpkms$X1) <- unlist(mclapply(levels(fpkms[['X1']]), function(x) names(markers[markers == x])))

    p <- ggplot(fpkms, aes(x=X1, y=value, color=Group)) + geom_point() + scale_color_manual(values=gg)
    p <- p + geom_errorbar(stat = "hline", yintercept = "mean", width=0.8,aes(ymax=..y..,ymin=..y..))
    p <- p + theme_bw() + xlab('') + ylab("FPKM") + theme(panel.grid.major.x = element_blank(),
                                                          panel.border = element_blank(),
                                                          axis.line=element_line())
    p <- p + theme(strip.text.x = element_text(size = 10, colour = "black", angle = 90),
                   axis.text.x=element_text(angle=90, vjust=1)) + ggtitle(name)
    print(p)
}# }}}

pdf('GO_enrichment.pdf')# {{{
for (x in names(GOres)) {
    go <- GOres[[x]]
    p <- ggplot(go$table,aes(y=Significant/Annotated,x=Term,fill=weight01Fisher)) + geom_bar(stat="identity") 
    p <- p + scale_fill_gradient(low="red",high="yellow")   + coord_flip() + ggtitle(x)
    print(p)

    condition <- unlist(mclapply(c('Epi'), function(i)grep(i,colnames(fpkms))))
    # for the 3 most significant terms return genes in category
    go_genes <- sapply(go$table[1,1], function(y) genesInTerm(go$GOdata, whichGO=y))
    names(go_genes) <- gsub("\\..*", '', names(go_genes))
    for (j in names(go_genes)){
        tmp_genes <- go_genes[[j]]
        names(tmp_genes) <- mclapply(tmp_genes, function(n) values(genes[values(genes)[['gene_id']] == n ])[['gene_name']] )
        f <- fpkms[tmp_genes, condition]
        plot_go_fpkm(f, tmp_genes, de, paste(x, go$table[go$table[,'GO.ID'] == j, 'Term'], sep="\n"))
    }
}
dev.off()
de <- deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi.txt")
up <- factor(as.integer(rownames(de$df) %in% rownames(de$up)))
names(up) <- rownames(de$df)
down <- factor(as.integer(rownames(de$df) %in% rownames(de$down)))
names(down) <- rownames(de$df)

go_up <- GOtermEnr(up, 'up-Epi', "BP", "org.Mm.eg.db")
go_down <- GOtermEnr(down, 'down-Epi', "BP", "org.Mm.eg.db")
save(go_up, go_down, file="GO-EpiSCs-up-down.Rdata")
# }}}

count_pseudogenes<-function(de, genes){
    de <- names(de[de==TRUE])
    table(values(genes[genes$gene_id %in% names(de[de==TRUE])])[['source']])
}

library(DESeq)
metadata <- read.delim(file="/nfs/research2/bertone/user/mxenoph/hendrich/hendrich_design.txt")

cds = newCountDataSetFromHTSeqCount(metadata, "")
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds.blind = estimateDispersions(cds, method = "blind")
vsd = varianceStabilizingTransformation(cds.blind)

pdf("PCA.pdf");# {{{
# Perform a principal components analysis on the transposed transformed expression data
pca = prcomp(t(exprs(vsd)))
# Group samples according to their "condition" class 

plot_data <- cbind(pca$x[,'PC1'], pca$x[,'PC2'])
group <- rep(NA, nrow(plot_data))
group[grep("2lox|3Flox", rownames(plot_data))] <- 'WT'
group[grep("spl2|3KO", rownames(plot_data))] <- 'KO'

condition <- rep(NA, nrow(plot_data))
condition[grep("2i", rownames(plot_data))] <- '2i + LIF'
condition[grep("_Lif", rownames(plot_data))] <- 'SL'
condition[grep("NoLif", rownames(plot_data))] <- 'No LIF'
condition[grep("N2B27", rownames(plot_data))] <- 'N2B27'
condition[grep("Epi", rownames(plot_data))] <- 'EpiSCs'

plot_data <- as.data.frame(plot_data)
colnames(plot_data) <- c('PC1','PC2')

plot_data$Group <- group
plot_data$Group <- factor(plot_data$Group)
plot_data$Condition <- condition
plot_data$Condition <- factor(plot_data$Condition)

p <- ggplot(plot_data, aes(x=PC1, y=PC2, colour=Group, shape=Condition))
#p <- p+ geom_text(aes(label=rownames(plot_data)), size=2.5, hjust=0.5, vjust=-0.5)
p <- p + geom_point(size=2.5)
p
dev.off()# }}}

# Calculate the Euclidean distances on the transposed exprs(vsd) matrix
dists = dist(t(exprs(vsd)));
# Load the library that includes the heatmap function
library(pheatmap);
# Plot the heatmap
pdf("EuclideanDistances_Heatmap.pdf");# {{{
library(RColorBrewer)
heat_cols <- colorRampPalette(rev(brewer.pal(9,"YlOrRd")))(255)

pheatmap(as.matrix(dists), #trace= "none", density.info= "none",
color=heat_cols,
cluster_cols=TRUE, cluster_rows=TRUE,
border_color=NA,
fontsize=10, fontsize_row=7, show_rownames=TRUE, show_colnames=TRUE)


heatmap.2(as.matrix(dists), # Convert the dists object into a matrix
          #trace="none", # can take values from "row", "column", "both" or "none". It indicates whether a solid trace line should be drawn.
          #density.info="none", # can take values "histogram", "density", or "none" and defines what type of plot should be plotted on the color-key.
          symm=TRUE, # a logical value indicating whether the matrix is symmetrical 
          cexRow=.5, cexCol=.5, # positive values defining the font size of row aand column labels 
          Colv=TRUE, # a logical value that determines if and how the "column dendrogram" should be reordered
          main="Distance Heatmap"); # the title of the plot
dev.off();# }}}



#read.table("../chip/config/modifications_update.conf", Header=T, sep="\t")
source("~/source/Rscripts/granges-functions.R")
mbd3 <-import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/macs/E12.M2.Epi.peaks.bed")
mbd3_gr <-macs2GRanges("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/macs/E12.M2.Epi_peaks.xls")

mi2b <-import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/macs/E12.Chd4.Epi.peaks.bed")
mi2b_gr <-macs2GRanges("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/macs/E12.Chd4.Epi_peaks.xls")

mi2b_ko <-import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/macs/G9.Chd4.Epi.peaks.bed")
mi2b_ko_gr <-macs2GRanges("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/macs/G9.Chd4.Epi_peaks.xls")
nurd <- mbd3[queryHits(findOverlaps(mbd3, mi2b))]
#nurd <- intersect(mbd3, mi2b)


promoters <- promoters(genes, upstream=2000, downstream=500)
genes_with_promoter <- resize(genes, width(genes)+2000, fix='end')

genes_width <- width(genes)
# genes that have genoby greater than 500 where promoter expands
genebody <- genes[(genes_width - 501) > 0]
genebody <- resize(genebody, width(genebody)-500, fix="end")
throw_out <- unique(queryHits(findOverlaps(genebody,promoters)))
keep <- !1:length(genebody) %in% throw_out
genebody <- genebody[keep]

load("../enhancers/enhancers.mm9.p300etal.Rdata")

throw_out <- unique(queryHits(findOverlaps(enh,promoters)))
throw_out <- unique(c(throw_out, unique(queryHits(findOverlaps(enh,genes)))))
keep <- !1:length(enh) %in% throw_out
intergenic_enhancers <- enh[keep]

nurd_annotation <- c(length(unique(subjectHits(findOverlaps(promoters, nurd)))),
                     length(unique(subjectHits(findOverlaps(genebody, nurd)))),
                     length(unique(subjectHits(findOverlaps(intergenic_enhancers, nurd))))
                     )
names(nurd_annotation) <- c('promoters','genebody','intergenic_enhancers')
nurd_annotation['intergenic'] <- length(nurd) - sum(nurd_annotation)
perc <- nurd_annotation/sum(nurd_annotation)

pdf('venn.pdf')
library(VennDiagram)
draw.pairwise.venn(# {{{
area1 = length(mbd3),
area2 = length(mi2b) - length(nurd),
cross.area = length(nurd),
scaled = TRUE,
category = c('Mbd3', 'Mi2b'),
fill = c("blue", "red"),
aplha = 0.5,
lty = "blank",
cex = 2.0,
cat.cex = 2.5,
cat.pos = c(285, 105),
cat.dist = 0.09,
cat.just = list(c(-1, -1), c(1, 1)),
ext.pos = 50,
ext.dist = 0.05,
ext.length = 0.85,
ext.line.lwd = 2,
ext.line.lty = "dashed")# }}}

nurd_annotation_df <- data.frame(nurd_annotation)
nurd_annotation_df$Feature <- c('Gene', 'Gene', 'Intergenic', 'Intergenic')
nurd_annotation_df$Type <- c('Promoter', 'Body', 'Enhancer', 'Intergenic')
nurd_annotation_df$percentage <- 100*(nurd_annotation_df$nurd_annotation/sum(nurd_annotation_df$nurd_annotation))

p <- ggplot(nurd_annotation_df, aes(x=Feature, y=percentage, fill=Type))
p <- p + geom_bar(stat="identity", width=.5) + scale_fill_grey(start = .4, end = .9, na.value = "grey50") + guides(fill= guide_legend(reverse=TRUE))
p <- p + ylab('% of NuRD bound regions (21765)')
p <- p + theme_bw() + theme(panel.grid.major.x = element_blank(),
                            panel.grid.major.y = element_blank(),
                            panel.background= element_blank(),
                            panel.border= element_blank(),
                            axis.line= element_line(),
                            axis.ticks= element_blank())
print(p)

dev.off()
# {{{
nurd_promoter <- values(promoters[unique(queryHits(findOverlaps(promoters, nurd)))])[['gene_id']]
nurd_genebody <- values(genes[unique(queryHits(findOverlaps(genes, nurd)))])[['gene_id']]

de <- deseq2vect("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi.txt")
de_genes <- rownames(de$de)
nurd_promoter_f <- factor(as.integer(values(genes)[['gene_id']] %in% nurd_promoter & values(genes)[['gene_id']] %in% de_genes))
names(nurd_promoter_f) <- values(genes)[['gene_id']]

go_nurd_promoter_de <- GOtermEnr(nurd_promoter_f, 'nurd_promoters', "BP", "org.Mm.eg.db")

nurd_promoter_f <- factor(as.integer(values(genes)[['gene_id']] %in% nurd_promoter ))
names(nurd_promoter_f) <- values(genes)[['gene_id']]
go_nurd_promoter <- GOtermEnr(nurd_promoter_f, 'nurd_promoters', "BP", "org.Mm.eg.db")

nurd_genebody <- factor(as.integer(values(genes)[['gene_id']] %in% nurd_genebody))
names(nurd_genebody) <- values(genes)[['gene_id']]
go_nurd_genebody <- GOtermEnr(nurd_genebody, 'nurd_promoters', "BP", "org.Mm.eg.db")

pdf('GO-bound-regions.pdf')
go <- go_nurd_promoter
p <- ggplot(go$table,aes(y=Significant/Annotated,x=Term,fill=weight01Fisher)) + geom_bar(stat="identity") 
p <- p + scale_fill_gradient(low="red",high="yellow")   + coord_flip() + ggtitle('NuRD bound promoters')
p
p %+% go_nurd_promoter_de$table + ggtitle('NuRD bound DE genes')
dev.off()

save(go_nurd_promoter, go_nurd_genebody, go_nurd_promoter_de, file="GO-EpiSCs-bound-pr-gb-de.Rdata")

mi2b_compensation <- nurd[unique(queryHits(findOverlaps(nurd,mi2b_ko_gr)))]
mi2b_ko_promoter <- values(promoters[unique(queryHits(findOverlaps(promoters, mi2b_compensation)))])[['gene_id']]
genes_with_promoter[genes_with_promoter$gene_id %in% mi2b_ko_promoter]


mi2b_ko_genebody <- values(promoters[unique(queryHits(findOverlaps(genebody, mi2b_compensation)))])[['gene_id']]
mi2b_ko_enhancer <- values(promoters[unique(queryHits(findOverlaps(intergenic_enhancers, mi2b_compensation)))])[['gene_id']]

throw_out <- unique(queryHits(findOverlaps(mi2b_ko_gr, mi2b)))
keep <- !1:length(mi2b_ko_gr) %in% throw_out
gain_mi2b <- mi2b_ko_gr[keep]
table(values(promoters[unique(queryHits(findOverlaps(promoters, gain_mi2b)))])[['gene_id']] %in% de_genes)

table(actively_expressed_genes %in% nurd_promoter)

# }}}

pdf("FC_EpiSCs.pdf")
plot(density(all_data[rownames(all_data) %in% nurd_promoter, 'log2FoldChange']),col="black",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="NuRD bound genes in Epi")
lines(density(all_data[rownames(all_data) %in% nurd_promoter & !is.na(all_data$padj) & all_data$padj < 0.05, 'log2FoldChange']),col="red3",lwd=2)
lines(density(all_data[rownames(all_data) %in% nurd_genebody,'log2FoldChange']),col="#666666",lwd=2)
lines(density(all_data[rownames(all_data) %in% nurd_genebody & !is.na(all_data$padj) & all_data$padj < 0.05, 'log2FoldChange']),col="orange",lwd=2)
legend("topright",legend=c("All","pval<0.05"),fill=c("blue","red"))
dev.off()


H3K27ac <- intersect(import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/chip/factor_2014/mm9/macs/mEpiSC_H3K27ac_peaks.bed"),
                     import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/chip/seki_2014/macs/ChIP_mEpi_H3K27ac.sort_peaks.bed"))

H3K4me3 <-intersect(import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/chip/factor_2014/mm9/macs/mEpiSC_H3K4me3_peaks.bed"),
                    import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/chip/seki_2014/macs/ChIP_mEpi_H3K4me3.sort_peaks.bed"))

active_genes_marks <- intersect(H3K27ac, H3K4me3)
nurd_at_active_sites <- reduce(nurd[queryHits(findOverlaps(nurd, active_genes_marks))])


active_genes_bound <- values(promoters[unique(subjectHits(findOverlaps(nurd_at_active_sites, promoters)))])[['gene_id']]

active_genes_bound_f <- factor(as.integer((rownames(average_fpkms) %in% active_genes_bound)))
names(active_genes_bound_f) <- rownames(average_fpkms)

active_bound <- data.frame(boolean=as.character(active_genes_bound_f), row.names=names(active_genes_bound_f))

epi_de <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/bowtie/test/WT_EpivsKO_Epi.fpkms.txt", 
           row.names='id', header=T, sep="\t")

out <- merge(epi_de, active_bound, by='row.names')
colnames(out) <- gsub('Row.names', 'id', colnames(out))


write.table(out, file="/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/bowtie/test/WT_EpivsKO_Epi.fpkms.binding.txt",
            row.names=F, col.names=T, sep="\t", quote=F)
