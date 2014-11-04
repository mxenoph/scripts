#!/homes/mxenoph/local/bin/Rscript --vanilla

library(GenomicRanges)
library(rtracklayer)
library(pheatmap)
library(RColorBrewer)
source("~/source/Rscripts/annotation-functions.R")
annotationInfo('mm9')

# This assumes it's mouse data. If not this has to be changed
chr <- c(paste0('chr', 1:19), 'chrX', 'chrY')

# get protein coding genes
genes <- import(genes_gtf, asRangedData=F)
genes_all <- import(genes_gtf, asRangedData=F)
genes <- subset(genes, seqnames(genes) %in% chr & source == 'protein_coding')

# slow -- this should be done once only
#txdb <- makeTranscriptDbFromGFF(gtf, format='gtf')
#test <- genes(txdb)

oct4 <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/chip/files-oct4-peaks.txt", header=F, col.names=c('sample', 'file'))
mbd3 <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/chip/files-mbd3-peaks.txt", header=F, col.names=c('sample', 'file'))

oct4_gr <-  mclapply(as.vector(oct4[,2]), function(file){
                     bed <- import.bed(file)
                     bed <- subset(bed, seqnames(bed) %in% chr)})
oct4_gr <- GRangesList(oct4_gr)
names(oct4_gr) <- as.character(oct4[,1])

mbd3_gr <-  mclapply(as.vector(mbd3[,2]), function(file){
                     bed <- import.bed(file)
                     bed <- subset(bed, seqnames(bed) %in% chr)})
mbd3_gr <- GRangesList(mbd3_gr)
names(mbd3_gr) <- as.character(mbd3[,1])

oct4_all <- reduce(unlist(oct4_gr))
names(oct4_all) <- paste(seqnames(oct4_all), start(oct4_all), end(oct4_all), sep=":")
mbd3_all <- reduce(unlist(mbd3_gr))
names(mbd3_all) <- paste(seqnames(mbd3_all), start(mbd3_all), end(mbd3_all), sep=":")

mat <- matrix(, nrow = length(oct4_all), ncol = length(oct4_gr) + length(mbd3_gr))
colnames(mat) <- c(names(oct4_gr), names(mbd3_gr))
rownames(mat) <- names(oct4_all)

for (x in names(oct4_gr)){
    # returns true/false list that we convert to integer
    mat[,x] <- as.integer(overlapsAny(oct4_all, oct4_gr[[x]]))
}

for (x in names(mbd3_gr)){
    mat[,x] <- as.integer(overlapsAny(oct4_all, mbd3_gr[[x]]))
}
oct4_gain <- mat[,'oct4-2i-2lox'] == 0 & mat[,'mbd3-2i'] == 1 & mat[,'mbd3-EpiSCs'] == 0  & mat[,'oct4-EpiLCs'] == 1
oct4_loss <- mat[,'oct4-2i-2lox'] == 1 & mat[,'mbd3-2i'] == 0 & mat[,'mbd3-EpiSCs'] == 1  & mat[,'oct4-EpiLCs'] == 0

# check if the peaks overlap genes etc 
gain <- genes[subjectHits(findOverlaps(oct4_all[oct4_gain], genes)), 'gene_id']
loss <- elementMetadata(genes[unique(subjectHits(findOverlaps(oct4_all[oct4_loss], genes)))])[['gene_id']]
gain <- elementMetadata(genes[unique(subjectHits(findOverlaps(oct4_all[oct4_gain], genes)))])[['gene_id']]

oct4_over_gene <- as.integer(overlapsAny(oct4_all, genes))
names(oct4_over_gene) <- names(oct4_all)

dev_time <- unlist(mclapply(c('0hr', '2i', 'SL', '16hr', 'EpiLCs', 'EpiSCs'), function(x){colnames(mat)[grep(x, colnames(mat))]}))
pdata <- mat[,dev_time]
pdata <- pdata[do.call(order, as.data.frame(pdata)),]
pdata <- cbind(pdata, 'overlaps_gene'=as.integer(oct4_over_gene[rownames(pdata)]))

cool_cols <- colorRampPalette(c('aliceblue','darkcyan'))(100)

pdf("plots/oct4matrix.pdf")
pheatmap(pdata, #trace= "none", density.info= "none",
         color=cool_cols,
         clustering_distance_cols= "binary", clustering_method= "ward",
         cluster_cols=FALSE, cluster_rows=FALSE,
         border_color=NA,
         fontsize=10, fontsize_row=7, show_rownames=FALSE, show_colnames=TRUE)
dev.off()

pdata_rev <- mat[,rev(dev_time)]
pdata_rev <- pdata_rev[do.call(order, as.data.frame(pdata_rev)),]
pdata_rev <- cbind(pdata_rev, 'overlaps_gene'=as.integer(oct4_over_gene[rownames(pdata_rev)]))


pdf("plots/oct4matrix-rev.pdf")
pheatmap(pdata_rev, #trace= "none", density.info= "none",
         color=cool_cols,
         clustering_distance_cols= "binary", clustering_method= "ward",
         cluster_cols=FALSE, cluster_rows=FALSE,
         border_color=NA,
         fontsize=10, fontsize_row=7, show_rownames=FALSE, show_colnames=TRUE)
dev.off()

source("~/local/functions.R")
library(topGO)
go <- factor(as.integer(elementMetadata(genes)[['gene_id']] %in% gain))
names(go) <- elementMetadata(genes)[['gene_id']]

go <- GOtermEnr(go, 'oct4-gain', "BP", "org.Mm.eg.db")

library(gridExtra)

pdf(file.path('plots',"oct-gain_GO.pdf"), paper='a4')
p <- ggplot(go$table,aes(y=Significant/Annotated,x=Term,fill=weight01Fisher)) + geom_bar(stat="identity") + scale_fill_gradient(low="red",high="yellow")   + coord_flip()
p
dev.off()

markers <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/markers.txt", header=T, stringsAsFactors=F)
markers <- lapply(as.list(markers), function(x){
                  x <- x[!is.na(x)]
                  id <- mclapply(x, function(n) values(genes[values(genes)[['gene_name']] == n ])[['gene_id']] )
                  names(id) <- x
                  unlist(id)
         })
