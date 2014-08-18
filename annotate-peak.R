require('GenomicRanges')
require('GenomicFeatures')
require('plyr') # for revaluing factors
require('ShortRead')
require('Rsamtools')

setwd('/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/macs')
peakFiles <- list.files(pattern="*Epi*.peaks.bed")

peaks <- lapply(peakFiles, read.table, col.names=c('chr', 'start', 'end', 'peakNo', 'score'))
names(peaks) <- gsub("\\.Epi.+", '', peakFiles)
peaks <- lapply(names(peaks), function(x){
                with(peaks[[x]], GRanges(chr, IRanges(start, end), strand="*", score=score))
})
names(peaks) <- gsub("\\.Epi.+", '', peakFiles)

#load annotation--gives my.annotation object with $genes, $exons etc
load("mm9.e67.ann.Rdata")
ext.up <- 5000
ext.down <- 500
ext.genes <- mm9.e67.ann$genes

#These produce warning messages that the ranges have been trimmed to be <=seqlengths
ext.genes <- suppressWarnings(resize(ext.genes, width(ext.genes) + ext.up, fix="end")) #For flanking the regions upstream -- works ok for genes on the minus strand
ext.genes <- suppressWarnings(resize(ext.genes, width(ext.genes) + ext.down, fix="start")) #For flanking downstream

bound <- lapply(names(peaks), function(x){
                ov <- findOverlaps(mm9.e67.ann$genes, peaks[[x]])
                id <- unique(names(mm9.e67.ann$genes)[queryHits(ov)])
})
names(bound) <- names(peaks)

bound.df <- do.call(cbind, lapply(names(bound), function(x){
    as.integer(names(mm9.e67.ann$genes) %in% bound[[x]])
}))
colnames(bound.df) <- paste(names(bound), 'bound', sep='.')
rownames(bound.df) <- names(mm9.e67.ann$genes)
bound.df<-as.data.frame(bound.df)


#Are the extended genes +500 -5000bp bound?
ext.bound <- lapply(names(peaks), function(x){
                ov <- findOverlaps(ext.genes, peaks[[x]])
                id <- unique(names(ext.genes)[queryHits(ov)])
})
names(ext.bound) <- names(peaks)

ext.bound.df <- do.call(cbind, lapply(names(ext.bound), function(x){
    as.integer(names(ext.genes) %in% ext.bound[[x]])
}))
colnames(ext.bound.df) <- paste(names(ext.bound), 'bound', sep='.')
rownames(ext.bound.df) <- names(ext.genes)
ext.bound.df<-as.data.frame(ext.bound.df)
######

#Is NuRD.bound bound= chd4 + mbd3
bound.df[, 'NuRD.bound'] <- 0
bound.df[bound.df$E12.Chd4.bound == 1 & bound.df$E12.M2.bound == 1, 'NuRD.bound'] <- 1

ext.bound.df[, 'NuRD.bound'] <- 0
ext.bound.df[ext.bound.df$E12.Chd4.bound == 1 & ext.bound.df$E12.M2.bound == 1, 'NuRD.bound'] <- 1

###
de.analysis <- read.table("/nfs/nobackup2/research/bertone/mxenoph/bh/forMM9/deseq/WT_EpivsKO_Epi.txt", header=TRUE, row.names=1, sep="\t")
my.data <- merge(bound.df, de.analysis[, c('log2FoldChange', 'padj')], by="row.names")
rownames(my.data) <- my.data[['Row.names']]
my.data <- my.data[, !names(my.data) %in% 'Row.names']

#set FC for not expressed genes to 0
my.data[is.na(my.data$log2FoldChange), 'log2FoldChange'] <- 0

###
ext.data <- merge(ext.bound.df, de.analysis[, c('log2FoldChange', 'padj')], by="row.names")
rownames(ext.data) <- ext.data[['Row.names']]
ext.data <- ext.data[, !names(ext.data) %in% 'Row.names']
ext.data[is.na(ext.data$log2FoldChange), 'log2FoldChange'] <- 0


#remco's data
serum.data <- read.table("/nfs/research2/bertone/user/remco/hendrich/chipseq/peaks/more/all.data.txt",header=T, row.names=1, as.is=TRUE,sep="\t")
twoI.data <- read.table("/nfs/research2/bertone/user/remco/hendrich/chipseq/peaks/2i/all.data_2i.txt",header=T, row.names=1, as.is=TRUE,sep="\t")

serum <- serum.data['nurd']
serum[serum == "not bound"] <- 0
serum[serum == "bound"] <- 1
twoI <- twoI.data['nurd']
twoI[twoI == "not bound"] <- 0
twoI[twoI == "bound"] <- 1

esc <- merge(serum, twoI, by="row.names")
rownames(esc) <- esc[['Row.names']]
esc <- esc[, !names(esc) %in% 'Row.names']
colnames(esc) <- c('NuRD.bound.serum', 'NuRD.bound.2i')

all.data <- merge(my.data, esc, by="row.names")
all.data[,'esc.not.bound'] <- 0
all.data[all.data$NuRD.bound == 1 & (all.data$NuRD.bound.serum ==0 | all.data$NuRD.bound.2i ==0), 'esc.not.bound'] <- 1
all.data[,'epi.not.bound'] <- 0
all.data[all.data$NuRD.bound == 0 & (all.data$NuRD.bound.2i ==1 | all.data$NuRD.bound.serum ==1),  'epi.not.bound'] <- 1
all.data[,'all.bound'] <- 0
all.data[all.data$NuRD.bound == 1 & (all.data$NuRD.bound.2i ==1 | all.data$NuRD.bound.serum ==1),  'all.bound'] <- 1

pdf("FCuponKO.Epi.NuRD.bound.pdf")
plot(density(my.data$log2FoldChange[my.data$NuRD.bound==1]),col="blue",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="NuRD bound genes")
lines(density(my.data$log2FoldChange[my.data$NuRD.bound==1 & !is.na(my.data$padj) & my.data$padj < 0.05]),col="red",lwd=2)
legend("topright",legend=c("All","pval<0.05"),fill=c("blue","red"))
dev.off()


pdf("FCuponKO.Epi.NuRD.bound.ext.pdf")
plot(density(ext.data$log2FoldChange[ext.data$NuRD.bound==1]),col="blue",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="NuRD bound genes")
lines(density(ext.data$log2FoldChange[ext.data$NuRD.bound==1 & !is.na(ext.data$padj) & ext.data$padj < 0.05]),col="red",lwd=2)
legend("topright",legend=c("All","pval<0.05"),fill=c("blue","red"))
dev.off()


pdf("FCuponKO.Epi.comp.NuRD.binding.pdf")
plot(density(all.data$log2FoldChange[all.data$all.bound==1]),col="blue",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="NuRD bound genes")
lines(density(all.data$log2FoldChange[all.data$all.bound==1 & !is.na(all.data$padj) & all.data$padj < 0.05]),col="red",lwd=2)

lines(density(all.data$log2FoldChange[all.data$esc.not.bound==1]), col="black",lwd=2)
lines(density(all.data$log2FoldChange[all.data$esc.not.bound==1 & !is.na(all.data$padj) & all.data$padj < 0.05]),col="green",lwd=2)

lines(density(all.data$log2FoldChange[all.data$epi.not.bound==1]), col="purple",lwd=2)
lines(density(all.data$log2FoldChange[all.data$epi.not.bound==1 & !is.na(all.data$padj) & all.data$padj < 0.05]),col="orange",lwd=2)
legend("topright",legend=c("All.bound","pval <0.05", "only.epi.bound", 'only.epi.pval', 'only.esc.bound', 'only.esc.pval.'),fill=c("blue","red", "black", "green", "purple", "orange"))
dev.off()




## first get gene information from UCSC. This will take 2 minutes or so.
## Or if you have these saved you can load them in using: refGene.hg18=loadFeatures("refGene.hg18.GR.sqlite")
#refGene.mm9 <- makeTranscriptDbFromUCSC(genom="mm9",tablename="refGene")
#tx <- transcripts(refGene.mm9)
#what=c("rname", "strand", "pos", "qwidth")
#TSS.counts=NULL
#param=ScanBamParam(what = what)
#
#strands=strand(tx)
#ix.gene.minus=which(strands=="-")
#TSS=start(tx)
#TSS[ix.gene.minus]=end(tx)[ix.gene.minus]
#ext=1000
#TSS.GRanges=GRanges(seqnames=seqnames(tx),
#                        ranges=IRanges(start=TSS-ext, end=TSS+ext))
#bam=scanBam("/nfs/nobackup2/research/bertone/mxenoph/bh/ChIP-Epi/mm9/bowtie/E12EpiM2.sorted.bam", param=param)[[1]]
#IRange.reads=GRanges(seqnames=Rle(bam$rname),
#                         ranges=IRanges(bam$pos, width=bam$qwidth))
#counts.m2=countOverlaps(tx, IRange.reads)
#

