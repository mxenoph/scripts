#!/usr/bin/env Rscript

# Parse arguments # {{{
library(argparse)
source("~/source/Rscripts/annotation-functions.R")

parser <-  ArgumentParser(description="Annotate any range. The gtf passed could be a default ensembl gtf or a custom one")
parser$add_argument('-b', '--bed', metavar= "file", required='True', type= "character", help= "BED file")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args <- parser$parse_args()

output_path <- file.path(args$out, 'annotated')
plot_path <- file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

# function in annotation-functions.R
# returns chr_size and genes_gtf directly in the current environment
annotationInfo(tolower(args$assembly))

# define variables
# }}}

# Import libraries & Turn off warning messages for loading the packages-globally# {{{
options(warn = -1)
# Required for file_path_sans_ext()
library(tools)
library(dplyr)
# IMPORTANT: if Genomic Features and hence AnnotatioDbi are needed the select 
# function from dplyr is masked and hence gives error

#library(GenomicRanges)
#library(GenomicFeatures)
#library(plyr) # for revaluing factors
#library(ShortRead)
#library(Rsamtools)
options(warn = 0)# }}}

map <- base::Map
# Cause something masks these functions and nothing works after
summarise <- dplyr::summarise
select <- dplyr::select

let <- function (.expr, ...)
    eval(substitute(.expr), list2env(list(...), parent = parent.frame()))

parseBedtools <- function(annotated){

    # Keep only gene name from attribute field
    extractAttribute <- function(name, header){
        gsub(';', '', gsub(name, '', let(match = regexpr(sprintf('%s\\S+', name), header, perl = TRUE),
                                         as.character(regmatches(header, match)))))
    }

    readout <- read.table(annotated,
                          col.names = c(paste('bed',
                                              c('chr', 'start', 'end', 'name', 'score'),
                                              sep = '_'),
                                        paste('gtf',
                                              c('chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
                                              sep = '_'),
                                        'closest_distance'
                                        ),
                          sep='\t')

#    genes <- unlist(map(extractAttribute, readout[['gtf_attribute']], name='^ gene_id '))
#    trx_name <- unlist(map(extractAttribute, readout[['gtf_attribute']], name=' transcript_name '))

    # Converted to dplyr dataframe to make things faster and easier
    readout_df <- tbl_df(readout)

    genes <- unname(unlist(map(extractAttribute, select(readout_df, gtf_attribute), name='^ gene_id ')))
    # might be good idea to remove this and make it more generic for all types of gtf provided to bedtoolsClosest
    gene_name <- unname(unlist(map(extractAttribute, select(readout_df, gtf_attribute), name=' gene_name ')))

    # %.% is a chain operator. For more detailed explanation see 
    # http://stackoverflow.com/questions/22314680/how-to-use-the-operator-in-r>
    readout_df <- readout_df %.%
                    select(bed_chr:bed_score, gtf_feature, gtf_strand, closest_distance) %.%
                    mutate(gene = genes,
                           gene_name = gene_name)
    # Mark as targets genes that are 2Kb upstream the features start and 1Kb downstream
    readout_df <- readout_df %.% mutate(target = ifelse(closest_distance >= -2000 & closest_distance <= 1000, 1, 0))

    readout_df %.%
    group_by(bed_name) %.%
    select(gene)



}

bedtoolsClosest <- function(bed, annotation, output_path){# {{{
    output <- paste0(basename(file_path_sans_ext(bed)),
                     '-annotated',
                     '.txt')
    
    system(sprintf('bedtools closest -D "b" -a %s -b %s > %s',
               bed,
               annotation, # coming from the annotationInfo function
               file.path(output_path, output)
               ))
    # Example output; last column is the distance from the feature -- here across multiple lines for readability
    # chr1	3373072 3373608	MACS_peak_1	128.12	\n
    # chr1	protein_coding  exon	3411783	3411982	.	-	.	\n
    # gene_id "ENSMUSG00000051951"; transcript_id "ENSMUST00000070533"; exon_number "2"; gene_name "Xkr4"; gene_biotype "protein_coding"; transcript_name "Xkr4-001";	38175
    # IMPORTANT: Negative distances report upstream features i.e. peak is upstream of gene
    parseBedtools(file.path(output_path, output))
}# }}}

bedtoolsClosest(args$bed, genes_gtf, output_path)


# Old code # {{{
old <- function(){
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
}
# }}}
