#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)

parser <-  ArgumentParser(description="Convert bed file to gtf")
parser$add_argument('-b', '--bed', metavar= "file", required='True', type= "character", help= "MACS generated bed file for called peaks")

args <- parser$parse_args()
# }}}}

# TODO: read file extension, check column number and whether have the right format and convert all plus bed files to gtf

x <- c('GenomicRanges', 'rtracklayer')
lapply(x, suppressMessages(library), character.only=T)

bed <- read.table(args$bed, sep="\t")
colnames(bed) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')[1:ncol(bed)]
levels(bed$strand) <- gsub('.', "*", levels(bed$strand))
bed$chromStart <- bed$chromStart - 1

bed <- with(bed, GRanges(chrom,
                         IRanges(chromStart, chromEnd),
                         strand,
                         peak = name,
                         fe = signalValue,
                         pValue = pValue,
                         qValue = qValue
                         ))



gtf <- paste0(gsub(file_ext(args$bed), '', args$bed), 'gtf')
export.gff(bed, gtf)

centered <- resize(bed,
                   ifelse(width(bed) < median(width(bed)), width(bed), median(width(bed)) ),
                   fix="center")

gtf <- paste0(gsub(file_ext(args$bed), '', args$bed), 'centered-median-window.gtf')
export.gff(bed, gtf)
