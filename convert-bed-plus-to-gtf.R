#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)

parser <-  ArgumentParser(description="Convert bed file to gtf")
parser$add_argument('-b', '--bed', metavar= "file", required='True', type= "character", help= "MACS generated bed file for called peaks")

args <- parser$parse_args()
# }}}}

# TODO: read file extension, check column number and whether have the right format and convert all plus bed files to gtf

x <- c('GenomicRanges', 'rtracklayer', 'plyr', 'tools')
lapply(x, suppressMessages(library), character.only=T)

parse_broadPeak <- function(broad){# {{{
    bed <- read.table(broad, sep="\t")
    colnames(bed) <- c('chrom', 'chromStart', 'chromEnd',
                       'name', 'score', 'strand',
                       'signalValue', 'pValue', 'qValue', 'peak')[1:ncol(bed)]
    
    # Change levels so that strand matches -. + and * for GRanges instead of -1, +1, .
    bed[['strand']] <- revalue(bed[['strand']], c("-1"="-", "+1" = "+", "." = "*"))
    bed[['chromStart']] <- bed[['chromStart']] - 1
    
    grange <- with(bed, GRanges(chrom,
                             IRanges(chromStart, chromEnd),
                             strand,
                             peak = name,
                             score = score,
                             fe = signalValue,
                             pValue = pValue,
                             qValue = qValue
                             ))

    if (ncol(bed) == 10) values(grange)[['summit']] <- with(bed, 'peak')
    
    gtf <- paste0(gsub(file_ext(broad), '', broad), 'gtf')
    export.gff(grange, gtf)

    centered <- resize(grange,
                       ifelse(width(grange) < median(width(grange)), width(grange), median(width(grange)) ),
                       fix="center")

    gtf <- paste0(gsub(file_ext(broad), '', broad), 'centered-median-window.gtf')
    export.gff(centered, gtf)
}# }}}

if ( file_ext(args$bed) %in% c('broadPeak', 'narrowPeak')) {
    parse_broadPeak(args$bed)
} else {
    warning("Can not convert bed file to gtf [unknown format].")
}

