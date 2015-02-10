#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)

parser <-  ArgumentParser(description="Convert bed file to gtf")
parser$add_argument('-b', '--bed', metavar= "file", required='True', type= "character", help= "MACS generated bed file for called peaks")
parser$add_argument('-a', '--assembly', required='True', type= "character", help= "assembly")

args <- parser$parse_args()
# }}}}

# TODO: read file extension, check column number and whether have the right format and convert all plus bed files to gtf

x <- c('GenomicRanges', 'rtracklayer', 'plyr', 'tools')
lapply(x, suppressMessages(library), character.only=T)
source("~/source/Rscripts/annotation-functions.R")

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

    # TODO: add a check and if narrowPeak with summit provided center around summit
    centered <- resize(grange,
                       ifelse(width(grange) < median(width(grange)), width(grange),
                              # round_any is in plyr package and will round down to the previous hundred
                              round_any(median(width(grange)), 100, f = floor)),
                       fix="center")

    gtf <- paste0(gsub(file_ext(broad), '', broad), 'centered-median-window.gtf')
    export.gff(centered, gtf)

    # centering around the left boundary of the peak
    left_boundary <- flank(resize(grange, 500, fix="start"), 500, "start")
    right_boundary <- resize(resize(grange, 500, fix="end"), 1000, "start")

    gtf <- paste0(gsub(file_ext(broad), '', broad), 'left-boundary-window.gtf')
    export.gff(left_boundary, gtf)
    
    gtf <- paste0(gsub(file_ext(broad), '', broad), 'right-boundary-window.gtf')
    export.gff(right_boundary, gtf)
}
# }}}

parse_gappedPeak <- function(gapped){# {{{
    # contains both the broad regiona nd narrow peak
    bed <- read.table(gapped, sep="\t")
    colnames(bed) <- c('chrom', 'chromStart', 'chromEnd',
                       'name', 'ucsc_qValue', 'strand',
                       'narrow_start', 'narrow_end', 'rgb_col',
                       'blocks', 'block_length', 'block_start',
                       'fe', 'pValue', 'qValue')
    
    # Change levels so that strand matches -. + and * for GRanges instead of -1, +1, .
    bed[['strand']] <- revalue(bed[['strand']], c("-1"="-", "+1" = "+", "." = "*"))
    bed[['chromStart']] <- bed[['chromStart']] - 1
    
    grange <- with(bed, GRanges(chrom,
                             IRanges(chromStart, chromEnd),
                             strand,
                             peak = name,
                             score = qvalue,
                             fe = fe,
                             pValue = pValue,
                             qValue = qValue
                             ))

    if (ncol(bed) == 10) values(grange)[['summit']] <- with(bed, 'peak')
    
    gtf <- paste0(gsub(file_ext(broad), '', broad), 'gtf')
    export.gff(grange, gtf)

    # TODO: add a check and if narrowPeak with summit provided center around summit
    centered <- resize(grange,
                       ifelse(width(grange) < median(width(grange)), width(grange),
                              # round_any is in plyr package and will round down to the previous hundred
                              round_any(median(width(grange)), 100, f = floor)),
                       fix="center")

    gtf <- paste0(gsub(file_ext(broad), '', broad), 'centered-median-window.gtf')
    export.gff(centered, gtf)

    # centering around the left boundary of the peak
    left_boundary <- flank(resize(grange, 500, fix="start"), 500, "start")
    right_boundary <- resize(resize(grange, 500, fix="end"), 1000, "start")

    gtf <- paste0(gsub(file_ext(broad), '', broad), 'left-boundary-window.gtf')
    export.gff(left_boundary, gtf)
    
    gtf <- paste0(gsub(file_ext(broad), '', broad), 'right-boundary-window.gtf')
    export.gff(right_boundary, gtf)
}
# }}}


parse_bed5 <- function(bed5){# {{{
    # contains both the broad regiona nd narrow peak
    bed <- read.table(bed5, sep="\t")
    colnames(bed) <- c('chrom', 'chromStart', 'chromEnd',
                       'name', 'pValue')
    
    # Change levels so that strand matches -. + and * for GRanges instead of -1, +1, .
    bed[['chromStart']] <- bed[['chromStart']] - 1
    
    grange <- with(bed, GRanges(chrom,
                             IRanges(chromStart, chromEnd),
                             strand = "*",
                             peak = name,
                             pValue = pValue))
    
    gtf <- paste0(gsub(file_ext(bed5), '', bed5), 'gtf')
    export.gff(grange, gtf)

    centered <- resize(grange,
                       ifelse(width(grange) < median(width(grange)), width(grange),
                              # round_any is in plyr package and will round down to the previous hundred
                              round_any(median(width(grange)), 100, f = floor)),
                       fix="center")

    gtf <- paste0(gsub(file_ext(bed5), '', bed5), 'centered-median-window.gtf')
    export.gff(centered, gtf)

    # centering around the left boundary of the peak
    left_boundary <- flank(resize(grange, 500, fix="start"), 500, "start")
    right_boundary <- resize(resize(grange, 500, fix="end"), 1000, "start")

    gtf <- paste0(gsub(file_ext(bed5), '', bed5), 'left-boundary-window.gtf')
    export.gff(left_boundary, gtf)
    
    gtf <- paste0(gsub(file_ext(bed5), '', bed5), 'right-boundary-window.gtf')
    export.gff(right_boundary, gtf)

    bound <- annotate_range(grange, genes, 'genes', return_att=T)
    x <- genes[seqnames(genes) %in% seqlevels(genes)[grep('chr', seqlevels(genes))]]
    df <- data.frame(gene_id = names(x), boolean = as.numeric(names(x) %in% bound$attribute))
    file <- gsub(file_ext(bed5), 'bound_genes.txt',  bed5)
    write.table(df, file, sep="\t", quote=F, row.names=F)

    # TODO: Function does not use label anywhere so update both function and calls maybe
    bound_promoter <- annotate_range(grange, tss_window, 'tss_window', return_att=T)
    x <- tss_window[seqnames(tss_window) %in% seqlevels(tss_window)[grep('chr', seqlevels(tss_window))]]
    df <- data.frame(gene_id = names(x), boolean = as.numeric(names(x) %in% bound_promoter$attribute))
    file <- gsub(file_ext(bed5), 'bound_promoters.txt',  bed5)
    write.table(df, file, sep="\t", quote=F, row.names=F)
}
# }}}

get_annotation(args$assembly)
if ( file_ext(args$bed) %in% c('broadPeak', 'narrowPeak')) {
    parse_broadPeak(args$bed)
} else {
    warning("Can not convert bed file to gtf [unknown format].")
}

