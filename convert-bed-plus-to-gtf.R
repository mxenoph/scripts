#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)

parser <-  ArgumentParser(description="Convert bed file to gtf")
parser$add_argument('-b', '--bed', metavar= "file", required='True', type= "character", help= "MACS generated bed file for called peaks")
parser$add_argument('-a', '--assembly', required='True', type= "character", help= "assembly")

args <- parser$parse_args()
# }}}}

x <- c('GenomicRanges', 'rtracklayer', 'plyr', 'tools')
lapply(x, suppressMessages(library), character.only=T)
source("~/source/Rscripts/annotation-functions.R")
# Get number of cores
processes <- as.numeric(system("echo $LSB_DJOB_NUMPROC", intern=T))

center <- function(grange) {# {{{
    # TODO: add a check and if narrowPeak with summit provided center around summit
    centered <- resize(grange,
                       ifelse(width(grange) < median(width(grange)), width(grange),
                              # round_any is in plyr package and will round down to the previous hundred
                              round_any(median(width(grange)), 100, f = floor)),
                       fix="center")
    return(centered)
}# }}}

center_around_boundary <- function(grange){# {{{
    # centering around the left boundary of the peak
    left_boundary <- flank(resize(grange, 500, fix="start"), 500, "start")
    right_boundary <- resize(resize(grange, 500, fix="end"), 500, "start")
    return(list('left'=left_boundary, 'right'=right_boundary))
}
# }}}

export_gtfs <- function(file, grange){# {{{
    gtf <- paste0(gsub(file_ext(file), '', file), 'gtf')
    export.gff(grange, gtf)

    gtf <- paste0(gsub(file_ext(file), '', file), 'centered-median-window.gtf')
    export.gff(center(grange), gtf)

    boundary_centered <- center_around_boundary(grange)
    gtf <- paste0(gsub(file_ext(file), '', file), 'left-boundary-window.gtf')
    export.gff(boundary_centered[['left']], gtf)
    
    gtf <- paste0(gsub(file_ext(file), '', file), 'right-boundary-window.gtf')
    export.gff(boundary_centered[['right']], gtf)
}# }}}

export_bound <- function(feature, label, grange, file){# {{{
    bound <- annotate_range(grange, feature, label, return_att=T)
    f_nocontig <- feature[seqnames(feature) %in% seqlevels(feature)[grep('chr', seqlevels(feature))]]
    
    bound_table <- data.frame(gene_id = names(f_nocontig), bound = as.numeric(names(f_nocontig) %in% bound$attribute))
    # Replacing bound with the marks name
    colnames(bound_table) <- gsub('bound', basename(args$bed), colnames(bound_table))

    compiled <- file.path(dirname(file), paste0('bound_', label, '.txt'))
    if (file.exists(compiled)) {

        existing_table <- read.table(compiled, header=T, sep="\t")
        flag <- FALSE
        if ( length(levels(existing_table$gene_id)) != length(levels(bound_table$gene_id)) ){
            warning("compiled of bound genes already exists and number of genes does not match.
                    \nReducing to genes in both lists. Proceed with caution"
                    )
            flag <- TRUE
        }

        matches <- match(bound_table$gene_id, existing_table$gene_id)
        m <- matches[!is.na(matches)]
        existing_table[m, eval(colnames(bound_table)[2])] <- as.character(bound_table[!is.na(matches), 2])
        # If match not found in existing table it will put NA and here any row with NA is excluded
        # In essence this keeps only the genes present in both 'annotations'
        existing_table <- existing_table[complete.cases(existing_table), ]
        
        write.table(existing_table, compiled, sep="\t", quote=F, row.names=F)
    } else {
        write.table(bound_table, compiled, sep="\t", quote=F, row.names=F)
    }

    return(bound_table)
    # TODO: Function does not use label anywhere so update both function and calls maybe
#    bound_promoter <- annotate_range(grange, tss_window, 'tss_window', return_att=T)
#    x <- tss_window[seqnames(tss_window) %in% seqlevels(tss_window)[grep('chr', seqlevels(tss_window))]]
#    df <- data.frame(gene_id = names(x), boolean = as.numeric(names(x) %in% bound_promoter$attribute))
#    file <- gsub(file_ext(file), 'bound_promoters.txt',  file)
#    write.table(df, file, sep="\t", quote=F, row.names=F)
}
# }}}

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
    return(grange)
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
    
    return(grange)
}
# }}}

get_annotation(args$assembly)
x <- list('genes'=genes, 'promoters'=tss_window)

if ( file_ext(args$bed) %in% c('broadPeak', 'narrowPeak')) {
    grange <- parse_broadPeak(args$bed)
    export_gtfs(args$bed, grange)
    mclapply(names(x), function(i) export_bound(x[[i]], i, grange, args$bed))
} else if(file_ext(args$bed) %in% c('bed')){
          grange <- parse_bed5(args$bed)
          export_gtfs(args$bed, grange)
          
          current <- gsub(file_ext(args$bed), 'bound.txt',  file)
          bound <- mclapply(names(x), function(i) export_bound(x[[i]], i, grange, args$bed))
          names(bound) <- names(x)
          df <- as.data.frame(bound)
}else {
    warning("Can not convert bed file to gtf [unknown format].")
}


