#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")

parser <-  ArgumentParser(description="Compare macs to macs2")
parser$add_argument('-x', '--peaks', metavar= "file", required='True', type= "character", help= "MACS generated excel file for called peaks")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args <- parser$parse_args()

plot_path <- file.path(args$out, 'plots')
dir.create(plot_path, recursive= TRUE)

# function in annotation-functions.R
get_annotation(tolower(args$assembly))

# }}}

# Import libraries & Turn off warning messages for loading the packages-globally# {{{
suppress <- base::suppressPackageStartupMessages
options(warn=-1)
source("~/source/Rscripts/granges-functions.R")
#Turn warnings back on
options(warn=0)
# }}}

# Functions # {{{
#}}}

args <- list()
args$out <- "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/macs"
setwd(args$out)
args$peaks <- list.files(pattern="macs2.*_peaks.xls")
args$assembly <- 'mm10'

macs2_peaks <- lapply(args$peaks, macs2_to_granges)
names(macs2_peaks) <- args$peaks

macs_summits_f <- list.files(pattern="_summits.bed")
macs_summits <- lapply(macs_summits_f, import)
names(summits) <- macs_summits_f

sharp_peaks <- macs2_peaks[grep("summits_peaks", names(macs2_peaks), perl=T)]
broad_peaks <- macs2_peaks[grep("summits_peaks", names(macs2_peaks), perl=T, invert=T)]

sharp_peaks <- lapply(sharp_peaks, function(x){ names(x) <- paste0('MACS_peak_', 1:length(x)); return(x)})
broad_peaks <- lapply(broad_peaks, function(x){ names(x) <- paste0('MACS_peak_', 1:length(x)); return(x)})

comparison <- annotate_range(broad_peaks[[1]], broad_peaks[[2]])
xy_intersection <- intersect(broad_peaks[[1]][1:100], broad_peaks[[2]][1:100])
names(xy_intersection) <- paste0('intersection_', 1:length(xy_intersection))


xy_union <- union(broad_peaks[[1]][1:100], broad_peaks[[2]][1:100])
names(xy_union) <- paste0('union_', 1:length(xy_union))
x <- broad_peaks[[1]][1:100]
y <- broad_peaks[[2]][1:100]
x <- annotate_range(x, xy_union)
y <- annotate_range(y, xy_union)

x <- GRanges(rep('chr1', 5), IRanges(c(1, 8, 15, 30, 45), c(5, 10, 20, 40, 50)), strand="*")
y <- GRanges(rep('chr1', 3), IRanges(c(1, 15, 100), c(9, 20, 200)), strand="*")
ia <- length(ov)- length(unique(queryHits(ov)))
ib <- length(ov)- length(unique(subjectHits(ov)))

files <- list.files(pattern="Epi.*bed", args$out)
dictionary <- list('7E12' = 'wt', '7G9' = 'ko', '_peaks.bed' = '')
names(files) <- gsub(names(dictionary)[1], dictionary[[1]], files)
names(files) <- gsub(names(dictionary)[2], dictionary[[2]], names(files))
names(files) <- gsub(names(dictionary)[3], dictionary[[3]], names(files))
bed <- mclapply(files, import.bed)
names(bed) <- names(files)

compare_mbd3_chd4 <- function(files){
    
}
