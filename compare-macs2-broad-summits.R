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
args$peaks <- list.files(pattern="macs2.*_peaks.xls")
args$assembly <- 'mm10'

gr <- lapply(args$peaks, macs2_to_granges)
names(gr) <- args$peaks

list.files(pattern="_summits.bed")
summits_files <- list.files(pattern="_summits.bed")
summits <- lapply(summits_files, import)
names(summits) <- summits_files
test[grep("_peak_\\d+$", test$name, perl=T, invert=T)]
