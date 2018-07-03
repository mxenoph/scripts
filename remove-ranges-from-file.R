#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Define enhancer elements")
parser$add_argument('-p', '--priority', metavar= "file", required='True', type= "character", help= "BED file in which to keep ranges common to both files")
parser$add_argument('-s', '--secondary', metavar= "file", required='True', type= "character", help= "BED file in which to remove ranges common to both files")
parser$add_argument('-o', '--outfile', metavar= "path", type= "character", help= "Output file to save the results from cleaning up the secondary file")

args = parser$parse_args()

#}}}

# script used to clean up poised enhancer list that contains active enhancers, if made with an older version of the script
library(rtracklayer)
library(GenomicRanges)

priority = import.bed(args$priority)
names(priority) = paste(seqnames(priority), start(priority), end(priority), sep =':')

secondary = import.bed(args$secondary)
names(secondary) = paste(seqnames(secondary), start(secondary), end(secondary), sep =':')
# remove active enhancers from poised ones

ov = findOverlaps(secondary, priority)
secondary = secondary[!names(secondary) %in% names(secondary[unique(queryHits(ov))])]
export.bed(secondary, args$outfile)

