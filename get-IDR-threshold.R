#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)
library(dplyr)

parser =  ArgumentParser(description="Get IDR thresholds")
parser$add_argument('-p', '--pre_IDR', metavar= "file", required='True', type= "character",
                    help= "TSV with number of peaks before IDR generated with get-IDR-threshold PHONY in Makefile")
parser$add_argument('-o', '--output', metavar= "file", required='True', type= "character",
                    help= "Output filename")
parser$add_argument('-t', '--type', action ="store_true", default = FALSE, help = "Set to true for returning thresholds for pooled consistency. Default = FALSE")

args = parser$parse_args()
#}}}

options(scipen = 999)

if(args$type){
    print('Reporting thresholds for pooled consistency')
    low = 0.005
    high = 0.01
} else {
    print('Reporting thresholds for self consistency')
    low = 0.02
    high = 0.05
}

df = read.delim(args$pre_IDR, head=F, stringsAsFactors=F)
df = df %>% select(V2,V3) %>%
    mutate(group = gsub("_[^-]*_peaks.*", "", V3)) %>%
    filter(group != "total") %>%
    group_by(group) %>% summarise(min = min(V2), max=max(V2)) %>% 
    mutate(top = ifelse(min < 100000, 100000, 150000),
           thr = ifelse(min < 100000, high, low), 
           tmp = thr) %>%
    tidyr::unite(output, group, tmp, sep="_conservative_") %>%
    mutate(output = basename(gsub("$", ".narrowPeak", output)))
    
write.table(df, args$output, col.names=T, row.names=F, quote=F, sep="\t")
