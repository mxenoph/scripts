#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
suppressMessages(library(argparse))
suppressMessages(library(tools))
source("~/source/Rscripts/annotation-functions.R")

parser <-  ArgumentParser(description="Define enhancer elements")
parser$add_argument('-c', '--config', metavar= "file", required='True', type= "character", help= "config file")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args <- parser$parse_args()
output_path <- file.path(args$out, 'enhancers')
plot_path <- file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

#}}}

# Load packages# {{{
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(biovizBase))
suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer))# }}}

#Functions
#Caution: different from previous script: take intersection for data coming from different scripts
parse_peaks <- function(config){# {{{
    source("~/source/Rscripts/granges-functions.R")
    # This assumes it's mouse data. If not this has to be changed
    chr <- c(paste0('chr', 1:19), 'chrX', 'chrY')
    
    y <- mclapply(as.vector(config[,'file']), function(files){
                  bed <- macs_to_granges(files)
                  bed <- subset(bed, seqnames(bed) %in% chr)
                  names(bed) <- name_gr(bed)
                  return(bed)
    })
    names(y) <- with(config, paste(system, factor, rownames(config), sep=':'))
    
    descriptors <- unique(gsub(":[0-9]{1,2}$", '', names(y)))
    pooled <- sapply(descriptors, function(x){
        i <- grep(x, names(y))
        if(length(i) == 1) {
            return(y[[i]])
        }else {
            gr <- intersect_multi(y[i])
            for(j in colnames(values(y[[1]])) ) values(gr)[j] <- rep(NA, length(gr))

            scores <- mclapply(y[i], function(x) get_hits_scores(gr, x, 'FE'))
            names(scores) <- names(y[i])
            scores <- do.call(cbind.data.frame, scores)
            scores[,'mean'] <- apply(scores,1,mean)
            values(gr)['FE'] <- scores[['mean']]
            #apply(scores,1,sd)
            return(gr)
        }
    })
    return(list(y, pooled))
}# }}}

#config <- read.table(args$config, sep=",", header=TRUE, stringsAsFactors=TRUE)
config <- read.table(args$config, sep="\t", header=TRUE, stringsAsFactors=TRUE)

# list of granges grL[[1]] it's all the data while grL[[2]] it's pooled experiments for same mark
grL <- parse_peaks(config)

active <- as.character(unique(
                              config[with(config, element=='enhancer' & type=='active'),
                                     'factor']))
poised <- as.character(unique(
                              config[with(config, element=='enhancer' & type=='poised'),
                                     'factor']))
promoter <- as.character(unique(
                              config[with(config, element=='promoter' & type=='active'),
                                     'factor']))

promoter <- mclapply(promoter, function(x){
                     intersect_multi(grL[[2]][unlist(mclapply(x, grep, names(grL[[2]])))])
    })

enhancers <- mclapply(list(active, poised), function(x){
                      intersect_multi(grL[[2]][unlist(mclapply(x, grep, names(grL[[2]])))])
    })

enhancers <- mclapply(enhancers, function(gr){
                      ov <- findOverlaps(gr, promoter[[1]])
                      # Removing elements that overlap with promoters by 1 bp
                      keep <- (1:length(gr))[! 1:length(gr) %in% unique(queryHits(ov))]
                      reduce(gr[keep])
   })
names(enhancers) <- c('active','poised')

enhancers[[1]][seqnames(enhancers[[1]]) == 'chr17' & ranges(enhancers[[1]]) %over% IRanges(start=35636289, end=35649778)]
promoter[[1]][seqnames(promoter[[1]]) == 'chr17' & ranges(promoter[[1]]) %over% IRanges(start=35636289, end=35649778)]

mclapply(names(enhancers), function(x) export.bed(enhancers[[x]], paste0(x, '-enhancers.bed')))

