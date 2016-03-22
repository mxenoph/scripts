#!/usr/bin/env Rscript

# Parse arguments # {{{
library = function (...) suppressMessages(base::library(...))
options(stringsAsFactors=F)
library(argparse)
library(tools)

parser =  ArgumentParser(description="Takes as input the file produced from check-chip-enrichment.pl and plots the cumulative distribution of tags, like CHANCE")
parser$add_argument('-b', '--bins', metavar= "file", required='True', type= "character", nargs="*", help= "output of check-chip-enrichment.pl, always chip should be given first")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args = parser$parse_args()

plot_path = file.path(args$out, 'plots')
dir.create(plot_path, recursive= TRUE)

library(dplyr)
library(ggplot2)
library(gridExtra)
# needed for the use of percent in x axis labels
library(scales)
# }}}

enrichment = do.call(rbind,
                     lapply(args$bins, function(bins){
                            reads = read.delim(bins, header=F)
                            reads = reads[order(reads[,2]), ]
                            reads =  cbind(bin = 1:dim(reads)[1], cum_reads = (cumsum(reads[, 2]) / sum(reads[, 2])))
                            reads = as.data.frame(reads) %>% mutate(ID = basename(bins))
                            return(reads)
                                }))

png(file.path(plot_path, paste0(basename(file_path_sans_ext(args$bin[1])), '.png')))

# Check if all files were generate for the same window so you can plot together
validity = sapply(basename(args$bins), function(x) nrow(enrichment[enrichment$ID==x,]) )

if (all(validity == validity[1])){
    breaks = (validity[1]*20)/100
    breaks = seq(from = breaks, to = validity[1],by = breaks)
    print(breaks)
    print(seq(from=0.2,to=1,by=0.2))

    p = ggplot(enrichment, aes(x=bin, y=cum_reads, colour=ID))
    p = p + geom_line(size=1) + scale_x_continuous(labels = seq(from=0.2,to=1,by=0.2), breaks = breaks)
    p = p + labs(x="Percentage of bins", y="Cumulative percentage of mapped reads") + theme_classic()
    p + theme(legend.justification=c(0,1), legend.position=c(0,1))
} else {
    p = ggplot(enrichment, aes(x=bin, y=cum_reads))
    p = p + geom_line(size=1) + scale_x_continuous(labels=percent) + facet_grid(. ~ ID)
    p = p + labs(x="bin", y="Cumulative percentage of mapped reads") + theme_classic()
    p
}
dev.off()
