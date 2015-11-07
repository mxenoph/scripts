#! /usr/bin/env bash

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)

parser =  ArgumentParser(description="Plot bowtie alignment statistics")
parser$add_argument('-f', '--file', required = 'True', nargs = '+',  metavar = "file", type= "character",  help= "Output file of count-bowtie-reads.sh")
parser$add_argument('-o', '--out', required = 'True', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
args = parser$parse_args()

output_path = file.path(args$out)
if(! grepl("plots", output_path)){
    plot_path = file.path(output_path, 'plots')
} else {
    plot_path = output_path
}
dir.create(plot_path, recursive= TRUE)# }}}

options(stringsAsFactors=F)

library(ggplot2)
library(reshape)
library(dplyr)

for(f in args$file){
    bowtie = read.table(f, sep = "\t", header = TRUE)
    print(head(bowtie))
    # dup_removed are all reads mapped_duplicate_free + unmapped reads (fail or more_mult)
    # so to get the uniquely mapped duplicate free read count one should subtract the failed alignments from 
    # dup_removed
    bowtie = bowtie %>% mutate(unique = (dup_removed - fail - more_mult) * 100 / all,
                               mapped = mapped * 100 / all,
                               duplicated = (mapped - unique),
                               fail = fail * 100 / all,
                               more_mult = more_mult * 100 / all
                               )
    bowtie = bowtie %>% select(filename, fail, more_mult, unique, duplicated)

    plot_data = melt(bowtie, id = 'filename')
    plot_data[['value']] = as.numeric(plot_data[['value']])
    
    cols = c('fail' = "#4A8F75", 'more_mult' = "#B6653F", 'unique' = "#B05BB3", 'duplicated' = "#667EA0")

    pdf(file.path(plot_path, 'bowtie-statistics.pdf'))
    q = ggplot(aes(x =filename, y = value, fill = variable), data = plot_data)
    q = q + geom_bar(stat = "identity", width = .3)
    q = q + theme_bw()
    q = q + coord_flip()
    q = q + theme(axis.line = element_line(colour = "black"),
                  panel.grid.major.x = element_line(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank())
    q = q + scale_y_continuous(limits=c(0, 100), name = "Percentage of total reads")
    q = q + scale_fill_manual(values = cols, name = "Alignment types")
    plot(q)

    dev.off()
}
