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
library(tidyr)

for(f in args$file){
    bowtie = read.table(f, sep = "\t", header = TRUE)
    # dup_removed are all reads mapped_duplicate_free + unmapped reads (fail or more_mult)
    # so to get the uniquely mapped duplicate free read count one should subtract the failed alignments from 
    # dup_removed
    bowtie = bowtie %>% mutate(unique = (dup_removed - fail - more_mult),
                               mapped = mapped,
                               duplicated = (mapped - unique),
                               fail = fail,
                               more_mult = more_mult
                               )

    # Saving calculations on stats (uniq without ddup counts included) and whether samples need # {{{
    # downsampling. 1 means no and any other number means downsample to that percentage to
    # establish same depth with other reps in the group condition-protein
    tmp = bowtie %>% mutate(tmp = file_path_sans_ext(filename)) %>%
                separate(tmp, into=c('condition','protein'), sep ='-') %>%
                mutate(condition = gsub('Remco', '', condition), protein = gsub("_.*", '', protein)) %>% 
                unite(tmp, condition, protein, sep='-') 
                
    summarised = tmp %>% group_by(tmp) %>%
                    summarise(group_median = median(unique))

    tmp = tmp %>% left_join(summarised, by = 'tmp') %>%
        mutate(index = group_median / unique,
               downsample = ifelse(round(index, digits = 2) < 0.7, round(index, digits = 1),
                                   ifelse(round(index, digits = 2) > 1.3, 'seq. depth < group median', 1))) %>%
        select(filename, fail, more_mult, unique, duplicated, all, group_median, index, downsample)

    write.table(tmp,
                file.path(output_path, 'alignment-summarised-stats.tsv'), sep="\t", row.names = F, quote = F)
    # }}}

    plot_data = melt((bowtie %>% select(filename, fail, more_mult, unique, duplicated, all)), id = 'filename')
    plot_data = plot_data %>% tidyr::separate(filename, into = c("condition", "protein"), sep="-" ) %>%
                mutate(protein = gsub(".stderr", "", protein),
                       ID = paste(condition, protein, sep="-"))
    plot_data[['value']] = as.numeric(plot_data[['value']])
    plot_data[['ID']] = as.factor(plot_data[['ID']])
    plot_data[['condition']] = as.factor(plot_data[['condition']])

    cols = c('all' = "#B2506F",
             'unique' = "#B05BB3",
             'duplicated' = "#667EA0",
             'more_mult' = "#B6653F",
             'fail' = "#4A8F75")
    plot_data[['variable']] = factor(plot_data[['variable']], levels = names(cols))

    pdf(file.path(plot_path, 'bowtie-statistics.pdf'), page = 'a4')
    for (cond in levels(plot_data$condition)){
        q = ggplot(aes(x =ID, y = value / 1000000, fill = variable), data = subset(plot_data, condition == cond))
        q = q + geom_bar(stat = "identity", width = .3, position = 'dodge')
        q = q + theme_bw()
        q = q + theme(axis.line = element_line(colour = "black"),
                      panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(), 
                      axis.text.x = element_text(angle = 25, hjust = 1))
        q = q + scale_fill_manual(values = cols, breaks = names(cols), name = "Alignment types")
        plot(q)
    }

    bowtie = bowtie %>% mutate(unique = unique * 100 / all,
                               mapped = mapped * 100 / all,
                               duplicated = duplicated * 100 /all,
                               fail = fail * 100 / all,
                               more_mult = more_mult * 100 / all
                               )
    bowtie = bowtie %>% select(filename, unique, duplicated, more_mult, fail)

    plot_data = melt(bowtie, id = 'filename')
    plot_data = plot_data %>% tidyr::separate(filename, into = c("condition", "protein"), sep="-" ) %>%
                mutate(protein = gsub(".stderr", "", protein),
                       ID = paste(condition, protein, sep="-"))
    plot_data[['value']] = as.numeric(plot_data[['value']])
    plot_data[['ID']] = as.factor(plot_data[['ID']])
    plot_data[['condition']] = as.factor(plot_data[['condition']])
    plot_data[['variable']] = factor(plot_data[['variable']], levels = names(cols))

    q = ggplot(aes(x =ID, y = value, fill = variable, order = variable), data = plot_data)
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
    q = q + scale_fill_manual(values = cols, breaks = names(cols), name = "Alignment types")
    plot(q)

    dev.off()
}
