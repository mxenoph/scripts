#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Contrast mean signal at peaks for ")
parser$add_argument('-f', '--files', metavar= "file", required='True', nargs = "+", type= "character",
                    help= "Output files of bigWigAverageOverBed.")
parser$add_argument('-l', '--label', type= "character",
                    default = format(Sys.time(), "%d%b%Y"), help= "Give a label for the output")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(),
                    help= "Output directory -- all subdirectories will be created here")

args = parser$parse_args()
plot_path = file.path(args$out, 'plots')
dir.create(plot_path, recursive= TRUE)
# }}}

# Libraries
library(dplyr)
library(ggplot2)
library(reshape2)

# Functions
add_rownames = function(df, var = 'rowname') {
    stopifnot(is.data.frame(df))
    rowname_df = setNames(data_frame(rownames(df)), var)
    cbind(rowname_df, df)
}

format_plot_data = function(x, y, labels) {
    if(nrow(x) != nrow(y) | !all(x[['name']] %in% y[['name']])){
        warning('Files do not have matching region names.')
    } else {
        # ensure order is the same by merging the dataframes
        plot_data = x %>% select(name, mean) %>%
                    rename_(.dots = setNames(names(.), gsub("mean", labels[1], names(.)))) %>%
                    left_join((y %>% select(name, mean) %>%
                               rename_(.dots = setNames(names(.), gsub("mean", labels[2], names(.))))),
                              by = 'name')
        return(plot_data)
    }
}

stats = lapply(args$files, function(x) {
                   if (file.exists(x)){
                       tmp = read.delim(x, sep = "\t", header = FALSE,
                                        col.names = c('name', 'size', 'covered', 'sum', 'mean0', 'mean'))
                       return(tmp)
                   } else {
                       command = "bigWigAverageOverBed in.bw out.tsv"
                       warning(paste('do not exist. Run `', command, '` first.', collapse = ' '))
                   }
                    })

names(stats) = gsub('\\.bigWigAverageOverBed', '', file_path_sans_ext(basename(args$file)))

pairwise_combn = expand.grid(combn(names(stats),2))
plot_data = lapply(pairwise_combn, function(x){
                       tmp_list = stats[levels(x)]
                       format_plot_data(tmp_list[[1]], tmp_list[[2]], levels(x))
                    })

pdf(file.path(plot_path, paste0(args$label, 'mean-signal.pdf')))
for (x in names(plot_data)) {
    tmp = plot_data[[x]] %>% mutate(group = as.factor(gsub("_[0-9]*$", '', name)))
    n = tmp %>% select(-name) %>% names()
    p = ggplot(tmp, aes_string(n[1], n[2])) + geom_point(alpha = 0.5, aes_string(colour = n[3]))
    p = p + theme_classic() + theme(legend.position="bottom")
    p = p + geom_smooth(method = "lm", se=F)
    p
}
dev.off()



