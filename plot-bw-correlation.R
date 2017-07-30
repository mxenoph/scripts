#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Get exons, introns, intergenic for given annotation")
parser$add_argument('-s', '--sample', metavar= "file", required='True', type= "character",
                    help= ".sample-correlation file generated in ChIP Makefile with bigWigCorrelate")

args = parser$parse_args()
# }}}
library(dplyr)
library(corrplot)

lines = readLines(args$sample)
files = lines[seq(1,length(lines),2)]
files = as.data.frame(files) %>% tidyr::separate(files, into = c('file_1', 'file_2'), sep=" ") %>%
        dplyr::mutate(file_1 = basename(file_1),
                      file_2 = basename(file_2))
files$corr = lines[seq(2,length(lines),2)]

len = length(unique(c(files[[1]], files[[2]])))
mat = matrix(nrow = len, ncol = len)
colnames(mat) = rownames(mat) = unique(c(files[[1]], files[[2]]))
head(mat)

for(i in colnames(mat)){
    for(j in rownames(mat)){
        if(i==j){
            tmp = 1
        } else{
            tmp = files %>% filter((file_1 == i & file_2 == j) | (file_1 == j & file_2 == i)) %>% 
                    dplyr::select(corr) %>% .[[1]] %>% as.numeric()
            if(length(tmp) != 1){
                print('Comparison found more than one time. Ensuring value is the same.')
                if(length(unique(tmp)) !=1) stop('Same comparisons do not give the same value.')
                tmp = unique(tmp)
            }
        }
        mat[j,i] = tmp
    }
}
col = rev(c("#67001F", "#B2182B", "#D6604D",
            "#F4A582", "#FDDBC7","#FFFFFF",
            "#D1E5F0", "#92C5DE", "#4393C3",
            "#2166AC", "#053061"))

pdf(gsub(".sample-correlation", ".pdf", args$sample))
corrplot(mat, method="color", type="lower",
         order="hclust", tl.col="black", tl.cex=1,
         tl.srt=90, number.digits = 2, number.cex = 0.5, cl.cex = 1,
         addCoef.col="grey",
         addrect=4, col=colorRampPalette(col)(200))
dev.off()
