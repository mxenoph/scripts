#!/usr/bin/env Rscript

library(argparse)
library(tools)

parser =  ArgumentParser(description="Explore peak co-localization by correlation")
parser$add_argument('-c', '--config', metavar= "file", required='True', type= "character", help= "TSV file with columns Condition (e.g. ESC-2i), factor (e.g. H3AC), peaks bed file")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
parser$add_argument('-l', '--label', type= "character", default=format(Sys.time(), "%d%b%Y"), help= "Give a label for the output")

args = parser$parse_args()

output_path = file.path(args$out)
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

if(FALSE){
    args = list()
    args$config = '/nfs/research2/bertone/user/mxenoph/hendrich/thesis/corr-plot-esc.config'
    args$out = '/nfs/research2/bertone/user/mxenoph/hendrich/thesis'
    args$label = 'esc-matrix'
}

# Import libraries & Turn off warning messages for loading the packages-globally# {{{
suppress = base::suppressPackageStartupMessages
options(warn=-1)
x = c('GenomicRanges', 'rtracklayer', 'plyr', 'tools', 'corrplot')
lapply(x, suppressMessages(library), character.only=T)
source("~/source/Rscripts/ggplot-functions.R")
source("~/source/Rscripts/granges-functions.R")
#Turn warnings back on
options(warn=0)

# Local functions# {{{
cor_mtest = function(mat, conf_level = 0.95) {# {{{
    # Initiate 3 matrices of same dimensions as the input mat.
    mat = as.matrix(mat)
    n = ncol(mat)
    p_mat = low_CI_mat = upp_CI_mat = matrix(NA, n, n)
    # Set the diagonal of the p values matrix to 0
    diag(p_mat) = 0
    diag(low_CI_mat) = diag(upp_CI_mat) = 1
    
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp = cor.test(mat[, i], mat[, j], conf.level = conf_level)
            p_mat[i, j] = p_mat[j, i] = tmp$p.value
            low_CI_mat[i, j] = low_CI_mat[j, i] = tmp$conf.int[1]
            upp_CI_mat[i, j] = upp_CI_mat[j, i] = tmp$conf.int[2]
        }
    }
    return(list(p_mat, low_CI_mat, upp_CI_mat))
}
# }}}

subset_by_cor = function(mat, corr_cutoff=0.6){# {{{
    m = abs(mat) >= corr_cutoff
    sets = c()

    for(i in 1:nrow(m)){
       for(j in i:ncol(m)){
              if(rownames(m)[i] == colnames(m)[j]){
                  next
              }else if(m[i,j]) {
                  sets = c(sets, paste(rownames(m)[i],colnames(m)[j], sep="-"))
              }
       }
    }
    return(sets)
}# }}}

# peak_universe= all_peaks, peaks = peaks, sets= subset_by_cor output
get_common_peaks = function(peak_universe, sets, peaks, output_path){# {{{
    x = unlist(strsplit(sets, '-'))
    if(length(x) != 2 )stop('Comparison can only be made for 2 sets')
    out = mclapply(peaks[x], function(p) findOverlaps(peak_universe, p))
    subset = peak_universe[queryHits(out[[1]])[queryHits(unique(out[[1]])) %in% queryHits(unique(out[[2]]))]]
    file = paste0(gsub(':', '_', sets), '.bed')
    export.bed(subset, file.path(output_path, file))
}

# }}}

# }}}

# Prepare the data -- make binary matrix of binding sites in each chip# {{{
conf = read.table(args$config, header=T, sep="\t")

# Running it interactively when conf is set manually requires as.character
# otherwise it gives you factors
#factors = as.character(conf[['factor']])
factors = as.character(apply( conf[ , c('Condition', 'factor') ] , 1 , paste , collapse = ":" ))
factors = toupper(factors)

peaks = GRangesList() # {{{
for (indx in 1:length(conf[['file']])) {
#    if (grepl('peaks', as.character(conf[['file']])[indx])){
#        print(as.character(conf[['file']])[indx])
#        ranges = import.bed(as.character(conf[['file']])[indx])
#    } else{
    # skipping ucsc track line
#    if (grepl('track', readLines(as.character(conf[['file']])[indx], n=1))){
#        tmp = read.table(as.character(conf[['file']])[indx], header=F, sep="\t", skip=1)
#    } else{
#        tmp = read.table(as.character(conf[['file']])[indx], header=F, sep="\t")
#    }
#        ranges = GRanges(seqnames=tmp[,1], IRanges(start= tmp[,2], end= tmp[,3]), strand="*")
#        head(ranges)
#    }
     if(grepl('.bed', as.character(as.character(conf[['file']])[indx]))){
         ranges = rtracklayer::import.bed(as.character(conf[['file']])[indx])
     } else{
         if(file_ext(as.character(conf[['file']])[indx]) == 'gappedPeak') {
             ranges = match.fun(paste('import', file_ext(as.character(conf[['file']])[indx]), sep = '_'))(as.character(conf[['file']])[indx])
             ranges = ranges$blocks[width(ranges$blocks) > 10]
         } else {
             ranges = match.fun(paste('import', file_ext(as.character(conf[['file']])[indx]), sep = '_'))(as.character(conf[['file']])[indx])
         }
     }
    
    # Getting rid of tags in contigs if present
    ranges = ranges[seqnames(ranges) %in% seqlevels(ranges)[grep('chr', seqlevels(ranges))]]

    if (duplicated(factors)[indx]) {
        # merging multiple data sets for the same factor
        peaks[[factors[indx]]] = c(peaks[[factors[indx]]], reduce(ranges))
    } else {
          peaks = c(peaks, GRangesList(reduce(ranges)))
          names(peaks)[length(peaks)] = factors[indx]
    }
}

# since for some factors several datasets got combined use reduce to remove overlapping peaks
peaks = endoapply(peaks, reduce)

# }}}

all_peaks = reduce(unlist(peaks))
#all_peaks = unlist(peaks)

# binary table with all binding sites for the TFs in the config file
all_peaks_binary = rep(0, length(all_peaks) * length(names(peaks)))
dim(all_peaks_binary) = c(length(all_peaks), length(names(peaks)))
colnames(all_peaks_binary) = names(peaks)

for (indx in 1:length(names(peaks))) {
    match = as.matrix(findOverlaps(all_peaks, peaks[[indx]]))
    all_peaks_binary[unique(match[,1]), indx] = rep(1, length(unique(match[,1])))
}
# }}}

# Set correlation matrix colours and plotting parameters# {{{
col = rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7","#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
#op = par(family="serif")# }}}

# Calculate correlation# {{{
M = cor(all_peaks_binary)
# Get hclust order
order.hc = corrMatOrder(M, order="hclust")
# Reorder correlation matrix
M.hc  = M[order.hc, order.hc ]

# Test significance of correlation
conf_level = 0.95
p_values = cor_mtest(all_peaks_binary, conf_level)
p_cutoff = 0.01
# }}}

# Colour the labels according to condition and format labels# {{{
x = factor(gsub(":.*", "", colnames(M.hc)))

tlcol = as.character(mapvalues(x,
                                from= levels(x),
#                                to= gg_color_hue(length(levels(x)))))
                                to= c("#336699", "#663399")))
#label = gg_color_hue(length(levels(x)))
label = c("#336699", "#663399")
names(label) = levels(x)

colnames(M.hc) = gsub(".*:", "", colnames(M.hc))
rownames(M.hc) = gsub(".*:", "", rownames(M.hc))# }}}

# Print plots# {{{
pdf.options(encoding = 'ISOLatin2.enc')
pdf(file.path(plot_path, paste0(args$label, "-CorrelationPlot.pdf")), paper='a4')

    corrplot(M, method="color", type="lower",# {{{
             order="hclust", tl.col="black", tl.cex=0.5,
             addrect=4, col=colorRampPalette(col)(200))# }}}

    corrplot(M.hc, method="color", type = "lower",# {{{
             order="original", tl.col=tlcol, tl.cex=0.5,
             addrect=5, col=colorRampPalette(col)(200))
    legend('topright', names(label), text.col=as.character(label), bty='n')# }}}

    colnames(M) = gsub(".*:", "", colnames(M))# {{{
    rownames(M) = gsub(".*:", "", rownames(M))
    corrplot(M, method="color",
             order="hclust", tl.col=tlcol, tl.cex=0.5,
             addrect=6, col=colorRampPalette(col)(200))
    legend('topright', names(label), text.col=as.character(label), bty='n')
    # }}}

    corrplot(M, p.mat= p_values[[1]], sig.level=p_cutoff, method = "color", # {{{
#             order="hclust", tl.cex=1, tl.col = "black", tl.srt=90, number.digits = 2,
             order="hclust", tl.col=tlcol, tl.cex=1, tl.srt=90, number.digits = 2,
             number.cex = 0.5, cl.cex = 1,
             addCoef.col="grey", addrect=6, col=colorRampPalette(col)(200))
    
    legend('topright', names(label), text.col=as.character(label), bty='n')
    # }}}

    # Combine with a significant test. Provide the matrix of correlation but also a matrix of p-val for all correlations computed# {{{# {{{
    # sig.level is the top value a p-val is considered significant
    corrplot(M, p.mat= p_values[[1]], sig.level=p_cutoff,
             order="hclust", method="color", #type="lower",
             tl.col=tlcol, tl.cex=0.5,
             insig="pch", pch.cex=0.7,
             addrect=6, col=colorRampPalette(col)(200))
    legend('topright', names(label), text.col=as.character(label), bty='n')
    #legend('topright', c(names(label), paste0('Confidence level: ', conf_level, ", p-val < ", p_cutoff)),
    #       text.col=c(as.character(label), 'black'), bty='n')
    # }}}# }}}

    #corrplot(M, method="color", type="lower",# {{{
    #         order="FPC", tl.col="black", tl.cex=0.5, 
    #         col=colorRampPalette(col)(200))
    #legend("topright", 'Ordered based on First Principal Component', col="black", bty='n')
    ## }}}

dev.off()
# }}}


