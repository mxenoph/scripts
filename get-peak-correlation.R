 #!/usr/bin/env Rscript --slave

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)

parser <-  ArgumentParser(description="Explore peak co-localization by correlation")
parser$add_argument('-c', '--config', metavar= "file", required='True', type= "character", help= "TSV file with columns Condition (e.g. ESC-2i), factor (e.g. H3AC), peaks bed file")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
parser$add_argument('-l', '--label', type= "character", default=format(Sys.time(), "%d%b%Y"), help= "Give a label for the output")

args <- parser$parse_args()
# }}}}

x <- c('GenomicRanges', 'rtracklayer', 'plyr', 'tools', 'corrplot')
lapply(x, suppressMessages(library), character.only=T)
source("~/source/Rscripts/ggplot-functions.R")

# Local functions# {{{
cor.mtest <- function(mat, conf.level = 0.95) {# {{{
    # Initiate 3 matrices of same dimensions as the input mat.
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
    # Set the diagonal of the p values matrix to 0
    diag(p.mat) <- 0
    diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
            lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
            uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
        }
    }
    return(list(p.mat, lowCI.mat, uppCI.mat))
}
# }}}

subset_by_cor <- function(mat, corr_cutoff=0.6){# {{{
    m <- abs(mat) >= corr_cutoff
    sets <- c()

    for(i in 1:nrow(m)){
       for(j in i:ncol(m)){
              if(rownames(m)[i] == colnames(m)[j]){
                  next
              }else if(m[i,j]) {
                  sets <- c(sets, paste(rownames(m)[i],colnames(m)[j], sep="-"))
              }
       }
    }
    return(sets)
}# }}}

# }}}

# read files and combine binding profiles of factors with several data sets
#args <- commandArgs(trailingOnly= TRUE)

# conf.txt colnames: factor\tfile
#conf <- read.table(args[1], header=T, sep="\t")
#name <- args[2]
#out_dir <- args[3]
#plot_dir <- file.path(out_dir, 'plots', '')
#

# This assumes it's mouse data. If not this has to be changed
#chr <- c(paste0('chr', 1:19), 'chrX', 'chrY')

conf <- read.table(args$config, header=T, sep="\t")

# Running it interactively when conf is set manually requires as.character
# otherwise it gives you factors
#factors <- as.character(conf[['factor']])
factors <- as.character(apply( conf[ , c('Condition', 'factor') ] , 1 , paste , collapse = ":" ))
factors <- toupper(factors)

grL <- GRangesList()# {{{
for (indx in 1:length(conf[['file']])) {
#    if (grepl('peaks', as.character(conf[['file']])[indx])){
#        print(as.character(conf[['file']])[indx])
#        ranges <- import.bed(as.character(conf[['file']])[indx])
#    } else{
        tmp <- read.table(as.character(conf[['file']])[indx], header=F, sep="\t")
        ranges <- GRanges(seqnames=tmp[,1], IRanges(start= tmp[,2], end= tmp[,3]), strand="*")
#        head(ranges)
#    }
    # Getting rid of tags in contigs if present
    ranges <- ranges[seqnames(ranges) %in% seqlevels(ranges)[grep('chr', seqlevels(ranges))]]

    if (duplicated(factors)[indx]) {
        # merging multiple data sets for the same factor
        grL[[factors[indx]]] <- c(grL[[factors[indx]]], ranges)
      } else {
          grL <- c(grL, GRangesList(ranges))
          names(grL)[length(grL)] <- factors[indx]
    }
}# }}}
# since for some factors several datasets got combined use reduce to remove overlapping peaks
grL = endoapply(grL, reduce)

#all_ranges <- reduce(unlist(grL))
all_ranges <- unlist(grL)
#print('length', length(all_ranges))

# binary table with all binding sites for the TFs in the config file
tab_all <- rep(0, length(all_ranges) * length(names(grL)))
dim(tab_all) <- c(length(all_ranges), length(names(grL)))
colnames(tab_all) <- names(grL)

for (indx in 1:length(names(grL))) {
  match = as.matrix(findOverlaps(all_ranges, grL[[indx]]))
  tab_all[unique(match[,1]), indx] = rep(1, length(unique(match[,1])))
}


# Set correlation matrix colours
col <- rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7","#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))

# Calculate correlation
M <- cor(tab_all)
# Get hclust order
order.hc <- corrMatOrder(M, order="hclust")
# Reorder correlation matrix
M.hc  <- M[order.hc, order.hc ]

# Colour the labels according to condition -- this needs gg_color_hue function
x <- factor(gsub(":.*", "", colnames(M.hc)))

tlcol <- as.character(mapvalues(x,
                                from= levels(x),
                                to= gg_color_hue(length(levels(x)))))
label <- gg_color_hue(length(levels(x)))
names(label) <- levels(x)

colnames(M.hc) <- gsub(".*:", "", colnames(M.hc))
rownames(M.hc) <- gsub(".*:", "", rownames(M.hc))

# Test significance of correlation
conf.lev <- 0.95
sig <- cor.mtest(tab_all, conf.lev)
sig.lev <- 0.01

op <- par(family="serif")

# Print plots# {{{
pdf(paste0(plot_dir, name, "-CorrelationPlot.pdf"), paper='a4')
corrplot(M, method="color", type="lower",
         order="hclust", tl.col="black", tl.cex=0.5,
         addrect=4, col=colorRampPalette(col)(200))

corrplot(M.hc, method="color", type="lower",
         order="original", tl.col=tlcol, tl.cex=0.5,
         addrect=5, col=colorRampPalette(col)(200))

legend('topright', names(label), text.col=as.character(label), bty='n')

colnames(M) <- gsub(".*:", "", colnames(M))
rownames(M) <- gsub(".*:", "", rownames(M))
corrplot(M, method="color",
         order="hclust", tl.col=tlcol, tl.cex=0.5,
         addrect=9, col=colorRampPalette(col)(200))

legend('topright', names(label), text.col=as.character(label), bty='n')

#corrplot.mixed(M, lower="color", upper="number",
#         #order="original", tl.col=tlcol, tl.cex=0.3, tl.srt=45,
#         order="hclust", tl.col=tlcol, tl.cex=0.3, tl.srt=45,
#         addrect=4, col=colorRampPalette(col)(200))
#
#legend('topleft', names(label), text.col=as.character(label), bty='n')

# Combine with a significant test. Provide the matrix of correlation but also a matrix of p-val for all correlations computed
# sig.level is the top value a p-val is considered significant
corrplot(M, p.mat= sig[[1]], sig.level=sig.lev,
         order="hclust", method="color", type="lower",
         tl.col=tlcol, tl.cex=0.5,
         insig="pch", pch.cex=0.7,
         addrect=4, col=colorRampPalette(col)(200))
legend('topright', c(names(label), paste0('Confidence level: ', conf.lev, ", p-val < ", sig.lev)),
       text.col=c(as.character(label), 'black'), bty='n')

#corrplot(M, method="color", type="lower",
#         order="FPC", tl.col="black", tl.cex=0.5, 
#         col=colorRampPalette(col)(200))
#legend("topright", 'Ordered based on First Principal Component', col="black", bty='n')
#

dev.off()# }}}


