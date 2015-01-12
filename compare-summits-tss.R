#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
x <- c('argparse', 'tools')
lapply(x, suppressMessages(library), character.only=T)
parser <-  ArgumentParser(description="Compare summits to TSS")
parser$add_argument('-p', '--path', metavar= "file", required='True', type= "character", help= "Path to look for MACS generated excel file for called peaks")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
parser$add_argument('-r', '--replotting', type= "character",  help= "To make the plots from an Rdata object")

args <- parser$parse_args()

output_path <- file.path(args$out, 'processed')
plot_path <- file.path(args$out, 'plots')
dir.create(output_path)
dir.create(plot_path)

# Load this before the other packages so that the import function from rtracklayer 
# is not masked
modules::import("annotation-functions", attach = T)
get_annotation(args$assembly)

source("~/source/Rscripts/granges-functions.R")
x <- c('rtracklayer', 'GenomicRanges', 'parallel')
lapply(x, suppressMessages(library), character.only=T)# }}}

# Functions
# x is the feature we want stuff centered to
get_position <- function(x, y){ # {{{
    # e.g. with x=tss=10, "+" directionality and summits 5 and 16 the positions will be -5 and 6 with
    # negative sign indicating upstream
    if(as.character(strand(x))=='+'){
        position <- start(y) - start(x)
    } else {
        # e.g. with x=tss=10, "-" directionality and summits 4 and 17
        # the positions will be 6 and -7
        position <- start(x) - start(y)
    }
    return(position)
}
# }}}

plot_scatter <- function(positions) {
    for(x in 1:length(positions)) {
        df <- positions[[x]]
        p <- ggplot(df, aes(x= summit, y= fe)) + geom_point(col="grey50")
        p <- p + scale_x_continuous(breaks=seq(-1000, 1000, by=200))
        p <- p + geom_rug()
        p <- p + ggtitle(names(positions)[x])
        p <- p + annotate("text", y= max(df$fe), x =1000,
                       # paste0 will not work cause it's not a plotmath function
                       label= paste("n[peaks_in_window]", "==", nrow(df)),
                       parse = TRUE,
                       col = "red",
                       hjust=1)
        print(p)
    }
}

#setwd("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/macs")

main <- function(){ #{{{
    # Preparing the data
    # Select the _peaks.xls files in the path and not the negative_peaks.xls
    peaks <- lapply(list.files(pattern=".*[^e]_peaks.xls"), macs2GRanges)
    names(peaks) <- list.files(pattern=".*[^e]_peaks.xls")

    summits <- mclapply(peaks, function(gr){# {{{
             gr <- GRanges(seqnames(gr),
                       IRanges(start = values(gr)[['summit']], end = values(gr)[['summit']]),
                       strand = "*",
                       fe = values(gr)[['FE']])
    })# }}}

    # tss and tss_window objects loaded from get_annotation(args$assembly) call

    # ranges of width 2Kb with center being the TSS
    #tss_window <- import("/nfs/research3/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/Mus_musculus.GRCm38.70_conv.tss_nocontig.gtf")
    #tss <- GRanges(seqnames(tss_window),
    #               IRanges(start = start(tss_window) + (width(tss_window)/2) -1,
    #                       width = 1),
    #               strand = strand(tss_window))# }}}

    # Calculate position of summit relative to tss and keep fe for peak# {{{
    positions <- mclapply(summits, function(summit) {
             ov <- findOverlaps(tss_window, summit)

             positions <- mclapply(1:length(ov), function(x, query, subj){
                                   pair <- ov[x]
                                   p <- get_position(query[queryHits(pair)], subj[subjectHits(pair)])
                                   d <- data.frame(summit = p, fe = values(subj[subjectHits(pair)])[['fe']])
                                   return(d)
                           }, tss, summit)

             positions <- do.call(rbind, positions)
             

             return(positions)
    })# }}}

}
# }}}


# Printing
if( !file.exists(file.path(output_path, 'peak-summitsVsTSS.Rda')) ) {
    main()
} else {
    if ()
}

pdf('peak-summitsVsTSS.Rda')
dev.off()
