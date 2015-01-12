modules::import("annotation-functions", attach = T)
get_annotation('mm10')

source("~/source/Rscripts/functions.R")
source("~/source/Rscripts/granges-functions.R")

x <- c('rtracklayer', 'GenomicRanges', 'parallel')
lapply(x, suppressMessages(library), character.only=T)

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

setwd("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm9/macs")
peaks <- lapply(list.files(pattern=".*[^e]_peaks.xls"), macs2GRanges)
names(peaks) <- list.files(pattern=".*[^e]_peaks.xls")
# ranges of width 2Kb with center being the TSS
tss_window <- import("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/Mus_musculus.GRCm38.70_conv.tss_nocontig.gtf")
tss <- GRanges(seqnames(tss_window),
               IRanges(start = start(tss_window) + (width(tss)/2) -1, 
                       width = 1),
               strand = strand(tss_window))

chd4_wt <- GRanges(seqnames(peaks[[1]]),
                   IRanges(start=values(peaks[[1]])[['summit']], end=values(peaks[[1]])[['summit']]),
                   strand="*", 
                   fe = values(peaks[[1]])[['FE']])
ov <- findOverlaps(tss_window, chd4_wt)

positions <- mclapply(1:length(ov), function(x, query, subj){# {{{
                      pair <- ov[x]
                      p <- get_position(query[queryHits(pair)], subj[subjectHits(pair)])
                      d <- data.frame(summit = p, fe = values(subj[subjectHits(pair)])[['fe']])
                      return(d)
               }, tss, chd4_wt)
positions <- do.call(rbind, positions)# }}}


pdf('test_bins.pdf')

ggplot(data=signal$sums, aes(x=range, y=value)) + geom_point() + scale_x_continuous(breaks=seq(-1000, 1000, by=200))
dev.off()


positions <- mclapply(1:100, function(x, query, subj, centers){ # {{{
                      ov <- findOverlaps(query[x], subj)

                     if ( length(queryHits(ov)) != 0 ){
                         # If more than 1 hits in subj then I get a GRanges object
                         # thus the function returns a vector of positions -- works fine
                         p <- get_position(centers[x], subj[subjectHits(ov)])
                         p <- p + 1000
#                         names(p) <- values(query[x])[['gene_id']]
#                         print(paste0('Position ', p, "  x=", x))

#                         sequence <- matrix(0, ncol=2000)
#                         sequence[1, p] <- 1
#                         rownames(sequence) <- values(query[x])[['gene_id']]
                     } else {
#                         sequence <- matrix(0, ncol=2000)
#                         rownames(sequence) <- values(query[x])[['gene_id']]
                         p <- 0
                     }
                     #return(sequence)
                     return(p)
}, tss_window, chd4_wt, tss)

#positions <- as.data.frame(do.call(rbind, positions))
#colnames(positions) <- 1:2000
# }}}

summarize_signal <- function(positions){# {{{
    x <- c('reshape')
    lapply(x, suppressMessages(library), character.only=T)

    sums <- colSums(positions)
    sums <- melt(sums)
    sums$range <- as.numeric(rownames(sums)) - 1000

    return(list(sums = sums))
}
# }}}
signal <- summarize_signal(positions)
