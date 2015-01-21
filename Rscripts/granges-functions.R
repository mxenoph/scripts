source("~/source/Rscripts/functions.R")

# chr should be either a file name (e.g. "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10_noHead.genome") or a dataframe
add_seqlengths_to_gr <- function(gr, chr){# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    # if not data.frame it should be a file pointing 
    if(class(chr) != "data.frame"){
        if( !file.exists(chr) ) { 
            stop(chr, " doesn't exist.")
        } else {
            chr <- read.table(chr, sep="\t", row.names=1)
            seqlengths(gr) <- chr[names(seqlengths(gr)),]
        }
    } else {
        name <- which(apply(chr, 2, function(y) any(grepl("chr\\d+", y, perl=T))))
        # if chr ar not a column
        if (length(name) == 0) {
            if( !any(grepl("chr\\d+", rownames(chr), perl=T)) ) stop("Dataframe does not contain chromosome names.")
            seqlengths(gr) <- chr[names(seqlengths(gr)),]
        } else {
            rownames(chr) <- chr[[name]]
            # removing the column with chr names as added the rownames in the above line
            chr[[name]] <- NULL
            seqlengths(gr) <- chr[names(seqlengths(gr)),]
        }
    }

    return(gr)
}# }}}

change_seqnames<-function(gr){# {{{
    # List of packages to load
    x <- c('GenomicRanges', 'gtools')
    lapply(x, suppressMessages(library), character.only=T)

    chr_to_num <- function(gr){# {{{
      level <- seqlevels(gr)[grepl("^chr[0-9XYM]", seqlevels(gr))]
      chr <- mixedsort(level)
      new <- mixedsort(gsub('chr', '', level))
      names(new) <- paste("chr", new, sep = "")
      new["chrM"] <- "MT"

      seqlevels(gr, force=TRUE) <- chr
      gr <- keepSeqlevels(gr, chr)
      gr <- renameSeqlevels(gr, new)
      return(gr)
    }# }}}

    num_to_chr <- function(gr){# {{{
      level <- seqlevels(gr)[grepl("^[0-9XY]|MT", seqlevels(gr))]
      level[match("MT", level)] <- "M"
      chr <- mixedsort(level)
      new <- mixedsort(gsub("^", 'chr', level))
      names(new) <- mixedsort(level)

      seqlevels(gr, force=TRUE) <- chr
      gr <- keepSeqlevels(gr, chr)
      gr <- renameSeqlevels(gr, new)
      return(gr)
    }# }}}

    if( any(grepl('chr', seqlevels(gr)) )){
      return(chr_to_num(gr))
    }
    else{
      return(num_to_chr(gr))
    }

}# }}}

macs_to_granges <- function(peaks){# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    peaks <- read.table(peaks, header=T, sep="\t")
    gr <- with(peaks, GRanges(
                 seqnames= chr,
                 range= IRanges(start=start, end=end),
                 strand= "*",
                 wd= length,
                 score= X.10.log10.pvalue.,
                 FE= fold_enrichment,
                 summit= start + summit,
                 tags= tags))
                 
    if('FDR...' %in% colnames(peaks)) values(gr)[['fdr']] <- with(peaks, FDR...) # fdr column is not available if the experiment had no control
    return(gr)
}
# }}}

macs2_to_granges <- function(peaks){# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    peaks <- read.table(peaks, header=T, sep="\t")
    gr <- with(peaks, GRanges(
                 seqnames= chr,
                 range= IRanges(start=start, end=end),
                 strand= "*",
                 wd= length,
                 tags= pileup,
                 score= X.log10.pvalue.,
                 FE= fold_enrichment))
                 
    if('X.log10.qvalue.' %in% colnames(peaks)) values(gr)[['fdr']] <- with(peaks, 'X.log10.qvalue.') # fdr column is not available if the experiment had no control
    return(gr)
}
# }}}

name_gr <- function(gr) paste(seqnames(gr), start(gr), end(gr), sep=':')
get_score <- function(gr, value) values(gr)[[value]]

intersect_multi <- function(grl){# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    if(length(grl) < 2) return(grl[[1]])
    gr <- grl[[1]]
    for (j in seq(2, length(grl), 1)){
        gr <- intersect(gr, grl[[j]])
    }
    names(gr) <- name_gr(gr)
    return(gr)
}# }}}

get_hits_scores <- function(query, subj, value) {# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    # returns a vector of length(common) containing index for hits in x
    ov <- findOverlaps(query, subj, select='first')
    names(ov) <- names(query)
    # subset based on those that have a hit
    ov[!is.na(ov)] <- values(subj[as.numeric(ov[!is.na(ov)])])[[value]]
    return(ov)
}# }}}


