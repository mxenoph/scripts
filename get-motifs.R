#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")

parser <-  ArgumentParser(description="Perform motif analysis")
parser$add_argument('-x', '--peaks', metavar= "file", required='True', type= "character", help= "MACS generated excel file for called peaks")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
parser$add_argument('-d', '--dist', type= "integer", default= 50,  help= "Distance around the summit to look for motifs. Default: 50bp")
parser$add_argument('-n', '--npeaks', type= "integer", default= 300,  help= "Number of peaks to look in for motifs. Default: 300")
parser$add_argument('-m', '--mode', type= "character", default= 'top',  help= "Look at the top/random/bottom peaks. Default: top")

args <- parser$parse_args()

parent_path <- file.path(args$out, 'meme', basename(file_path_sans_ext(args$peaks)))
output_path <- file.path(parent_path, paste0(args$mode, args$npeaks))
plot_path <- file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

# function in annotation-functions.R
annotationInfo(tolower(args$assembly))

# define variables
memeDatabasePath <- '/nfs/research2/bertone/user/mxenoph/common/meme/motif_databases'

# }}}

# Import libraries & Turn off warning messages for loading the packages-globally# {{{
options(warn=-1)
suppressMessages(library(BSgenome))
# Load BSgenome based on assembly
suppressMessages(library(paste0("BSgenome.Mmusculus.UCSC.", args$assembly), character.only=TRUE))
source("~/local/granges.functions.R")
#Turn warnings back on
options(warn=0)
# }}}

# Functions # {{{

# Functions by konrad Rudolph# {{{
# Used like that in Konrad's parser of pspm
map <- base::Map

# This uses R's peculiarities in argument matching explained here:
# <http://stat.ethz.ch/R-manual/R-devel/doc/manual/R-lang.html#Argument-matching>
# `.expr` starts with a dot to allow `expr` being used in the actual
# expression.
let <- function (.expr, ...)
        eval(substitute(.expr), list2env(list(...), parent = parent.frame()))

#' Create a closure over a given environment for the specified formals and body.
closure <- function (formals, body, env)
        eval(call('function', as.pairlist(formals), body), env)
#' Create a list of empty symbols, with names set
symlist <- function (names)
        setNames(Map(function (p) quote(expr = ), names), names)

#' A shortcut to create a function
#'
#' @note Using \code{.(args = body)} is analogous to using
#' \code{function (args) body} with one exception: \code{.} arguments do not
#' support defaults.
#' Since its purpose is mainly for lambdas in higher-order list functions, this
#' functionality is not needed.
. <- function (...) {
    args <- match.call(expand.dots = FALSE)$...
    last <- length(args)
    params <- symlist(c(args[-last], names(args)[[last]]))
    if (length(args) > 1 && length(params) != length(args))
        stop('Must be of the form `fun(a, b = expr)`')
    for (arg in args[-last])
        if (! is.name(arg))
            stop('Invalid argument specifier: ', arg)

    closure(params, args[[last]], parent.frame())
}# }}}

# Read in peaks and get summits# {{{
getLocus <- function(file){
    peaks <- macs2GRanges(file)
    peaks <- peaks[order(-values(peaks)$score)]
    peaks <- add.seqlengths(peaks, chr_size)

    #why just looking at + strand?
    summits <- GRanges(seqnames= seqnames(peaks),
                       range= IRanges(start= values(peaks)$summit - args$dist,
                                      end= values(peaks)$summit + args$dist),
                       strand= "+",
                       seqlengths= seqlengths(peaks))

    values(summits) <- values(peaks)
    #trimming ranges that exceed the chromosome length
    summits <- trim(summits)

    #Because of trimming I might end up with ranges of 0 length so get rid of those
    summits <- summits[width(summits) != 0]

    summits <- summits[order(-values(summits)$score)]
    names(summits) <- paste(seqnames(summits), start(summits), end(summits), sep=":")

    return(summits)
} # }}}

# Get sequences for peaks # {{{
seq_motifs <- function(summits, n, output_path){
    seqs <- getSeq(Mmusculus, summits[1:n], as.character=FALSE)
    # make a dataframe holding information on FE, FDR, sequence
    names(seqs) <- paste(names(summits)[1:n],
                            as.vector(values(summits)[['FE']][1:n]),
                            as.vector(values(summits)[['fdr']][1:n]), sep='|')

    fasta <- file.path(output_path, paste0('top', n, '.fa'))
    writeXStringSet(seqs, file= fasta, "fasta", append=FALSE)
}# }}}

# Parse motifs and E-values from MEME output. From Konrad Rudolph  # {{{
parsePspm <- function (results) {
    # from Konrad what the hell is the output?
    # Parse output raw text file and retrieve PSSMs for all motifs.
    # Note: we could also parse the XML file using XPath but the XML output file
    # is quite frankly not very nice. Parsing the raw text is easier. Fail.
    extractHeader <- function (name, header)
        let(match = regexpr(sprintf('(?<=%s= )\\S+', name), header, perl = TRUE),
            as.numeric(regmatches(header, match)))

    file <- file.path(results, 'meme.txt')
    lines <- readLines(file)
    start <- grep('^\tMotif \\d+ position-specific probability matrix', lines) + 2
    length <- extractHeader('w', lines[start])
    nsites <- extractHeader('nsites', lines[start])
    evalue <- extractHeader('E', lines[start])
    matrix <- map(.(start, end = read.table(text = paste(lines[start : end],
                                                         collapse = '\n'),
                                            col.names = c('A', 'C', 'G', 'T'))),
                  start + 1, start + length)

    map(.(matrix, nsites, evalue =
          let(str = list(matrix = matrix, nsites = nsites, e = evalue),
              structure(str, class = 'pspm'))),
        matrix, nsites, evalue)
}# }}}

# Adapted from Konrad# {{{
# No need to filter the sequence
memeFileContent <- function (pspm, fasta) {
    # Format documented at
    # <http://meme.nbcr.net/meme/doc/meme-format.html#min_format>
    version <- 'MEME version 4.9\n'

    inputFasta <- readLines(fasta)

    # Remove Fasta headers and merge lines
    inputFasta <- paste(inputFasta[! grepl('^>', inputFasta)], collapse = '')
    # Normalise casing and remove Ns
    inputFasta <- gsub('N', '', toupper(inputFasta))
    freqs <- table(strsplit(inputFasta, ''))

    frequencies <- sprintf('Background letter frequencies:\n%s\n',
    paste(mapply(paste, names(freqs), freqs), collapse = ' '))

    motif <- capture.output(write.table(pspm$matrix, col.names = FALSE,
    row.names = FALSE))

    c(version,
      # Alphabet omitted (can be inferred)
      # Strand omitted (can be inferred)
      frequencies,
      pspmHeader(pspm),
      motif)
}# }}}

consensusSequence <- function (pspm){# {{{
    with(pspm, paste(colnames(matrix)[
                     vapply(as.data.frame(t(matrix)), which.max,
                            numeric(1))], collapse = ''))
}# }}}

pspmHeader <- function (pspm){# {{{
    paste(sprintf('MOTIF %s', consensusSequence(pspm)),
          sprintf('letter-probability matrix: alenght= %s w= %s nsites= %s E= %e',
                  ncol(pspm$matrix), nrow(pspm$matrix), pspm$nsites, pspm$e),
          sep = '\n\n')
}# }}}

# Run MEME# {{{
meme <- function(output) {
    memeBin <- 'meme'
    # nmotifs = max number of motifs to find
    # minsites= min number of sites for each motif
    # minw = min width of motif, maxw = max n of motifs
    # revcomp =  allow sites on + or - strands
    # maxsize = minimum dataset size in characters
    # -evt = stop if motif E-value greater than <evt>
    # possibly don't limit to 3 motifs but filter by E-value
    fasta <- list.files(output, pattern = ".fa", full.name = TRUE)

    system(sprintf('%s %s -dna -oc %s -maxsize %s -mod zoops -nmotifs 3 -evt 0.05 -minw 6 -maxw 35 -revcomp',
                   memeBin,
                   fasta,
                   file.path(output),
                   round(as.numeric(system(paste0("wc -c ", fasta, " | awk -F' ' '{print $1}'"), intern=TRUE)) -2)
                   ))


    parsePspm(file.path(output))
} # }}}

# Run tomtom# {{{
tomtom <- function (pspm, outputPath, databases) {
    motif <- consensusSequence(pspm)
    inputFile <- file.path(outputPath, paste0(motif, '.meme'))
    fasta <- list.files(pattern = "*.fa", output_path, full.names = TRUE)
    writeLines(memeFileContent(pspm, fasta), inputFile)

    if (missing(databases))
        databases <- file.path(memeDatabasePath,
                               c('JASPAR_CORE_2014_vertebrates.meme',
                                 'uniprobe_mouse.meme'))
    
    system(sprintf('%s -no-ssc -oc %s -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10 %s %s',
                   'tomtom',
                   file.path(outputPath, motif),
                   inputFile,
                   databases))

    extractMotifName <- function (header){
        match <- regexpr(sprintf('\\t(\\S+)\\t'), header, perl = TRUE)
        id <- gsub('\t', '', regmatches(header, match))
        name <- system(sprintf('grep %s %s | cut -d \' \' -f 3 | tr -d \'\n\'',
                               id,
                               file.path(memeDatabasePath, c('JASPAR_CORE_2014_vertebrates.meme',
                                                             'uniprobe_mouse.meme'))
                               ), intern = TRUE)
        return(name)
    }
   lines <- readLines(file.path(outputPath, motif, 'tomtom.txt'), n= -1L)
   lines <- lines[2:length(lines)]
   factors <- as.vector(unlist(map(extractMotifName, lines)))
   combined <- as.vector(unlist(map(function(x, y)
                                    paste(x, y, sep = '\t'),
                                    lines, factors)))

   write.table(file = file.path(outputPath, motif, 'tomtom.txt'), combined,
               col.names = paste(c('QueryID', 'TargetID', 'Offset', 'pvalue',
                                   'Evalue', 'qvalue', 'Overlap', 'QueryConsensus',
                                   'TargetConsensus', 'Orientation', 'Factor'),
                                 collapse = '\t'),
               row.names = FALSE,
               quote = FALSE)

   motifId <- read.table(pipe(sprintf('cut -f 2 %s',
                                     file.path(outputPath, motif, 'tomtom.txt'))),
                         header = TRUE)
   return(list(motif= motif, motif_id= motifId, motif_name = factors))
}# }}}

# Get sequences for peaks # {{{
seq_general <- function(regions){
    seqs <- getSeq(Mmusculus, regions, as.character=TRUE)
    # make a dataframe holding information on FE, FDR, sequence
    names(seqs) <- names(regions)

#    fasta <- file.path(output_path, paste0('top', n, '.fa'))
#    writeXStringSet(seqs, file= fasta, "fasta", append=FALSE)
    return(seqs)
}# }}}

# Parse motifs and E-values from MEME formatted databases. Adapted from Konrad Rudolph  # {{{
parsePspmDb <- function (databasePath) {
    # from Konrad what the hell is the output?
    # Parse output raw text file and retrieve PSSMs for all motifs.
    # Note: we could also parse the XML file using XPath but the XML output file
    # is quite frankly not very nice. Parsing the raw text is easier. Fail.
    extractHeader <- function (name, header)
        let(match = regexpr(sprintf('(?<=%s= )\\S+', name), header, perl = TRUE),
            as.numeric(regmatches(header, match)))
    extractName <- function (header)
        gsub('MOTIF ', '', let(match = regexpr('^MOTIF \\S+', header, perl = TRUE),
                              as.character(regmatches(header, match))))

    lines <- readLines(databasePath)
    start <- grep('^MOTIF \\S+', lines) + 2
    name <- extractName(lines[start-2])
    length <- extractHeader('w', lines[start])
    evalue <- extractHeader('E', lines[start])
    matrix <- map(.(start, end = read.table(text = paste(lines[start : end],
                                                         collapse = '\n'),
                                            col.names = c('A', 'C', 'G', 'T'))),
                  start + 1, start + length)

    pspm <- map(.(matrix, evalue =
                  let(str = list(name = name, matrix = matrix, e = evalue),
                      structure(str, class = 'pspm'))),
                matrix, evalue)
    names(pspm) <- name
    return(pspm)
}# }}}

# Scan a string/sequence for PWM matches# {{{
extractPspmMatches <- function(pspm, sequence, threshold, summit){
    # Get all matches for pspm
    # needs a numeric matrix with rows A,C,T,G, the minimum score for counting a match
    # that could be a percentage of the highest possible score or a number
   match <- matchPWM(pspm, sequence, min.score = threshold)
#   rmatch <- matchPWM(reverseComplement(pspm), sequence, min.score = threshold)

    if (length(match) == 0)
        found <- cbind(NA, NA, 0)
    else
        found <- cbind(start(match),
                       as.character(match),
                       length(start(match)))
    colnames(found) <- c('start', 'seq', 'hits')

    # Keep only the match that is closest to the summit, that's why we use abs()
    found <- found[order(abs(as.integer(found[,'start']) - summit), decreasing = FALSE)[1], ]
    return(found)
}# }}}

# Scan peaks for a motif # {{{
scanPeaks <- function(peaks, pspm, threshold){
    # Get sequences for all peaks
    sequences <- elementMetadata(peaks)$seqs
    motifs <- sapply(1:length(sequences),
                     function(indx){
                         extractPspmMatches(pspm, sequences[indx], threshold, elementMetadata(peaks)$summit[indx])})
    motifs <- t(motifs)
    rownames(motifs) <- names(peaks)
    return(motifs)
}# }}}

# Get the motif profile around the peak summits# {{{
getPeakProfile <- function(peaks, pspm, window){
    summits <- GRanges(seqnames = seqnames(peaks),
                       range = IRanges(start = elementMetadata(peaks)[['summit']] - window,
                                       end = elementMetadata(peaks)[['summit']]) + window,
                       strand = '+',
                       summit = elementMetadata(peaks)[['summit']])
    names(summits) <- paste(seqnames(summits), start(summits), end(summits), sep=':')

    elementMetadata(summits)$seqs <- seq_general(summits)
    scanned <- scanPeaks(summits, pspm, "80%")
    # Keep only those ranges for which a motif was matched
    scanned <- scanned[!is.na(scanned[,'start']), ]
    # Get position covered by the motif
    positions <- sapply(as.integer(scanned[,'start']), function(x) x + ncol(pspm)-1)
    positions <- data.frame(table(unlist(positions)))
    names(positions) <- c('position','frequency')
    # shift position relative to summit
    positions$position <- as.integer(as.character(positions$position)) - window
    # ensure all positions are present in the output
    profile <- data.frame('position' = c(-window:window),
                          'frequency' = positions$frequency[match(-window:window, positions$position)])
    profile[is.na(profile)] <- 0
    return(profile)
}# }}}

runPeakProfileOnAll <- function(peaks, motifs){# {{{
    pspmDb <- parsePspmDb(file.path(memeDatabasePath, "JASPAR_CORE_2014_vertebrates.meme"))
    pspmDb <- c(pspmDb, parsePspmDb(file.path(memeDatabasePath, "uniprobe_mouse.meme")))

    library(gridExtra)
    pspms <- sapply(1:length(motifs), function(indx, db){
                    ids <- levels(unlist(motifs[[indx]]$motif_id))
                    mnames <- levels(unlist(motifs[[indx]]$motif_name))
                    results <- lapply(ids, function(x){
                                      pspm <- t(db[[x]]$matrix)
                                      scanned <- scanPeaks(peaks, pspm, "80")
                                      profile <- getPeakProfile(peaks, pspm, 200)

                                      result <- list(hits = scanned, profile = profile)
                                      return(result)
                                       })
                    names(results) <- paste(ids, mnames, sep=':')

                    p <- vector("list", length(results)*2)
                    for (i in seq(1, length(results)*2, 2)) {
                        x <- ceiling(i/2)
                        tmp <- as.data.frame(results[[x]][['hits']])
                        # Convert factor to numeric values so that x-axis is sorted
                        tmp$hits <- as.numeric(as.character(tmp$hits))
                        a <- ggplot(tmp, aes(hits))
                        a <- a + geom_bar(stat = "bin", binwidth = 1) + geom_hline(yintercept=0, colour="white", size=0.5)
                        a <- a + labs(title = paste0('Motif', indx, ':', names(results)[x]))

                        b <- ggplot(as.data.frame(results[[x]][['profile']]), aes(x=position, y=frequency))
                        b <- b + geom_point(aes(colour = cut(position, breaks=c(-Inf,-50,50, Inf)),
                                                alpha = frequency))
                        b <- b + guides(colour = FALSE) + theme(legend.position='none')
                        p[[i]] <- a
                        p[[i+1]] <- b
                    }
                    # call marrangeGrob instead of grid.arrange so that the plots span
                    # multiple pdf pages as in answer in 
                    # <http://stackoverflow.com/questions/19059826/multiple-graphs-over-multiple-pages-using-ggplot>
                    # top=NULL removes page numbers
                    multipage <- do.call(marrangeGrob, c(p, nrow = 3, ncol = 2, top = NULL))
                    ggsave(file.path(plot_path, paste0('summary-motif', indx, '.pdf')), multipage, width = 8.3, height = 11.7)
                          }, pspmDb)
}# }}}

# end of fold sections # }}}

peaks <- macs2GRanges(args$peaks)
peaks <- add.seqlengths(peaks, chr_size)
peaks <- peaks[order(-values(peaks)$score)]
names(peaks) <- paste(seqnames(peaks), start(peaks), end(peaks), sep=":")
elementMetadata(peaks)$seqs <- seq_general(peaks)

loci <- getLocus(args$peaks)
seq_motifs(loci, args$npeaks, output_path)
pspm <- meme(output_path)
motifs <- map(tomtom, pspm, output_path)

runPeakProfileOnAll(peaks, motifs)

# TODO: seq_motifs on random and bottom peaks {{{
# should be a new separate function
# print(paste0("Getting sequences for bottom ", n, " peaks"))
# writeXStringSet(seqs[names(seqs) %in% names(summits)[length(summits)-(n-1):length(summits)]], file= paste(out, prefix, ".bottom", n, ".fa", sep=""), "fasta", append=FALSE)
#
# print(paste0("Getting sequences for random ", n, " peaks"))
# index <- sample(1:length(summits), size=n, replace=FALSE)
# writeXStringSet(seqs[names(seqs) %in% names(summits)[index]], file= paste(out, prefix, ".random", n, ".fa", sep=""), "fasta", append=FALSE)
## }}}

# ggplot specific# {{{
require(grid)
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
    dots <- list(...)
    n <- length(dots)
    if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
    if(is.null(nrow)) { nrow = ceiling(n/ncol)}
    if(is.null(ncol)) { ncol = ceiling(n/nrow)}
    ## NOTE see n2mfrow in grDevices for possible alternative
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
    ii.p <- 1
    for(ii.row in seq(1, nrow)){
        ii.table.row <- ii.row  
        if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
        for(ii.col in seq(1, ncol)){
            ii.table <- ii.p
            if(ii.p > n) break
            print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
            ii.p <- ii.p + 1
        }
    }
}# }}}
