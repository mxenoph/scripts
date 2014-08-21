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
dir.create(output_path, recursive= TRUE)

# function in annotation-functions.R
annotationInfo(tolower(args$assembly))

# define variables
memeDatabasePath <- '/nfs/research2/bertone/user/mxenoph/common/meme/motif_databases'

# }}}

# Import libraries & Turn off warning messages for loading the packages-globally# {{{
options(warn=-1)
suppressMessages(library(BSgenome))
suppressMessages(library("BSgenome.Mmusculus.UCSC.mm9"))
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

    #system(sprintf('%s %s -dna -oc %s -maxsize %s -mod zoops -nmotifs 3 -evt 0.05 -minw 6 -maxw 35 -revcomp',
    #               memeBin,
    #               paste0(output, '.fa'),
    #               file.path(output),
    #               round(as.numeric(system(paste0("wc -c ", paste0(output, '.fa'), " | awk -F' ' '{print $1}'"), intern=TRUE)) -2)
    #               ))
    #

    parsePspm(file.path(output))
} # }}}

# Run tomtom# {{{
tomtom <- function (pspm, outputPath, databases) {
    motifName <- consensusSequence(pspm)
    inputFile <- file.path(outputPath, paste0(motifName, '.meme'))
    fasta <- list.files(pattern = "*.fa", output_path, full.names = TRUE)
    writeLines(memeFileContent(pspm, fasta), inputFile)

    if (missing(databases))
        databases <- file.path(memeDatabasePath,
                               c('JASPAR_CORE_2014_vertebrates.meme',
                                 'uniprobe_mouse.meme'))
    
    system(sprintf('%s -no-ssc -oc %s -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10 %s %s',
                   'tomtom',
                   file.path(outputPath, motifName),
                   inputFile,
                   databases))

    extractMotifId <- function (header){
        match <- regexpr(sprintf('\\t(\\S+)\\t'), header, perl = TRUE)
        id <- gsub('\t', '', regmatches(header, match))
        name <- system(sprintf('grep %s %s | cut -d \' \' -f 3 | tr -d \'\n\'',
                               id,
                               file.path(memeDatabasePath, c('JASPAR_CORE_2014_vertebrates.meme',
                                                             'uniprobe_mouse.meme'))
                               ), intern = TRUE)
        return(name)
    }
   lines <- readLines(file.path(outputPath, motifName, 'tomtom.txt'), n= -1L)
   lines <- lines[2:length(lines)]
   factors <- as.vector(unlist(map(extractMotifId, lines)))
   combined <- as.vector(unlist(map(function(x, y)
                                    paste(x, y, sep = '\t'),
                                    lines, factors)))

   write.table(file = file.path(outputPath, motifName, 'tomtom.txt'), combined,
               col.names = paste(c('QueryID', 'TargetID', 'Offset', 'pvalue',
                                   'Evalue', 'qvalue', 'Overlap', 'QueryConsensus',
                                   'TargetConsensus', 'Orientation', 'Factor'),
                                 collapse = '\t'),
               row.names = FALSE,
               quote = FALSE)
   print('tomtom function')
}# }}}


# end of fold sections # }}}

#loci <- getLocus(args$peaks)
#seq_motifs(loci, args$npeaks, output_path)
pspm <- meme(output_path)
map(tomtom, pspm, output_path)


# Get sequences for peaks # {{{
seq_general <- function(region){
    seqs <- getSeq(Mmusculus, region, as.character=FALSE)
    # make a dataframe holding information on FE, FDR, sequence
    names(seqs) <- names(regions)

#    fasta <- file.path(output_path, paste0('top', n, '.fa'))
#    writeXStringSet(seqs, file= fasta, "fasta", append=FALSE)
    return(seqs)
}# }}}

peaks <- macs2GRanges(args$peaks)
peaks <- add.seqlengths(peaks, chr_size)
peaks <- peaks[order(-values(peaks)$score)]
names(peaks) <- paste(seqnames(peaks), start(peaks), end(peaks), sep=":")


# TODO: seq_motifs on random and bottom peaks {{{
# should be a new separate function
# print(paste0("Getting sequences for bottom ", n, " peaks"))
# writeXStringSet(seqs[names(seqs) %in% names(summits)[length(summits)-(n-1):length(summits)]], file= paste(out, prefix, ".bottom", n, ".fa", sep=""), "fasta", append=FALSE)
#
# print(paste0("Getting sequences for random ", n, " peaks"))
# index <- sample(1:length(summits), size=n, replace=FALSE)
# writeXStringSet(seqs[names(seqs) %in% names(summits)[index]], file= paste(out, prefix, ".random", n, ".fa", sep=""), "fasta", append=FALSE)
## }}}

#editfunc <- function(a){# {{{
## Motif scan
## read the NF-$\kappa$B PWM from a text file
#NFKB.pwm <- read.table("/nfs/training/PeakCalling/NFKB_pwm_meme.txt",
#                                              sep="\t", header=FALSE, row.names=1)
#NFKB.pwm <- t(NFKB.pwm)
#rownames(NFKB.pwm) <- c("A","C","G","T")
## view the motif matrix
#NFKB.pwm
## plot the motif logo
#library(seqLogo)
#seqLogo(NFKB.pwm, ic.scale=TRUE)
#
## obtain the sequences under the peaks
#elementMetadata(macspeaks)$seqs <- getSeq(Hsapiens, macspeaks, as.character=TRUE)
#
#
## function to scan a string for PWM matches at the specified threshold by calling matchPWM()
## keeps only the closest of multiple hits
#mymatchPWM <- function (pwm, myseq, threshold, summit) {
#        # get all matches of PWM
#        mymatch <- matchPWM(pwm, myseq, min.score=threshold)
#    # collect starts/seqs into matrix (if any)
#if (length(mymatch)==0) {
#        found <- cbind(NA,NA,0)
#} else {
#        found <- cbind(start(mymatch), as.character(mymatch), length(start(mymatch)))
#}
#colnames(found) <- c("start","seq","nr")
## keep only the match that is closest to the summit
#found <- found[order(abs(as.integer(found[,1])-summit),decreasing=FALSE)[1],]
## return matrix
#return(found)
## function to call mymatchPWM() on each peak in a set
## returns a matrix with motif information per peak
#ScanPeaks <- function(peak.GR, pwm, threshold) {
#        # get all peak sequences
#        myseqs <- elementMetadata(peak.GR)$seqs
#    # get summit positions (relative to peak coordinates)
#    summits <- elementMetadata(peak.GR)$maxpos - start(peak.GR)
#        # apply mymatchPWM() to all peaks in the set
#        motifmatrix <- sapply(1:length(myseqs), function(x) mymatchPWM(pwm, myseqs[x], threshold, summits[x]))
#        motifmatrix <- t(motifmatrix)
#            # set peak IDs as rownames
#            rownames(motifmatrix) <- names(peak.GR)
#            # return motif matrix
#            return(motifmatrix)
#}
#
## scan all peak sequences for the motif (using a cutoff of 80%)
#### This step is likely to take ~ 10 min, so at this point it's convenient to have a coffee :P ###
#macs.motifs <- ScanPeaks(macspeaks, as.matrix(NFKB.pwm), "80%")
#useq.motifs <- ScanPeaks(useqpeaks, as.matrix(NFKB.pwm), "80%")
#chipseq.motifs <- ScanPeaks(chippeaks, as.matrix(NFKB.pwm), "80%")
## add motif information (number of motifs per peaks) to the peak GRanges object
#elementMetadata(macspeaks)$motif.no <- as.integer(macs.motifs[,3])
#elementMetadata(useqpeaks)$motif.no <- as.integer(useq.motifs[,3])
#elementMetadata(chippeaks)$motif.no <- as.integer(chipseq.motifs[,3])
## plot the proportion of peaks with 0, 1 or more motifs in the three peak sets
#par(mfrow=c(1,3))
#pie( c(sum(as.numeric(elementMetadata(macspeaks)$motif.no==0)),
#              sum(as.numeric(elementMetadata(macspeaks)$motif.no==1)),
#                     sum(as.numeric(elementMetadata(macspeaks)$motif.no>1))),
#         labels=c("0","1","2+"), main=list("MACS v2",cex=1.5),
#              col=c("white","lightgrey","darkgrey")
#         )
#pie( c(sum(as.numeric(elementMetadata(useqpeaks)$motif.no==0)),
#              sum(as.numeric(elementMetadata(useqpeaks)$motif.no==1)),
#                     sum(as.numeric(elementMetadata(useqpeaks)$motif.no>1))),
#         labels=c("0","1","2+"), main=list("USeq",cex=1.5),
#              col=c("white","lightgrey","darkgrey")
#         )
#pie( c(sum(as.numeric(elementMetadata(chippeaks)$motif.no==0)),
#         sum(as.numeric(elementMetadata(chippeaks)$motif.no==1)),
#           sum(as.numeric(elementMetadata(chippeaks)$motif.no>1))),
#    labels=c("0","1","2+"), main=list("chipseq",cex=1.5),
#    col=c("white","lightgrey","darkgrey")
#    )
#
#
## Motif localisation
## set window size
#mydist <- 200
## function to determine the motif profile around peak summits
#getProfile <- function(peaks.GR, pwm, window.size){
#    # get regions around summit
#    summits.GR <- GRanges(
#                                  seqnames=seqnames(peaks.GR),
#                                          range=IRanges( start=elementMetadata(peaks.GR)$maxpos - window.size,
#                                                                        end=elementMetadata(peaks.GR)$maxpos + window.size),
#                          strand="+"
#                          )
#    # create unique names for all peaks
#    names(summits.GR) <- paste( seqnames(summits.GR), start(summits.GR), end(summits.GR), sep=":")
#    # get sequence
#    elementMetadata(summits.GR)$seqs <- getSeq(Hsapiens, summits.GR, as.character=TRUE)
#    # scan sequences with the PWM
#    summit.motifs <- ScanPeaks(summits.GR, pwm, "80%")
#    summit.motifs <- summit.motifs[!is.na(summit.motifs[,1]),]
#    # get all covered positions
#    motif.pos <- sapply(as.integer(summit.motifs[,1]), function(x) seq(x, x+ncol(pwm)-1))
#    motif.pos <- table(unlist(motif.pos))
#    # convert to data.frame
#    motif.pos <- data.frame(motif.pos)
#    names(motif.pos) <- c("position","frequency")
#    # shift positions relative to summit
#    motif.pos$position <- as.integer(as.character(motif.pos$position)) - window.size
#    # ensure all positions are present in the output
#    profile <- data.frame(
#                                  position=-window.size:window.size,
#                                          frequency=motif.pos$frequency[match(-window.size:window.size, motif.pos$position)]
#                                  )
#    profile[is.na(profile)] <- 0
#    return(profile)
#    # get the profiles for the three different peak sets
#    macs.profile <- getProfile(macspeaks, as.matrix(NFKB.pwm), mydist)
#    useq.profile <- getProfile(useqpeaks, as.matrix(NFKB.pwm), mydist)
#    chipseq.profile <- getProfile(chippeaks, as.matrix(NFKB.pwm), mydist)
#    # generate plots
#    pl1 <- xyplot(frequency~position, data=macs.profile,
#                          type="l", main="motifs around MACS v2 peak summits",
#                                  aspect=0.8
#                          )
#    pl2 <- xyplot(frequency~position, data=useq.profile,
#                          type="l", main="motifs around USeq peak summits",
#                                  aspect=0.8
#                          )
#    pl3 <- xyplot(frequency~position, data=chipseq.profile,,
#                          type="l", main="motifs around chipseq peak summits",
#                                  aspect=0.8
#                          )
#    # print plots
#    print(pl1, split=c(1,1,1,3), more=TRUE)
#    print(pl2, split=c(1,2,1,3), more=TRUE)
#    print(pl3, split=c(1,3,1,3))
#
## }}}
