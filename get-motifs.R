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

outputPath <- file.path(args$out, 'meme', basename(file_path_sans_ext(args$peaks)))
dir.create(outputPath, recursive= TRUE)
output <- file.path(outputPath, paste0(args$mode, args$npeaks))

# function in annotation-functions.R
annotationInfo(tolower(args$assembly))

# }}}

# Import some libraries & Turn off warning messages for loading the packages-globally# {{{
options(warn=-1)
suppressMessages(library(BSgenome))
suppressMessages(library("BSgenome.Mmusculus.UCSC.mm9"))
source("~/local/granges.functions.R")
#Turn warnings back on
options(warn=0)# }}}

# Get sequences for peaks # {{{
seq_motifs <- function(summits, n){
    seqs <- getSeq(Mmusculus, summits[1:n], as.character=FALSE)
    # make a dataframe holding information on FE, FDR, sequence
    names(seqs) <- paste(names(summits)[1:n],
                            as.vector(values(summits)[['FE']][1:n]),
                            as.vector(values(summits)[['fdr']][1:n]), sep='|')

    writeXStringSet(seqs, file= paste0(output, '.fa'), "fasta", append=FALSE)
}# }}}

# Read in peaks and get summits# {{{
peaks <- macs2GRanges(args$peaks)
peaks <- peaks[order(-values(peaks)$score)]
peaks <- add.seqlengths(peaks, chr.sizes)

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
names(summits) <- paste(seqnames(summits), start(summits), end(summits), sep=":")# }}}

seq_motifs(summits, args$npeaks)

#n <- 300# {{{
#print(paste0("Getting sequences for top ", n, " peaks"))
#writeXStringSet(seqs[names(seqs) %in% names(summits)[1:n]], file= paste(out, prefix, ".top", n, ".fa", sep=""), "fasta", append=FALSE)
#
#print(paste0("Getting sequences for bottom ", n, " peaks"))
#writeXStringSet(seqs[names(seqs) %in% names(summits)[length(summits)-(n-1):length(summits)]], file= paste(out, prefix, ".bottom", n, ".fa", sep=""), "fasta", append=FALSE)
#
#print(paste0("Getting sequences for random ", n, " peaks"))
#index <- sample(1:length(summits), size=n, replace=FALSE)
#writeXStringSet(seqs[names(seqs) %in% names(summits)[index]], file= paste(out, prefix, ".random", n, ".fa", sep=""), "fasta", append=FALSE)
## }}}

memeBin <- 'meme'
# possibly don't limit to 3 motifs but filter by E
system(sprintf('%s %s -dna -oc %s -maxsize %s -mod zoops -nmotifs 3 -evt 0.05 -minw 6 -maxw 35 -revcomp',
               memeBin,
               paste0(output, '.fa'),
               round(as.numeric(system(paste0("wc -c ", paste0(output, '.fa'), " | awk -F' ' '{print $1}'"), intern=TRUE)) -2),
               file.path(output, 'result')
               ))

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
}

parsePspm(file.path(output, 'result'))

#fa.f <- list.files(gsub("/$", '', out), pattern=paste0(prefix, ".*", n, ".fa"), full.names=TRUE)# {{{
#maxsize <- max(sapply(fa.f, function(f){
#                  char <- system(paste0("wc -c ", f, " | awk -F' ' '{print $1}'"), intern=TRUE)
#                  return(round(as.numeric(char), -2))
#                   }))
#
# nmotifs = max number of motifs to find
# minsites= min number of sites for each motif
# minw = min width of motif, maxw = max n of motifs
# revcomp =  allow sites on + or - strands
# maxsize = minimum dataset size in characters
# -evt = stop if motif E-value greater than <evt>
#meme.param <- paste0(" -nmotifs 3 -minsites 100 -minw 6 -maxw 35 -revcomp -dna -maxsize ", maxsize, " -oc ")
#cluster.param <- 'bsub -M 20000 -R "rusage[mem=10000]"'
#
#for (fa in fa.f){
#   print(paste0("Running meme on ", fa, ". Meme options: ", meme.param))
#   #removing fa extension
#   meme.out <- substr(fa, 1, nchar(fa)-3)
#   meme.command <- paste(cluster.param, " meme ", fa, meme.param, meme.out, sep="")
#   print(meme.command)
#   system(meme.command)
#}
#
## }}}
