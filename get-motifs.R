#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
parser <-  ArgumentParser(description="Perform motif analysis")
parser$add_argument('-x', '--peaks', metavar= "file", required='True', type= "character", help= "MACS generated excel file for called peaks")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
parser$add_argument('-d', '--dist', type= "integer", default= 50,  help= "Distance around the summit to look for motifs. Default: 50bp")
parser$add_argument('-n', '--npeaks', type= "integer", default= 300,  help= "Number of peaks to look in for motifs. Default: 300")

args <- parser$parse_args()

meme <- file.path(args$out, 'meme', '')
plots <- file.path(args$out, 'plots', '')
dir.create(meme)
dir.create(plots)

#args <- commandArgs(trailingOnly=TRUE)
#This has to be the peaks.xls
#peaks.f <- args[1]
#assembly <- args[2]
#out <- args[3]

prefix <- gsub("_.*",'',gsub(".*/",'', args$peaks))

#assembly= 'mm9' # convert to a config that will have assembly\tchromfilepath\tannotationfilepath
if(grepl('mm10', args$assembly, ignore.case=TRUE)){
      chr.sizes <- "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10.chrom.sizes"
}

if(grepl('mm9', args$assembly, ignore.case=TRUE)){
      chr.sizes <- "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/mm9.chrom.sizes"
}# }}}

# Import some libraries & Turn off warning messages for loading the packages-globally# {{{
options(warn=-1)
suppressMessages(library(BSgenome))
suppressMessages(library("BSgenome.Mmusculus.UCSC.mm9"))
source("~/local/granges.functions.R")
#Turn warnings back on
options(warn=0)# }}}

#peaks <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/macs/E12.M2_peaks.xls", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#peaks <- read.table(args$peaks, header=TRUE, sep="\t", stringsAsFactors=FALSE)

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

# Get sequences for peaks
seqMotifs <- function(summits, n){
    tmp <- getSeq(Mmusculus, summits[1:n], as.character=FALSE)
    seqs <- as.data.frame(tmp)
    colnames(seqs) <- 'sequence'
    values(summits[1:n]) <- cbind(values(summits[1:n]), seqs)
    print(summits)
    stop(-1)

    #add seq info to grange
#    s <- as.data.frame(seqs)
#    colnames(s) <- 'seq'
#    values(summits) <- cbind(values(summits), s)
    names(seqs) <- names(summits)
    writeXStringSet(seqs, file= paste(out, prefix, "-top", n, ".fa", sep=""), "fasta", append=FALSE)
}

seqMotifs(summits, args$npeaks)

n <- 300
print(paste0("Getting sequences for top ", n, " peaks"))
writeXStringSet(seqs[names(seqs) %in% names(summits)[1:n]], file= paste(out, prefix, ".top", n, ".fa", sep=""), "fasta", append=FALSE)

print(paste0("Getting sequences for bottom ", n, " peaks"))
writeXStringSet(seqs[names(seqs) %in% names(summits)[length(summits)-(n-1):length(summits)]], file= paste(out, prefix, ".bottom", n, ".fa", sep=""), "fasta", append=FALSE)

print(paste0("Getting sequences for random ", n, " peaks"))
index <- sample(1:length(summits), size=n, replace=FALSE)
writeXStringSet(seqs[names(seqs) %in% names(summits)[index]], file= paste(out, prefix, ".random", n, ".fa", sep=""), "fasta", append=FALSE)

fa.f <- list.files(gsub("/$", '', out), pattern=paste0(prefix, ".*", n, ".fa"), full.names=TRUE)
maxsize <- max(sapply(fa.f, function(f){
                  char <- system(paste0("wc -c ", f, " | awk -F' ' '{print $1}'"), intern=TRUE)
                  return(round(as.numeric(char), -2))
                   }))

# nmotifs = max number of motifs to find
# minsites= min number of sites for each motif
# minw = min width of motif, maxw = max n of motifs
# revcomp =  allow sites on + or - strands
# maxsize = minimum dataset size in characters
# -evt = stop if motif E-value greater than <evt>
meme.param <- paste0(" -nmotifs 3 -minsites 100 -minw 12 -maxw 35 -revcomp -dna -maxsize ", maxsize, " -oc ")
cluster.param <- 'bsub -M 20000 -R "rusage[mem=10000]"'

for (fa in fa.f){
   print(paste0("Running meme on ", fa, ". Meme options: ", meme.param))
   #removing fa extension
   meme.out <- substr(fa, 1, nchar(fa)-3)
   meme.command <- paste(cluster.param, " meme ", fa, meme.param, meme.out, sep="")
   print(meme.command)
   system(meme.command)
}


