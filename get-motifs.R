#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")

parser =  ArgumentParser(description="Perform motif analysis")
parser$add_argument('-x', '--peaks', metavar= "file", required='True', type= "character", help= "MACS generated excel file for called peaks")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
parser$add_argument('-d', '--dist', type= "integer", default= 50,  help= "Distance around the summit to look for motifs. Default: 50bp")
parser$add_argument('-n', '--npeaks', type= "integer", default= 300,  help= "Number of peaks to look in for motifs. Default: 300")
parser$add_argument('-m', '--mode', type= "character", default= 'top',  help= "Look at the top/random/bottom peaks. Default: top")

args = parser$parse_args()

if(FALSE){
    args = list()
    args$peaks = '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/macs2/sharp/7E12_EpiSC-M2_filtered_peaks.narrowPeak'
    args$assembly = 'mm10'
    args$out = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013"
    args$dist = 50
    args$npeaks = 300
    args$mode = 'top'
}

parent_path = file.path(args$out, args$assembly, 'meme', basename(file_path_sans_ext(args$peaks)))
output_path = file.path(parent_path, paste0(args$mode, args$npeaks))
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

# define variables #{{{
genome_path = '/nfs/research2/bertone/user/mxenoph/common/genome'
meme_db_path = '/nfs/research2/bertone/user/mxenoph/common/meme/motif_databases'

# Read and assign to global the chromosme lengths
get_chr_length = function(){
    chr_file = paste0(file.path(genome_path, toupper(args$assembly), toupper(args$assembly)),
                      '.genome')
    if(file.exists(chr_file)){
        chr_lengths = read.table(chr_file, row.names = 1)
        colnames(chr_lengths) = 'size'
        assign('chr_lengths', chr_lengths, envir = globalenv())
    } else {
        assert(paste0(chr_file, " file does not exist."), !file.exists(chr_file))
    }

}
get_chr_length()

# }}}

# }}}

# Import libraries & Turn off warning messages for loading the packages-globally# {{{
suppress = base::suppressPackageStartupMessages
options(warn=-1)
source("~/source/Rscripts/granges-functions.R")
#Turn warnings back on
options(warn=0)
# }}}

# Functions # {{{

# Functions by konrad Rudolph# {{{
# Used like that in Konrad's parser of pspm
map = base::Map

# This uses R's peculiarities in argument matching explained here:
# <http://stat.ethz.ch/R-manual/R-devel/doc/manual/R-lang.html#Argument-matching>
# `.expr` starts with a dot to allow `expr` being used in the actual
# expression.
let = function (.expr, ...)
        eval(substitute(.expr), list2env(list(...), parent = parent.frame()))

#' Create a closure over a given environment for the specified formals and body.
closure = function (formals, body, env)
        eval(call('function', as.pairlist(formals), body), env)
#' Create a list of empty symbols, with names set
symlist = function (names)
        setNames(Map(function (p) quote(expr = ), names), names)

#' A shortcut to create a function
#'
#' @note Using \code{.(args = body)} is analogous to using
#' \code{function (args) body} with one exception: \code{.} arguments do not
#' support defaults.
#' Since its purpose is mainly for lambdas in higher-order list functions, this
#' functionality is not needed.
. = function (...) {
    args = match.call(expand.dots = FALSE)$...
    last = length(args)
    params = symlist(c(args[-last], names(args)[[last]]))
    if (length(args) > 1 && length(params) != length(args))
        stop('Must be of the form `fun(a, b = expr)`')
    for (arg in args[-last])
        if (! is.name(arg))
            stop('Invalid argument specifier: ', arg)

    closure(params, args[[last]], parent.frame())
}# }}}


# Read in peaks and get summits# {{{
getLocus = function(file){

    if(grepl('.narrowPeak', as.character(file))) {
        peaks = import_narrowPeak(as.character(file))
    } else if(grepl('.bed', as.character(file))){
        peaks = import_bed(as.character(x))
    } else{
        stop("Can not handle this file format.")
    }

    peaks = peaks[order(-values(peaks)$qvalue)]

    # When calling narrow peaks with subpeaks, peaks i reported twice
    # Keep only the highest summit within the peak
    peaks = as.data.frame(peaks) %>% group_by(seqnames, start, end) %>%
        # making sure that the rows are ordered high to low for qvalue and fe
        arrange(desc(qvalue), desc(fe)) %>%
        top_n(1) %>% ungroup()
    peaks = with(peaks, GRanges(seqnames, IRanges(start,end), strand, qvalue=qvalue, fe=fe, summit=summit))
    
    peaks = add_seqlengths_to_gr(peaks, chr_lengths)

    #why just looking at + strand?
    summits = GRanges(seqnames= seqnames(peaks),
                       range= IRanges(start= values(peaks)$summit - args$dist,
                                      end= values(peaks)$summit + args$dist),
                       strand= "+",
                       seqlengths= seqlengths(peaks))

    values(summits) = values(peaks)
    #trimming ranges that exceed the chromosome length
    summits = trim(summits)

    #Because of trimming I might end up with ranges of 0 length so get rid of those
    summits = summits[width(summits) != 0]

    summits = summits[order(-values(summits)$qvalue)]
    names(summits) = paste(seqnames(summits), start(summits), end(summits), sep=":")

    names(peaks) = paste(seqnames(peaks),
                          start(peaks),
                          end(peaks),
                          sep=":")
#    return(summits)
    return(list('peaks'= peaks, 'summits'= summits))
} # }}}

# Get sequences for peaks # {{{
seq_motifs = function(regions, prefix, output_path, assembly=tolower(args$assembly), report=F){
    # Load BSgenome based on assembly
    suppress(library(BSgenome))
    suppress(library(paste0("BSgenome.Mmusculus.UCSC.", assembly), character.only=TRUE))

    seqs = getSeq(Mmusculus, regions, as.character=FALSE)
    # keep information on FE, FDR, sequence
    names(seqs) = paste(names(regions),
                            as.vector(values(regions)[['fe']]),
                            as.vector(values(regions)[['qvalue']]), sep='|')

    fasta = file.path(output_path, paste0(prefix, '.fa'))
    writeXStringSet(seqs, filepath = fasta, format = "fasta", append = FALSE)
    if (report == TRUE) return(as.character(seqs))
}# }}}

# Parse motifs and E-values from MEME output. From Konrad Rudolph  # {{{
parsePspm = function (results) {
    # from Konrad what the hell is the output?
    # Parse output raw text file and retrieve PSSMs for all motifs.
    # Note: we could also parse the XML file using XPath but the XML output file
    # is quite frankly not very nice. Parsing the raw text is easier. Fail.
    extractHeader = function (name, header)
        let(match = regexpr(sprintf('(?<=%s= )\\S+', name), header, perl = TRUE),
            as.numeric(regmatches(header, match)))

    file = file.path(results, 'meme.txt')
    lines = readLines(file)
    start = grep('^\tMotif \\d+ position-specific probability matrix', lines) + 2
    length = extractHeader('w', lines[start])
    nsites = extractHeader('nsites', lines[start])
    evalue = extractHeader('E', lines[start])
    matrix = map(.(start, end = read.table(text = paste(lines[start : end],
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
memeFileContent = function (pspm, fasta) {
    # Format documented at
    # <http://meme.nbcr.net/meme/doc/meme-format.html#min_format>
    version = 'MEME version 4.9\n'

    inputFasta = readLines(fasta)

    # Remove Fasta headers and merge lines
    inputFasta = paste(inputFasta[! grepl('^>', inputFasta)], collapse = '')
    # Normalise casing and remove Ns
    inputFasta = gsub('N', '', toupper(inputFasta))
    freqs = table(strsplit(inputFasta, ''))

    frequencies = sprintf('Background letter frequencies:\n%s\n',
    paste(mapply(paste, names(freqs), freqs), collapse = ' '))

    motif = capture.output(write.table(pspm$matrix, col.names = FALSE,
    row.names = FALSE))

    c(version,
      # Alphabet omitted (can be inferred)
      # Strand omitted (can be inferred)
      frequencies,
      pspmHeader(pspm),
      motif)
}# }}}

consensusSequence = function (pspm){# {{{
    with(pspm, paste(colnames(matrix)[
                     vapply(as.data.frame(t(matrix)), which.max,
                            numeric(1))], collapse = ''))
}# }}}

pspmHeader = function (pspm){# {{{
    paste(sprintf('MOTIF %s', consensusSequence(pspm)),
          sprintf('letter-probability matrix: alength= %s w= %s nsites= %s E= %e',
                  ncol(pspm$matrix), nrow(pspm$matrix), pspm$nsites, pspm$e),
          sep = '\n\n')
}# }}}

plot_motif = function(meme_output){# {{{
    
    get_pfm = function(x){
        library(ggplot2)
        library(ggseqlogo)
        to_skip = as.numeric(system(paste0("awk '/letter-probability/ {print NR}' <",
                                           x), intern = T))
        ppm = read.table(x, header=F, skip=last(to_skip))
        info_line = readLines(x)[to_skip]
        nsites = as.integer(gsub(' E.*', '', gsub('.*nsites= ', '', info_line)))
        pfm = round(t(ppm)*nsites)
        rownames(pfm) = c('A', 'C', 'G', 'T')

        pdf(file.path(dirname(x), 'plots', paste0(file_path_sans_ext(basename(x)), '.pdf')))
        p = ggplot() + geom_logo(pfm, seq_type = 'dna') + theme_logo(base_size = 18, base_family='Helvetica')
        p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        plot(p)
        dev.off()
        return(pfm)
    }
    meme_output = list.files(meme_output, '*.meme', full.name = T)
    names(meme_output) = basename(file_path_sans_ext(meme_output))
    lapply(meme_output, get_pfm)
}
# }}}

# Run MEME# {{{
meme = function(output) {
    memeBin = 'meme'
    # nmotifs = max number of motifs to find
    # minsites= min number of sites for each motif
    # minw = min width of motif, maxw = max n of motifs
    # revcomp =  allow sites on + or - strands
    # maxsize = minimum dataset size in characters
    # -evt = stop if motif E-value greater than <evt>
    # possibly don't limit to 3 motifs but filter by E-value
    # -p ncores to use for parallel processing
    fasta = list.files(output, pattern = ".fa", full.name = TRUE)

    system(sprintf('%s %s -dna -oc %s -maxsize %s -mod zoops -nmotifs 3 -evt 0.05 -minw 6 -maxw 35 -revcomp',
                   memeBin,
                   fasta,
                   file.path(output),
                   round(as.numeric(system(paste0("wc -c ", fasta, " | awk -F' ' '{print $1}'"), intern=TRUE)) -2)
                   ))

    parsePspm(file.path(output))
} # }}}

# Run tomtom# {{{
tomtom = function (pspm, outputPath, databases) {
    motif = consensusSequence(pspm)
    inputFile = file.path(outputPath, paste0(motif, '.meme'))
    fasta = list.files(pattern = "*.fa", output_path, full.names = TRUE)
    writeLines(memeFileContent(pspm, fasta), inputFile)

    if (missing(databases))
        databases = file.path(meme_db_path,
                               c('JASPAR_CORE_2014_vertebrates.meme',
                                 'uniprobe_mouse.meme'))
   
    system(sprintf('%s -no-ssc -oc %s -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 0.1 %s %s',
                   'tomtom',
                   file.path(outputPath, motif),
                   inputFile,
                   databases))

    extractMotifName = function (header){
        match = regexpr(sprintf('\\t(\\S+)\\t'), header, perl = TRUE)
        id = gsub('\t', '', regmatches(header, match))
        name = system(sprintf('grep %s %s | cut -d \' \' -f 3 | tr -d \'\n\'',
                               id,
                               file.path(meme_db_path, c('JASPAR_CORE_2014_vertebrates.meme',
                                                             'uniprobe_mouse.meme'))
                               ), intern = TRUE)
        return(name)
    }
   lines = readLines(file.path(outputPath, motif, 'tomtom.txt'), n= -1L)
   if ( length(lines) != 1 ){
       lines = lines[2:length(lines)]
       factors = as.vector(unlist(map(extractMotifName, lines)))
       combined = as.vector(unlist(map(function(x, y)
                                        paste(x, y, sep = '\t'),
                                        lines, factors)))

       write.table(file = file.path(outputPath, motif, 'tomtom.txt'), combined,
                   col.names = paste(c('QueryID', 'TargetID', 'Offset', 'pvalue',
                                       'Evalue', 'qvalue', 'Overlap', 'QueryConsensus',
                                       'TargetConsensus', 'Orientation', 'Factor'),
                                     collapse = '\t'),
                   row.names = FALSE,
                   quote = FALSE)

       motifId = read.table(pipe(sprintf('cut -f 2 %s',
                                         file.path(outputPath, motif, 'tomtom.txt'))),
                             header = TRUE)
   } else {
       motifId = c(NA)
       factors = c(NA)
   }
   return(list(motif= motif, motif_id= motifId, motif_name = factors))
}# }}}

# Get sequences for peaks # {{{
seq_general = function(regions){
    seqs = getSeq(Mmusculus, regions, as.character=TRUE)
    # make a dataframe holding information on FE, FDR, sequence
    names(seqs) = names(regions)

#    fasta = file.path(output_path, paste0('top', n, '.fa'))
#    writeXStringSet(seqs, file= fasta, "fasta", append=FALSE)
    return(seqs)
}# }}}

# Parse motifs and E-values from MEME formatted databases. Adapted from Konrad Rudolph  # {{{
parsePspmDb = function (databasePath) {
    # Parse output raw text file and retrieve PSSMs for all motifs.
    # Note: we could also parse the XML file using XPath but the XML output file
    # is quite frankly not very nice. Parsing the raw text is easier. Fail.
    extractHeader = function (name, header)
        let(match = regexpr(sprintf('(?<=%s= )\\S+', name), header, perl = TRUE),
            as.numeric(regmatches(header, match)))
    extractName = function (header)
        gsub('MOTIF ', '', let(match = regexpr('^MOTIF \\S+', header, perl = TRUE),
                              as.character(regmatches(header, match))))

    lines = readLines(databasePath)
    start = grep('^MOTIF \\S+', lines) + 2
    name = extractName(lines[start-2])
    length = extractHeader('w', lines[start])
    evalue = extractHeader('E', lines[start])
    matrix = map(.(start, end = read.table(text = paste(lines[start : end],
                                                         collapse = '\n'),
                                            col.names = c('A', 'C', 'G', 'T'))),
                  start + 1, start + length)

    pspm = map(.(matrix, evalue =
                  let(str = list(name = name, matrix = matrix, e = evalue),
                      structure(str, class = 'pspm'))),
                matrix, evalue)
    names(pspm) = name
    return(pspm)
}# }}}

# Scan a string/sequence for PWM matches# {{{
extractPspmMatches = function(pspm, sequence, threshold, summit){
    # Get all matches for pspm
    # needs a numeric matrix with rows A,C,T,G, the minimum score for counting a match
    # that could be a percentage of the highest possible score or a number
   match = matchPWM(pspm, sequence, min.score = threshold)
   rmatch = matchPWM(reverseComplement(pspm), sequence, min.score = threshold)

   parseMatches = function(match){
       if (length(match) == 0)
           found = cbind(NA, NA, 0)
       else
           found = cbind(start(match),
                          as.character(match),
                          length(start(match)))
       colnames(found) = c('start', 'seq', 'hits')
       return(found)
   }
   fwd = parseMatches(match)
   rev = parseMatches(rmatch)
   found = rbind(fwd, rev)
   
   # Keep only the match that is closest to the summit, that's why we use abs()
   found = found[order(abs(as.integer(found[,'start']) - summit), decreasing = FALSE)[1], ]
   return(found)
}# }}}

# Scan peaks for a motif # {{{
scanPeaks = function(peaks, pspm, threshold){
    # Get sequences for all peaks
    sequences = elementMetadata(peaks)$seqs
    motifs = sapply(1:length(sequences),
                     function(indx){
                         extractPspmMatches(pspm, sequences[indx], threshold, elementMetadata(peaks)$summit[indx])})
    motifs = t(motifs)
    rownames(motifs) = names(peaks)
    return(motifs)
}# }}}

# Get the motif profile around the peak summits# {{{
getPeakProfile = function(peaks, pspm, window){
    summits = GRanges(seqnames = seqnames(peaks),
                       range = IRanges(start = elementMetadata(peaks)[['summit']] - window,
                                       end = elementMetadata(peaks)[['summit']]) + window,
                       strand = '+',
                       summit = elementMetadata(peaks)[['summit']])
    names(summits) = paste(seqnames(summits), start(summits), end(summits), sep=':')

    elementMetadata(summits)$seqs = seq_general(summits)
    scanned = scanPeaks(summits, pspm, "80%")
    # Keep only those ranges for which a motif was matched
    scanned = scanned[!is.na(scanned[,'start']), ]
    # Get position covered by the motif
    positions = sapply(as.integer(scanned[,'start']), function(x) x + ncol(pspm)-1)
    # Frequency on the plot is how many peaks have a motif at that position
    positions = data.frame(table(unlist(positions)))
    names(positions) = c('position','frequency')
    # shift position relative to summit
    positions$position = as.integer(as.character(positions$position)) - window
    # ensure all positions are present in the output
    profile = data.frame('position' = c(-window:window),
                          'frequency' = positions$frequency[match(-window:window, positions$position)])
    profile[is.na(profile)] = 0
    return(profile)
}# }}}

runPeakProfileOnAll = function(peaks, motifs, output_path){# {{{
    pspmDb = parsePspmDb(file.path(meme_db_path, "JASPAR_CORE_2014_vertebrates.meme"))
    pspmDb = c(pspmDb, parsePspmDb(file.path(meme_db_path, "uniprobe_mouse.meme")))

    library(gridExtra)

    getResults = function(peaks, motif_id, db){
        pspm = t(db[[motif_id]]$matrix)
        scanned = scanPeaks(peaks, pspm, "80")
        profile = getPeakProfile(peaks, pspm, 600)

        relative_summit = peaks$summit - start(peaks)
        positions = as.integer(scanned[,'start']) - relative_summit
        positions = positions + ncol(pspm)-1
#
        result = list(hits = scanned, profile = profile, positions = positions)
        return(result)
    }

    plotResult = function(results){# {{{
        p = vector("list", length(results)*2)
        for (i in seq(1, length(results)*2, 2)) {
            x = ceiling(i/2)
            tmp = as.data.frame(results[[x]][['hits']])
            # Convert factor to numeric values so that x-axis is sorted
            tmp$hits = as.numeric(as.character(tmp$hits))
            a = ggplot(tmp, aes(hits))
            a = a + geom_bar(stat = "bin", binwidth = 1) + geom_hline(yintercept=0, colour="white", size=0.5)
            a = a + labs(title = paste0('Motif', indx, ':', names(results)[x]))

            b = ggplot(as.data.frame(results[[x]][['profile']]), aes(x=position, y=frequency))
            b = b + geom_point(aes(colour = cut(position, breaks=c(-Inf,-50,50, Inf)),
                                    alpha = frequency))
            b = b + guides(colour = FALSE) + theme(legend.position='none')
            p[[i]] = a
            p[[i+1]] = b

        }
        # call marrangeGrob instead of grid.arrange so that the plots span
        # multiple pdf pages as in answer in 
        # <http://stackoverflow.com/questions/19059826/multiple-graphs-over-multiple-pages-using-ggplot>
        # top=NULL removes page numbers
        multipage = do.call(marrangeGrob, c(p, nrow = 3, ncol = 2, top = NULL))
        i = unlist(strsplit(names(results), ":"))[1]
        ggsave(file.path(plot_path, paste0('summary-', i, '.pdf')), multipage, width = 8.3, height = 11.7)
    }# }}}

    runAllMotifs = function(motifs, peaks, db){
        # Get all Jaspar/Uniprobe IDs matched to each predicted motif
        results = mclapply(motifs, function(m){
                        ids = as.character(unlist(m[['motif_id']]))
                        if(any(is.na(ids))){
                            return()
                        } else {
                            results = mclapply(ids, getResults, peaks=peaks, db=db)
                            return(results)
                        }
                          })

        mnames = mclapply(motifs, function(m)
                        m[['motif_name']])



    }
#    pspms = sapply(1:length(motifs), function(indx, db){# {{{
#    pspms = lapply(1:length(motifs), function(indx, db){
#                    ids = as.character(unlist(motifs[[indx]]$motif_id))
#                    mnames = motifs[[indx]]$motif_name
#
#                    # Only run if tomtom found TFs matching the motif at FDR 0.1
#                    if ( !any(is.na(ids)) ) {
#                        results = lapply(ids, function(x){
#                                          pspm = t(db[[x]]$matrix)
#                                          scanned = scanPeaks(peaks, pspm, "80")
#                                          profile = getPeakProfile(peaks, pspm, 600)
#
#                                          relative_summit = peaks$summit - start(peaks)
#                                          positions = as.integer(scanned[,'start']) - relative_summit
#                                          positions = positions + ncol(pspm)-1
#
#                                          result = list(hits = scanned, profile = profile, positions = positions)
#                                          return(result)
#                                           })
#                        names(results) = paste(ids, mnames, sep=':')
#
#                        p = vector("list", length(results)*2)
#                        for (i in seq(1, length(results)*2, 2)) {
#                            x = ceiling(i/2)
#                            tmp = as.data.frame(results[[x]][['hits']])
#                            # Convert factor to numeric values so that x-axis is sorted
#                            tmp$hits = as.numeric(as.character(tmp$hits))
#                            a = ggplot(tmp, aes(hits))
#                            a = a + geom_bar(stat = "bin", binwidth = 1) + geom_hline(yintercept=0, colour="white", size=0.5)
#                            a = a + labs(title = paste0('Motif', indx, ':', names(results)[x]))
#
#                            b = ggplot(as.data.frame(results[[x]][['profile']]), aes(x=position, y=frequency))
#                            b = b + geom_point(aes(colour = cut(position, breaks=c(-Inf,-50,50, Inf)),
#                                                    alpha = frequency))
#                            b = b + guides(colour = FALSE) + theme(legend.position='none')
#                            p[[i]] = a
#                            p[[i+1]] = b
#
#                        }
#                        # call marrangeGrob instead of grid.arrange so that the plots span
#                        # multiple pdf pages as in answer in 
#                        # <http://stackoverflow.com/questions/19059826/multiple-graphs-over-multiple-pages-using-ggplot>
#                        # top=NULL removes page numbers
#                        multipage = do.call(marrangeGrob, c(p, nrow = 3, ncol = 2, top = NULL))
#                        ggsave(file.path(plot_path, paste0('summary-motif-new', indx, '.pdf')), multipage, width = 8.3, height = 11.7)
#
#                        return(results)

    #                    print('in motif')
    #                    print(indx)
#                    } else {
#                        return(c(NA))
#                    }
#                          }, pspmDb)
#
#    names(pspms) = paste0('Motif', 1:length(motifs))
#    return(pspms)# }}}
}# }}}

# end of fold sections # }}}

main = function(){
    loci = getLocus(args$peaks)
    seq_motifs(loci$summits[1:300], paste0(args$mode, args$npeaks), output_path)

    if (!exists(paste0(output_path, '.Rda'))) {
        print('Run meme')
        pspm = meme(output_path)
        motifs = map(tomtom, pspm, output_path)
    } else {
        load(file=paste0(output_path, '.Rda'))
        plot_motif(output_path)
    }

    save(loci, pspm, motifs, file=paste0(output_path, '.Rda'))
    elementMetadata(loci$peaks)$seqs = seq_motifs(loci$peaks, paste0('all-peaks'), output_path, report= TRUE)
}

main()

#peaks = macs_to_granges(args$peaks)
#peaks = add.seqlengths(peaks, chr_size)
#peaks = peaks[order(-values(peaks)$score)]
#names(peaks) = paste(seqnames(peaks), start(peaks), end(peaks), sep=":")
#elementMetadata(peaks)$seqs = seq_motifs(peaks, paste0('all-peaks'), output_path, report= TRUE)
##
#loci = getLocus(args$peaks)
#seq_motifs(loci, args$npeaks, output_path)
#
#if (!exists(file.path(paste0(output_path, '.Rda')))) {
#    pspm = meme(output_path)
#    motifs = map(tomtom, pspm, output_path)
#} else {
#    load(file=paste0(output_path, '.Rda'))
#}
#
#tfs = runPeakProfileOnAll(peaks, motifs, output_path)


# TODO: seq_motifs on random and bottom peaks {{{
# should be a new separate function
# print(paste0("Getting sequences for bottom ", n, " peaks"))
# writeXStringSet(seqs[names(seqs) %in% names(summits)[length(summits)-(n-1):length(summits)]], file= paste(out, prefix, ".bottom", n, ".fa", sep=""), "fasta", append=FALSE)
#
# print(paste0("Getting sequences for random ", n, " peaks"))
# index = sample(1:length(summits), size=n, replace=FALSE)
# writeXStringSet(seqs[names(seqs) %in% names(summits)[index]], file= paste(out, prefix, ".random", n, ".fa", sep=""), "fasta", append=FALSE)
## }}}

## ggplot specific# {{{
#require(grid)
#vp.layout = function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
#arrange_ggplot2 = function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
#    dots = list(...)
#    n = length(dots)
#    if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
#    if(is.null(nrow)) { nrow = ceiling(n/ncol)}
#    if(is.null(ncol)) { ncol = ceiling(n/nrow)}
#    ## NOTE see n2mfrow in grDevices for possible alternative
#    grid.newpage()
#    pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
#    ii.p = 1
#    for(ii.row in seq(1, nrow)){
#        ii.table.row = ii.row  
#        if(as.table) {ii.table.row = nrow - ii.table.row + 1}
#        for(ii.col in seq(1, ncol)){
#            ii.table = ii.p
#            if(ii.p > n) break
#            print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
#            ii.p = ii.p + 1
#        }
#    }
#}# }}}
