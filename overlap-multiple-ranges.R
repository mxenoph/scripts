#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Get exons, introns, intergenic for given annotation")
parser$add_argument('-b', '--bed', metavar= "file", nargs = "+", required='True', type= "character", help= "BED files, narrowPeak, broadPeak etc")
parser$add_argument('-i', '--overlap', type= "character", default = -1,
                    help= "Keep regions that overlap by this percentage. If not provided default is -1, to use 1bp overlap.")
parser$add_argument('-p', '--threads', default = 4, help= "Number of processors to use")
parser$add_argument('-o', '--output_path', required='True', help= "Output path")
parser$add_argument('-l', '--label', default = "NA", help= "Useful when comparing complex in conditions. e.g 2i-NuRD:EpiSC-NuRD")
parser$add_argument('-a', '--anno', default = "NA", help= "Rename the merged ranges to anno in output bed")

args = parser$parse_args()

if(FALSE){
    args = list()
    args$bed = c("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/macs2/sharp/7E12_EpiSC-M2_filtered_peaks.narrowPeak",
                 "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/macs2/sharp/7E12_EpiSC-Chd4_filtered_peaks.narrowPeak")
    args$threads = 4
    args$output_path = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/merged-ranges"
    args$label = 'NA'
    args$anno = 'NA'
    args$overlap = -1
}

plot_path = file.path(args$output_path, "plots")
dir.create(plot_path, recursive= TRUE)

ncores = args$threads
#}}}

import_narrowPeak_df = function(narrow) {# {{{
    x = c('dplyr','GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)
    if (grepl('track', readLines(narrow, n=1))){
        # File contains track line for ucsc which we skip
        df = read.table(narrow, header=F, skip=1)
    } else {
        df = read.table(narrow, header=F)
    }
    names(df) = c('chr', 'start', 'end', 'name', 'score', 'strand', 'fe', 'pvalue', 'qvalue', 'summit')
    return(df)
}

import_broadPeak_df = function(broad){# {{{
    x = c('dplyr','GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)
    
    if (grepl('track', readLines(broad, n=1))){
        to_skip = as.numeric(system(paste0("awk '/track/ {print NR}' <", broad), intern = T))
        # check if vector is sequential and increasing!
        if(! all(diff(to_skip) == 1)) stop("File contains more than one track line. Lines not sequential so don't know how to skip them. Inspect the file. ")

        # File contains track line for ucsc which we skip
        df = read.table(broad, header=F, skip=last(to_skip))
    } else {
        # contains both the broad region and narrow peak
        df = read.table(broad, header=F)
    }
    
    names(df) = c('chr', 'start', 'end', 'name', 'score', 'strand', 'fe', 'pvalue', 'qvalue')
    return(df)
}
# }}}

import_bed_df = function(bed) {# {{{
    x = c('dplyr','GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)
    if (grepl('track', readLines(bed, n=1))){
        # File contains track line for ucsc which we skip
        df = read.table(bed, header=F, skip=1)
    } else {
        df = read.table(bed, header=F)
    }
    names(df) = c('chr', 'start', 'end', 'name', 'score', 'strand')
    return(df)
}
#}}}
#}}}

# Extract clusters from Hits object.# {{{
# from https://support.bioconductor.org/p/68021/
extract_clusters_from_self_hits = function(hits){
    library(GenomicRanges)
    stopifnot(is(hits, "Hits"))
    stopifnot(queryLength(hits) == subjectLength(hits))
    # not quite sure why he is doing this
    hits = BiocGenerics::union(hits, t(hits))
    qh = queryHits(hits)
    sh = subjectHits(hits)
    # Assume all ranges are a cluster of there own to start with
    cid = seq_len(queryLength(hits))  # cluster ids
    cnt = 0
    while (TRUE) {
        # create a hits object. This changes every time as cid changes
        # until it's the same with cid2
        h = Hits(qh, cid[sh],
                 queryLength(hits), subjectLength(hits))
        cid2 = pmin(cid, selectHits(h, "first"))
        if (identical(cid2, cid))
            break
        cid = cid2
    }
    # always returns an integer list where the groups are defined in seq_len(queryLength(hits)
    # and the values for each group come from cid
    unname(splitAsList(seq_len(queryLength(hits)), cid))
}
# }}}

# Merge ranges that are "connected" (directly or indirectly)# {{{
# via a hit (or several hits) in 'hits'.
merge_connected_ranges = function(x, hits){
    library(GenomicRanges)
    stopifnot(is(x, "GenomicRanges"))
    stopifnot(is(hits, "Hits"))
    stopifnot(queryLength(hits) == subjectLength(hits))
    stopifnot(queryLength(hits) == length(x))
    clusters = extract_clusters_from_self_hits(hits)
    # will merge all ranges in the cluster e.g. chr1 [2, 8] and chr1 [3,12] will 
    # become chr1 [2,12]
    ans = range(extractList(x, clusters))
    if (any(elementLengths(ans) != 1L))
        stop(wmsg("some connected ranges are not on the same ",
                  "chromosome and strand, and thus cannot be ",
                  "merged"))
    ans = unlist(ans)
    mcols(ans)$revmap = clusters
    ans
}
# }}}

# Input is a grange that contains all ranges for TFs to be compared. Metadata are important;# {{{
# the name of the TF must be set in gr0
get_merged = function(input_granges, min_ov = 0){
    library(GenomicRanges)
    hits = findOverlaps(input_granges)
    
    ## subset 'hits' to keep only hits that achieve 60% overlap.
    x = input_granges[queryHits(hits)]
    y = input_granges[subjectHits(hits)]
    relative_overlap = width(pintersect(x, y)) / pmin(width(x), width(y))
    hits = hits[relative_overlap > min_ov]
    
    ## merge the ranges in 'input_granges' that are connected via one or more hits in 'hits'.
    merged_ranges = merge_connected_ranges(input_granges, hits)
    groups = CharacterList(mclapply(mcols(merged_ranges)$revmap,
                                    function(x){
                                        x = mixedsort(unique(values(input_granges[x])[['protein']]))
                                        if (length(x) != 1) {
                                            x = paste(x, collapse = "_AND_")
                                        }
                                        return(x)
                                    }, mc.cores = 4))
    mcols(merged_ranges)$name = groups
    scores = unlist(mclapply(mcols(merged_ranges)$revmap,
                                    function(x){
                                        x = max(values(input_granges[x])[['score']])
                                        return(x)
                                    }, mc.cores = 4))
    mcols(merged_ranges)$score = scores

    # if group name contains _and_ then range is merged from connected tf ranges
    multiple = grepl("*_AND_*", groups)
    multiple_peaks = merged_ranges[multiple]
    single_peaks = merged_ranges[!multiple]
    return(GRangesList('annotated' = merged_ranges, 'merged' = multiple_peaks, 'single' = single_peaks))
} # }}}

# Make venn count table# {{{
make_venn_cnt = function(granges){
    stopifnot(is(granges, "GenomicRanges"))
    stopifnot(any(grepl('name', colnames(mcols(granges)))))

    sets = levels(as.factor(unlist(mcols(granges)$name)))
    columns = sets[grep('_AND_', sets, invert=T)]
    x = lapply(1:length(columns), function(x) c(0,1))
    names(x) = columns
    mat = expand.grid(x)
    # expand.grid always puts the (0 0 0) case first
    Counts = c(0)
    for(i in 2:nrow(mat)){
        current_TF = columns[as.logical(unlist(mat[i,]))]
        if(length(current_TF) == 1) {
            regex = paste0("^", current_TF, "$")
            counts_table = table(grepl(regex, mcols(granges)$name))
            # if no single peaks for that TF, table will not have a TRUE element
            # so set the count to 0 like this
            if (! TRUE %in% names(counts_table)) {
                Counts = c(Counts, 0)
            } else {
                Counts = c(Counts, counts_table[['TRUE']])
            }
        } else {
            tmp = paste(mixedsort(current_TF), collapse = "_AND_")
            regex = paste0("^", tmp, "$")
            counts_table = table(grepl(regex, mcols(granges)$name))
            if (! TRUE %in% names(counts_table)) {
                Counts = c(Counts, 0)
            } else {
                Counts = c(Counts, counts_table[['TRUE']])
            }
        }
    }
    mat = cbind(mat, Counts)
    return(mat)
}# }}}

main = function(){# {{{
    source("/homes/mxenoph/source/Rscripts/granges-functions.R")
    source("/homes/mxenoph/source/Rscripts/plotting-functions.R")
    library(GenomicRanges)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(gtools)
    library(rtracklayer)

    # Prepare the data# {{{
    if(length(args$bed) < 2) stop("Two or more sets of regions need to be provided.")

    protein_data = as.data.frame(args$bed)
    colnames(protein_data) = 'protein_data'

    if(args$label == "NA"){
        tmp = protein_data %>% .[[1]] %>% as.character() %>%
                basename() %>% file_path_sans_ext %>% str_replace_all("_peaks", "") %>% strsplit("-")
        tmp = do.call(rbind.data.frame, tmp)
        colnames(tmp) = c('condition', 'protein')
        protein_data = cbind(protein_data, tmp)
    } else {
        if(!grepl(':', args$label)) stop("Labels must be ordered as the input files and separated by :")
        tmp = data.frame(tmp = unlist(strsplit(args$label, ":"))) %>% separate(tmp, into = c('condition', 'protein'), sep='-')
        protein_data = cbind(protein_data, tmp)
    }

    proteins_df = lapply(protein_data[['protein_data']], function(x) {
                    tmp = protein_data %>% filter(protein_data == x) %>% unite(id, condition, protein, sep="-") %>% select(id) %>% .[[1]]
                    if(grepl('.bed', as.character(x))){
                        x = import_bed_df(as.character(x))
                    } else if(grepl('.narrowPeak', as.character(x))) {
                        x = import_narrowPeak_df(as.character(x))
                    } else{
                        x = import_broadPeak_df(as.character(x))
                    }
                    x = x[,1:6]
                    tmp = rep(tmp, nrow(x))
                    # Handles peaks that are reported multiple times due to multiple summits
                    x %>% mutate(protein = tmp) %>% group_by(chr, start, end) %>% slice(1:1) %>% ungroup })
# }}}

    # keep all peaks in one GRange object and find all overlaps# {{{
    gr0 = do.call(rbind, proteins_df)
    print(head(gr0))
    if(! 'summit' %in% colnames(gr0)){
        gr0 = with(gr0 %>% mutate(strand=gsub('.', '*', strand)),
                   GRanges(chr, IRanges(start,end), strand, score=score, protein=protein))
    } else {
        # make the grange outside of function as it maybe from broadPeak and not narrowPeak
        gr0 = with(gr0 %>% mutate(strand=gsub('.', '*', strand), summit=start+summit),
                   GRanges(chr, IRanges(start,end), strand, qvalue=qvalue, fe=fe, summit=summit, protein=protein))
    }

    # means that user has defined a minimum percentage of overlap
    if(args$overlap != -1){
        merged = get_merged(gr0, min_ov = args$overlap)
        overlap = paste0(args$overlap, 'perc')
    } else {
        overlap = '1bp'
        merged = get_merged(gr0)
    }# }}}

# Computing and plotting venn diagram# {{{
#    granges = merged$annotated
#    venn_cnt = make_venn_cnt(granges)
    venn_cnt = make_venn_cnt(merged$annotated)
    v = venn_counts2venn(venn_cnt)

    target = paste0("Condition_", 
                    paste0((protein_data %>% dplyr::select(condition) %>%
                           .[[1]] %>% as.character() %>% unique()), collapse = "_AND_"),
                    "_Proteins_",
                    paste0((protein_data %>% dplyr::select(protein) %>%
                           .[[1]] %>% as.character() %>% unique()), collapse = "_AND_"), 
                    "_Ov", overlap)

    pdf(file.path(plot_path, paste0(target, ".pdf")))
    # Draw circles for 2 sets (ncol will be 3 because of the count column)
    # e.g venn_cnt
    # 2lox_2i-Esrrb 2lox_2i-Klf4 Counts
    # 0            0      0
    # 1            0  24778
    # 0            1   9120
    # 1            1   3621
    if(ncol(venn_cnt) == 3){
        plot(v[['vn']], doWeights=TRUE, type="circles", gp = v[['theme']])
    } else if(ncol(venn_cnt) == 5) {
        # ellipses only work for 4 sets and easier to read than chowRusjey so draw that too
        plot(v[['vn']], doWeights=FALSE, type="ellipses", gp = v[['theme']])
    } else {
         library(UpSetR)
#        plot(v[['vn']], doWeights=T, type="ChowRuskey",  gp = v[['theme']])
        upsetR_cnt = do.call(rbind, apply(venn_cnt, 1, function(x) {
                                              if(x['Counts'] != 0){
                                                  y = as.data.frame(t(x)) %>% dplyr::select(-Counts)
                                                  cnt = x['Counts']
                                                  y[1:cnt,] = y[1,]
                                                  return(y)
                                              }
                                    }))
        upset(upsetR_cnt, order.by="freq", mb.ratio=c(0.7, 0.3))
    }
    dev.off()# }}}

    tmp = merged$merged
    # if user provides a label use that to rename the metadata of the merged file
    if(args$anno != "NA"){
        print('Label provided - using that in BED file.')
        values(tmp) = values(tmp)[,grep('name', colnames(values(tmp)), invert=T)]
        tmp$name = args$anno
    }
    export.bed(tmp, con = file.path(args$output_path, paste0(target, '.bed')))
    export.bed(merged$annotated, con = file.path(args$output_path, paste0(target, '_complete_annotated.bed')))
}# }}}

main()
