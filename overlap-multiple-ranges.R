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
parser$add_argument('-o', '--output_path', help= "Output path")

args = parser$parse_args()

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
    groups = characterlist(mclapply(mcols(merged_ranges)$revmap,
                                    function(x){
                                        x = mixedsort(unique(values(input_granges[x])[['protein']]))
                                        if (length(x) != 1) {
                                            x = paste(x, collapse = "_and_")
                                        }
                                        return(x)
                                    }, mc.cores = 4))
    mcols(merged_ranges)$group = groups
    # if group name contains _and_ then range is merged from connected tf ranges
    multiple = grepl("*_and_*", groups)
    multiple_peaks = merged_ranges[multiple]
    single_peaks = merged_ranges[!multiple]
    return(GRangesList('annotated' = merged_ranges, 'merged' = multiple_peaks, 'single' = single_peaks))
} # }}}

main = function(){
    print('in main')
    source("/homes/mxenoph/source/Rscripts/granges-functions.R")
    source("/homes/mxenoph/source/Rscripts/plotting-functions.R")
    library(GenomicRanges)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(gtools)


    # Make venn count table# {{{
    make_venn_cnt = function(granges){
        stopifnot(is(granges, "GenomicRanges"))
        stopifnot(any(grepl('Group', colnames(mcols(granges)))))

        sets = levels(as.factor(unlist(mcols(granges)$Group)))
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
                counts_table = table(grepl(regex, mcols(granges)$Group))
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
                counts_table = table(grepl(regex, mcols(granges)$Group))
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

    #chip_path = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/before_22.9/macs2/macs2/sharp"
    #TF_data = list.files(pattern = "*[2i]*.narrowPeak", path = chip_path, full.name = T)

    # Prepare the data# {{{
    protein_data = as.data.frame(args$bed)
    colnames(protein_data) = 'protein_data'

    tmp = protein_data %>% .[[1]] %>% as.character() %>%
        basename() %>% file_path_sans_ext %>% str_replace_all("_peaks", "") %>% strsplit("-")
    tmp = do.call(rbind.data.frame, tmp)
    colnames(tmp) = c('condition', 'protein')
    protein_data = cbind(protein_data, tmp)

    proteins_df = lapply(protein_data[['protein_data']], function(x) {
                    tmp = protein_data %>% filter(protein_data == x) %>% unite(id, condition, protein, sep="-") %>% select(id) %>% .[[1]]
                    x = import_narrowPeak_df(as.character(x))
                    tmp = rep(tmp, nrow(x))
                    x %>% mutate(protein = tmp) })
# }}}

    # keep all peaks in one GRange object and find all overlaps# {{{
    gr0 = do.call(rbind, proteins_df)
    # make the grange outside of function as it maybe from broadPeak and not narrowPeak
    gr0 = with(gr0 %>% mutate(strand=gsub('.', '*', strand), summit=start+summit),
               GRanges(chr, IRanges(start,end), strand, qvalue=qvalue, fe=fe, summit=summit, protein=protein))

    # means that user has defined a minimum percentage of overlap
    if(args$overlap != -1){
        merged = get_merged(gr0, min_ov = args$overlap)
    }# }}}
    
    esc = a$merged
    mcols(esc) = rep('NuRD_ESC', length(esc))
    colnames(mcols(esc)) = 'TF'
    hits = findOverlaps(gr0)

    # Merge the ranges in 'gr0' that are connected via one or more hits in 'hits'.
    # {{{
    gr1 = mergeConnectedRanges(gr0, hits)

    groups = CharacterList(mclapply(mcols(gr1)$revmap,
                                    function(x){
                                        x = mixedsort(unique(values(gr0[x])[['protein']]))
                                        if (length(x) != 1) {
                                            x = paste(x, collapse = "_AND_")
                                        }
                                        return(x)
                                    }, mc.cores = ncores))# }}}

    mcols(gr1)$Group = groups
    # if group name contains _AND_ then range is merged from connected protein ranges
    multiple = grepl("*_AND_*", groups)
    multiple_peaks = gr1[multiple]
    single_peaks = gr1[!multiple]

    # Computing and plotting venn diagram# {{{
    venn_cnt = make_venn_cnt(gr1)
    v = venn_counts2venn(venn_cnt)

    target = paste0("Condition_", 
                    paste0((protein_data %>% dplyr::select(condition) %>%
                           .[[1]] %>% as.character() %>% unique()), collapse = "_AND_"),
                    "_Proteins_",
                    paste0((protein_data %>% dplyr::select(protein) %>%
                           .[[1]] %>% as.character() %>% unique()), collapse = "_AND_"))

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
    } else {
        plot(v[['vn']], doWeights=FALSE, type="ellipses", gp = v[['theme']])
        plot(v[['vn']], doWeights=T, type="ChowRuskey",  gp = v[['theme']])
    }
    dev.off()# }}}


    chip_path = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/new/macs2/sharp"
    epi_data = list.files(pattern = "7E12_EpiSC.*filtered_peaks.narrowPeak", path = chip_path, full.name = T, ignore.case = T)
    epi_data = data.frame(file = epi_data)
    idr_path = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/new/IDR"
    esc_data = list.files(pattern = "*2i.*filtered.*narrowPeak", path = idr_path, full.name = T)
    esc_data = data.frame(file = esc_data)
    df = rbind(esc_data, epi_data)

    tmp = df %>% .[[1]] %>% as.character() %>%
        basename() %>% file_path_sans_ext %>% str_replace_all("_filtered_peaks", "") %>%
        str_replace_all("_conservative.*", "") %>% strsplit("-")
    tmp = do.call(rbind.data.frame, tmp)
    colnames(tmp) = c('condition', 'TF')
    df = cbind(df, tmp)

    esc_df = lapply(df[['file']], function(x) {
                    tmp = df %>% filter(file == x) %>% unite(id, condition, TF, sep="-") %>% dplyr::select(id) %>% .[[1]]
                    x = import_narrowPeak_df(as.character(x))
                    tmp = rep(tmp, nrow(x))
                    x %>% mutate(TF = tmp) })
    names(esc_df) = df %>% unite(id, condition, TF, sep = '-') %>% dplyr::select(id) %>% .[[1]]



    gr0 = do.call(rbind, esc_df[grepl('2i', names(esc_df))])
    # make the grange outside of function as it maybe from broadPeak and not narrowPeak
    gr0 = with(gr0 %>% mutate(strand=gsub('.', '*', strand), summit=start+summit),
               GRanges(chr, IRanges(start,end), strand, qvalue=qvalue, fe=fe, summit=summit, TF=TF))
    # NuRD in ESCs
    a = get_merged(gr0, min_ov = 0)

    gr0 = do.call(rbind, esc_df[!grepl('2i', names(esc_df))])
    # make the grange outside of function as it maybe from broadPeak and not narrowPeak
    gr0 = with(gr0 %>% mutate(strand=gsub('.', '*', strand), summit=start+summit),
               GRanges(chr, IRanges(start,end), strand, qvalue=qvalue, fe=fe, summit=summit, TF=TF))
    b = get_merged(gr0, min_ov = 0)

    esc = a$merged
    mcols(esc) = rep('NuRD_ESC', length(esc))
    colnames(mcols(esc)) = 'TF'

    epi = b$merged
    mcols(epi) = rep('NuRD_EpiSC', length(epi))
    colnames(mcols(epi)) = 'TF'

    gr0 = c(esc, epi)
    esc_epi = get_merged(gr0, min_ov = 0)
    epi_only = esc_epi$annotated[grepl("^NuRD_EpiSC$", mcols(esc_epi$annotated)[['Group']])]

    venn_cnt = make_venn_cnt(esc_epi$annotated)
    v = venn_counts2venn(venn_cnt)
    pdf('test-ESC-Epi.pdf')
    plot(v[['vn']], doWeights=TRUE, type="circles", gp = v[['theme']])
    dev.off()

    source("~/source/Rscripts/annotation-functions.R")
    get_annotation('mm10')
    tmp = epi_only
    mcols(tmp) = paste(seqnames(tmp), start(tmp), end(tmp), sep = '_')
    tss_with_name = as.data.frame(mcols(tss_window)['gene_id'])
    tss_with_name = tss_with_name %>% left_join(as.data.frame(mcols(genes)[c('gene_id','gene_name')]), by = 'gene_id')
    mcols(tss_window) = tss_with_name
    epi_only_annotated = annotatePeakInBatch(epi_only, AnnotationData = tss_window)

    keep =  as.data.frame(epi_only_annotated) %>%
            filter(!insideFeature %in% c("upstream", "downstream")) %>%
            mutate(gene_id = feature) %>%
            left_join(as.data.frame(mcols(genes)[c('gene_id', 'gene_name')]), by = "gene_id")
    write.table(keep, file = "NuRD_EpiSC_only_tss_overlap", col.names= TRUE, row.names = FALSE, sep = "\t", quote = F)
    library(rtracklayer)
}

main()
