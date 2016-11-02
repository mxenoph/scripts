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
parser$add_argument('-l', '--label', required = 'True', help= "This script always requires the labels for the different regions provided e.g NuRD:active-enhancers")
parser$add_argument('-k', '--keep', required = 'True', help= "Label given to the ranges we want to return classified")

args = parser$parse_args()

if(FALSE){
    args = list()
    args$bed = c("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/merged-ranges/Condition_7E12_2i_Proteins_M2_pooled_filtered_AND_Chd4_pooled_filtered_Ov1bp.bed", "/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-active_2015-08-10.bed", "/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-poised_2015-08-10-clean-for-use-on-the-fly.bed", "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.75.rtracklayer-2000-gene_start-500.bed")
    args$threads = 4
    args$output_path = getwd()
    args$label = '2i-NuRD:active-enhancers:poised-enhancers:promoters'
    args$keep  ='2i-NuRD'
}

plot_path = file.path(args$output_path, "plots")
dir.create(plot_path, recursive= TRUE)

ncores = args$threads
#}}}

# import bed files# {{{
proper_import = function(x, type="BED"){
    
    # local function reading any BED formatted file# {{{
    read_bed = function(x){
        # possible bed formats not recognised by rtracklayer without tweaking
        extraCols_gff2bed = c(source = "character", type = "character", phase = "numeric", attributes = "character")
        extraCols_narrowPeak = c(fe = "numeric", pvalue = "numeric", qvalue = "numeric", summit = "integer")
        extraCols_broadPeak = c(fe = "numeric", pvalue = "numeric", qvalue = "numeric")

        if(grepl('.narrowPeak', x)) {
            gr = tryCatch(import(x, format = "BED", extraCols = extraCols_narrowPeak),
                         error=function(e) stop(paste0('Cannot read in ', x, '. Please make sure it is in narrowPeak format\n', e)),
                         warning=function(w) w)
        } else if(grepl('.broadPeak', x)){
            gr = tryCatch(import(x, format = "BED", extraCols = extraCols_broadPeak),
                         error=function(e) stop(paste0('Cannot read in ', x, '. Please make sure it is in narrowPeak format\n', e)),
                         warning=function(w) w)
        } else {
            gr = tryCatch(import(x, format = "BED"),
                         error=function(e) e,
                         warning=function(w) w)
            if(inherits(gr, 'error')){
                gr = tryCatch(import(x, format = "BED", extraCols = extraCols_gff2bed),
                              error=function(e) e,
                              warning=function(w) w);
                if(inherits(gr, 'error')) stop(paste0('Cannot read in ', x, '. Please make sure it is in BED format'));
            }
        }
        return(gr)
    }# }}}
    
    # should be updated to include a read_gff function anf then use switch()
#    switch(as.character(x), 
#           BED = read_bed())
    return(read_bed(as.character(x)))
}# }}}

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

main = function(){# {{{
    source("/homes/mxenoph/source/Rscripts/granges-functions.R")
    library(GenomicRanges)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(gtools)
    library(rtracklayer)

    # Prepare the data# {{{
    if(length(args$bed) < 2) stop("Two or more sets of regions need to be provided.")

    sites = as.data.frame(args$bed)
    colnames(sites) = 'sites'

    # args$label is always required so no need to check it's there
    if(!grepl(':', args$label)) stop("Labels must be ordered as the input files and separated by :")
    tmp = data.frame(Class = unlist(strsplit(args$label, ":")))
    sites = cbind(sites, tmp)

    sites_grl = lapply(1:nrow(sites), function(x){
                           gr = proper_import(sites[x, 'sites'])
                           values(gr)$Class = sites[x, 'Class']
                           return(gr)
                    })
    names(sites_grl) = as.character(sites$Class)

    metadata_cols = lapply(sites_grl, function(x) colnames(values(x)))
    # keeping common values in metadata
    keep_metadata =  Reduce(intersect, metadata_cols)
    gr0 = unlist(GRangesList(lapply(sites_grl, function(x) x[,keep_metadata])))

    hits = findOverlaps(gr0)
    stopifnot(queryLength(hits) == subjectLength(hits))
    stopifnot(queryLength(hits) == length(gr0))
    clusters = extract_clusters_from_self_hits(hits)
    in_clusters = extractList(gr0, clusters[sapply(clusters, length) > 1])

    in_clusters_ranges = lapply(in_clusters, function(x) {
                                    tmp = paste(sort(unique(as.character(values(x)$Class))), collapse = "_AND_")
                                    x$feature = tmp
                                    x
                    })

    overlapping_keep = unlist(GRangesList(lapply(in_clusters_ranges, function(x) x[values(x)$Class == args$keep])))
    ov = findOverlaps(sites_grl[[args$keep]], overlapping_keep)
    # intergenic - not overlapping any of the features given
    non_overlapping = sites_grl[[args$keep]][! 1:length(sites_grl[[args$keep]]) %in% queryHits(ov)]
    non_overlapping$feature = non_overlapping$Class
    
    results = unlist(GRangesList(overlapping_keep, non_overlapping))
    # dropping unused levels
    results$feature = droplevels(results$feature)

    output_base = file.path(args$output_path,
                            file_path_sans_ext(basename(sites %>% filter(Class == args$keep) %>% .[['sites']] %>% as.character())))
    res = split(results, results$feature)
    lapply(names(res), function(x) {
               export.bed(res[[x]], paste0(output_base, '.', x, '.bed'))
                            })
    lapply(names(res), function(x) {
               export.bed(resize(res[[x]], width = 1, fix = 'center'),
                          paste0(output_base, '.', x, '-centers.bed'))
                            })
    
}#}}}


main()

#extra
#    in_clusters_ranges = GRangesList(lapply(in_clusters, function(x) x[values(x)$Class == '2i-NuRD']))# {{{
#    # contains 2 or more NuRD ranges
#    indirect = which(sapply(in_clusters_ranges, length) > 1)
#    indirect_no = sapply(in_clusters_ranges[indirect], length)
#
#    x = extractList(values(gr0)$Class, clusters[sapply(clusters, length) > 1])
#    in_clusters_names = sapply(x, function(i) {paste(unique(i), collapse = "_AND_")})
#    knames = unlist(sapply(1:length(indirect), function(i){
#               if( i == 1 ){
#                    stop_at = indirect[i] - 1
#                    f = in_clusters_names[1:stop_at]
#                    repeat_from = indirect[i]
#                    cat(paste0('1-', stop_at,'...', repeat_from, '(', indirect_no[i]), ')...')
#                    f = c(f, rep(in_clusters_names[repeat_from], indirect_no[i]))
#               } else {
#                   start_at = indirect[i - 1] + 1
#                   if(start_at != indirect[i]){
#                   stop_at = indirect[i] - 1
#                   cat(paste0(start_at, '-', stop_at))
#                   f = c(f, in_clusters_names[start_at:stop_at])
#                   }
#                   repeat_from = indirect[i]
#                    cat(paste0('...', repeat_from, '(', indirect_no[i], ')...'))
#                    f = c(f, rep(in_clusters_names[repeat_from], indirect_no[i]))
#                    if( i == length(indirect)) {
#                        from = indirect[i] + 1
#                        to = length(in_clusters_names)
#                        if(from <= to){
#                            f = c(f, in_clusters_names[from:to])
#                            cat(paste0(from, '-', to, '\n'))
#                        }
#                    }
#               }
#               return(f)
#               }))
#    in_clusters_ranges_0 = unlist(in_clusters_ranges)
#    values(in_clusters_ranges)$Class = knames# }}}
