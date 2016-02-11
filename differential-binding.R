#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library = function (...) suppressMessages(base::library(...))
select = dplyr::select
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")
source("~/source/Rscripts/granges-functions.R")

parser =  ArgumentParser(description="Perform differential binding analysis")
parser$add_argument('-s', '--sheet', metavar= "file", required='True', type= "character", help= "Sample sheet must have SampleID,Condition,Treatment,Replicate,bamReads,bamControl,Peaks (last 3 are the path to the respective bed files)")
parser$add_argument('-a', '--assembly', type= "character", default='mm10', help= "Give preferred assembly e.g. mm10. Default: mm10")
parser$add_argument('-d', '--de', metavar = "path", type= "character", default='', help= "Path to look for the list of de expressed genes")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
parser$add_argument('-g', '--gtf', metavar= "file", type= "character",
                    default = "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.75.gtf",
                    help= "GTF file. Looks for the rtracklayer generated gtfs for that annotation. If not found it runs get-annotation-rscript.R")
parser$add_argument('-e', '--active_enhancers', metavar= "file", type= "character",
                    default = "/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-active_2015-08-10.bed",
                    help= "BED file of active enhancers")
parser$add_argument('-p', '--poised_enhancers', metavar= "file", type= "character",
                    default = "/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-poised_2015-08-10.bed",
                    help= "BED file of active enhancers")

args = parser$parse_args()

label = file_path_sans_ext(basename(args$sheet))
output_path = file.path(args$out, 'DiffBind', label)
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)
#}}}

# Load Packages
library(DiffBind)
library(dplyr)
library(tidyr)
# For str_replace_all function
library(stringr)
# For mixedsort
library(gtools)
library(ggplot2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
# Color scheme for divergent colors.
divergent_colors = colorRampPalette(c('#603D71', 'white', '#A4B962'))(30)
heat_cols = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(100)

# Less intrusive defaults for heatmap.2.
heatmap.2 = function (...) 
            gplots::heatmap.2(..., trace = 'none', density.info = 'none',
                              col = heat_cols)
pheatmap = function (...)
            pheatmap::pheatmap(..., trace = 'none', density.info = 'none',
                               #color = divergent_colors, border_color = NA,
                               color = heat_cols, border_color = NA,
                               fontsize = 10, fontsize_row = 6,
                               show_colnames = TRUE)
add_rownames = function(df, var = 'rowname') {
    stopifnot(is.data.frame(df))
    rowname_df = setNames(data_frame(rownames(df)), var)
    cbind(rowname_df, df)
}
## Functions

# Read in and format data# {{{
import_data = function(cnf_df){

    sapply(cnf_df[['Peaks']],
           function(x){
               cmd = sprintf("awk -F \"\t\" \'NR==1 {if($0 ~ /^track/){print 1}else{print 0}}\' %s", x)
               # If track=1 file contains ucsc track line
               track = pipe(cmd, open = "r")
               # close connection
               on.exit(close(track))

               flag = scan(track, quiet=TRUE)
               if(flag==1){
                   cmd = sprintf("sed -i \'1d\' %s", x)
                   scan(pipe(cmd), quiet=TRUE)
               }
           })

    # by default dba plots the correltaion heatmap. Set it to False and print on call
    peakset = dba(sampleSheet = cnf_df, bCorPlot = FALSE)
#    print(peakset)
    return(list(cnf_df = cnf_df, peakset = peakset))
}
# }}}

# Differential Binding Analysis: doing all the work # {{{
differential_binding = function (cnf, protein){
    tmp = import_data(cnf)
    config = tmp[['cnf_df']]
    peakset = tmp[['peakset']]

    # Using only this peak caller data, a correlation heatmap can be generated 
    # which gives an initial clustering of the samples
    pdf(file.path(plot_path, paste0(protein, '.pdf')))
    # if only using one peak file then dba.plotHeatmap will throw error for soem reason
    # as correlations will all be 1
    if (length(levels(factor(cnf[['Peaks']]))) > 1 ) {
        # For this plot clustering is based on peak scores. Peaks identified somewhere
        # in the experiment, but not called for a specific sample, will have a missing 
        # score (set to -1), so the score vectors being correlating may be quite sparse
        dba.plotHeatmap(peakset, main = 'Correlation on peak scores (from peak caller)')
        # Adding consensus peaks for all conditions-proteins that are the same and only
        # differ in terms of replicates
        peakset = dba.peakset(peakset, consensus = -DBA_REPLICATE, minOverlap = 0.5)
        dba.count = function (...)
            DiffBind::dba.count(..., minOverlap = 0.5)
    } else {
        x = unique(cnf[['Peaks']])
        if (file_ext(x) == 'narrowPeak') {
            peaks_gr = import_narrowPeak(x)
        } else if (file_ext(x) == 'broadPeak') {
            peaks_gr = import_broadPeak(x)
        }
        peakset = dba.peakset(peakset, consensus = -DBA_REPLICATE, minOverlap = 1)
        dba.count = function (...)
            DiffBind::dba.count(..., minOverlap = 0)
    }

    ov_rate = dba.overlap(peakset, mode=DBA_OLAP_RATE)
    # If this curve drops off too quickly (worse than geometric), that indicates that there
    # is little agreement between the peaksets called for each sample
    p = ggplot(data.frame(common = ov_rate, sets = 1:length(ov_rate)), aes(y=common, x=sets)) + geom_point(shape=19) + geom_line()
    p + xlab('# peaksets') + xlab("# common peaks") + ggtitle(protein) + theme_classic()

    # https://support.bioconductor.org/p/63466/
    # DBA_SCORE_TMM_READS_EFFECTIVE normalizes to the number of reads actually overlapping peaks 
    # (should only be used if you expect most of the peaks to have similar binding levels)
    # DBA_SCORE_TMM_READS_FULL, which will normalize to the overall depth of sequencing in each library
    counts = dba.count(peakset, score = DBA_SCORE_TMM_MINUS_FULL,
                       # only include peaks in at least this many peaksets
#                       minOverlap = 2,
#                       peaks = peakset$masks$Consensus,
                       # Assuming bam file is ddup
                       # bRemoveDuplicates = TRUE,
                       bScaleControl = TRUE, bParallel = TRUE,
                       bCorPlot = FALSE)
    dba.plotHeatmap(counts, main = 'Correlation on peak read count score')
    dba.plotHeatmap(counts, main = 'Correlation on peak read count score', correlations=F)
    dev.off()

    contrast = dba.contrast(counts, categories = DBA_CONDITION, minMembers = 2)
    
    diff_bound = dba.analyze(contrast, method = DBA_DESEQ2,
                             # indicating control reads are subtracted
                             bSubControl= TRUE,
                             # full library size used for normalization, effective library size
                             # normalization is preferred if binding levels are expected to be similar
                             # between samples
                             bFullLibrarySize = TRUE,
                             bTagwise = TRUE,
                             bParallel = TRUE,
                             bCorPlot = FALSE)

    for (index in 1:length(contrast$contrasts)){
        comparison = paste(protein,
                           paste(contrast$contrasts[[index]]$name1,
                                 contrast$contrasts[[index]]$name2, sep='-'),
                           sep = '_')
        
        # write report with all sites irrespective if differentially bound or not
        results =  dba.report(diff_bound, method = DBA_DESEQ2,
                              contrast = index,
                              # add count data for individual samples
                              bCounts = TRUE,
                              # only include sites with an absolute Fold value greater than equal
                              # fold= 2
                              # Threshold will be FDR and not pval
                              bUsePval = FALSE,
                              th = 1)
        write.table(as.data.frame(results),
                    file = file.path(output_path,
                                     paste0(comparison, '.tsv')), sep="\t", quote=FALSE, row.names=FALSE)
        
        results =  dba.report(diff_bound, method = DBA_DESEQ2,
                              contrast = index,
                              # add count data for individual samples
                              bCounts = TRUE,
                              # only include sites with an absolute Fold value greater than equal
                              # fold= 2
                              # Threshold will be FDR and not pval
                              bUsePval = FALSE,
                              th = 0.05)
        
        # So that the function doesn't crush for cases where no db detected 
        if( length(results) != 0 ) {
        
            write.table(as.data.frame(results),
                        file = file.path(output_path,
                                         paste0(comparison, '-db.tsv')), sep="\t", quote=FALSE, row.names=FALSE)
            
            pdf(file.path(plot_path, paste0(comparison, '.pdf')))

            dba.plotMA(diff_bound, contrast = index, method = DBA_DESEQ2,
                       th = 0.05, bUsePval = FALSE )

            dba.plotPCA(diff_bound, contrast = index, method = DBA_DESEQ2,
                        th = 0.05, bUsePval = FALSE,
                        label = DBA_ID)

            groups_pvals = dba.plotBox(diff_bound, contrast = index, method = DBA_DESEQ2, bAll = TRUE, pvalMethod = NULL)

            dba.plotHeatmap(diff_bound, method = DBA_DESEQ2,
                            contrast = index,
                            # Threshold will be FDR and not pval
                            bUsePval = FALSE,
                            score = DBA_SCORE_TMM_MINUS_FULL,
                            th = 0.05)

            dba.plotHeatmap(diff_bound, method = DBA_DESEQ2,
                            contrast = index,
                            # Threshold will be FDR and not pval
                            bUsePval = FALSE,
                            score = DBA_SCORE_TMM_MINUS_FULL,
                            correlations=FALSE,
                            th = 0.05)

            dba.plotHeatmap(diff_bound, method = DBA_DESEQ2,
                            contrast = index, mask = diff_bound$masks$All,
                            # Threshold will be FDR and not pval
                            bUsePval = FALSE,
                            score = DBA_SCORE_TMM_MINUS_FULL,
                            correlations=FALSE,
                            th = 0.05)
            dev.off()
        }
    }
}
# }}}

aggregate_counts = function(fs){ # {{{
    first = TRUE
    for (x in 1:length(fs)) {
        tmp = read.delim(fs[x], head = T, sep = "\t", stringsAsFactors = F)

        rename_map = c(paste("FC", names(fs)[x], sep = "."),
                       paste("FDR", names(fs)[x], sep = "."))
        combine_by = c('seqnames', 'strand', 'start', 'end')
        tmp = tmp %>% dplyr::select(seqnames:end, strand, Fold, FDR:ncol(.)) %>%
              rename_(.dots = setNames(list("Fold", "FDR"), rename_map))

        if(first) {# {{{
            ofs = tmp
            first = FALSE
        } else {
            ofs = inner_join(ofs, tmp, by = combine_by)
            # since samples at h0 will be in all comparisons expression columns for those will be
            # duplicated and dplyr add '.x' or '.y' to them depending from which df they came from
            # remove extension, find duplicates and remove. This would work nicer with rename() but I
            # can't get it to work on selected columns
            colnames(ofs) = gsub("\\.x|\\.y$", "", colnames(ofs))
            discard = ofs %>% names(.) %>% duplicated()
            ofs = ofs[,!discard]
        }
    }# }}}
    return(ofs)
}
# }}}

# Prepare a data.frame containing counts/TMM for multiple conditions/samples for plotting# {{{
# ca is the name for column annotation
format_cdf = function(cdf, statistic = "FDR", ca = "Time", threshold = 0.05){
    source("~/source/Rscripts/ggplot-functions.R")
    # Significant changes between any samples in the cdf
    significant = dplyr::select(cdf, contains(statistic)) %>% apply(1, function(x) min(x[!is.na(x)]) < threshold)

    if(any(significant)) {
        plot_data = cdf[significant,]
        # Save because it is lost with dplyr
        row_names = plot_data[['rowname']]
        
        # Row annotation needs to be defined before selecting columns
        row_annotation = dplyr::select(plot_data, contains(statistic)) %>%
                         mutate_each(funs(ifelse(. > 0.05 | is.na(.), 'non significant', 'significant')))
        row_annotation = row_annotation %>% dplyr::select(one_of(rev(mixedsort(colnames(row_annotation)))))
        rownames(row_annotation) = row_names

        plot_data = plot_data %>% dplyr::select(-starts_with('FC'),
                                                -starts_with(statistic),
                                                -seqnames, -strand, -start, -end,
                                                -rowname)

        # Formatting time in row_annotation and plot_data and consequenctly column_annotation and ac# {{{
        if(ca == "Time") {
            convert_time = function(x) {
                index = grepl('m', x)
                values = x[index]
                # R introduces X in front of the string if starting with number
                values = gsub("^X", '', values)
                values = as.numeric(gsub("m", '', values))
                values = gsub("^", "h", as.character(values / 60))
                return(list(index = index, values = values))
            }
            timepoints = data.frame(x =  gsub("^.*\\.", "", colnames(row_annotation))) %>%
                            tidyr::separate(x, into = c("timepoint1", "timepoint2"), sep = '-')
            x = convert_time(timepoints[['timepoint1']])
            timepoints[['timepoint1']][x$index] = x$values

            x = convert_time(timepoints[['timepoint2']])
            timepoints[['timepoint2']][x$index] = x$values
            colnames(row_annotation) = (timepoints %>% tidyr::unite(x, timepoint1, timepoint2, sep = "-") %>% .[[1]])

            timepoints = data.frame(x = colnames(plot_data)) %>% tidyr::separate(x, into = c("timepoint", "factor", "rep"), sep="_")
            x = convert_time(timepoints[['timepoint']])
            timepoints[['timepoint']][x$index] = x$values
            colnames(plot_data) = (timepoints %>% tidyr::unite(x, timepoint, factor, rep, sep = "_") %>% .[[1]])
        }# }}}
        
        # Reordering row_annotation columns so that all comparison involving zero show up first on heatmap
        row_annotation = row_annotation %>% select(one_of(mixedsort(names(.))))
        plot_data = plot_data %>% dplyr::select(match(mixedsort(colnames(plot_data)), colnames(plot_data)))

        plot_data = as.matrix(plot_data)
        rownames(plot_data) = row_names

        # Column names should be sample_rep with no '_' in "sample"
        column_annotation = data.frame( x = gsub("_.*", '', colnames(plot_data)))
        colnames(column_annotation) = ca
        rownames(column_annotation) = colnames(plot_data)

        color = c('slategrey', 'violetred')
        anno_colors = lapply(1:ncol(row_annotation), function(x) c("non significant" = color[1],
                                                                   "significant" = color[2]))
        names(anno_colors) = colnames(row_annotation)

        cool_colors = colorRampPalette(brewer.pal(9,"YlGnBu"))(8)
        tmp = cool_colors[1:length(unique(column_annotation[[ca]]))]
        names(tmp) = mixedsort(as.character(unique(column_annotation[[ca]])))
        anno_colors[[ca]] = tmp

        return(list(plot_data = plot_data,
                    ra = row_annotation, ca = column_annotation,
                    ac = anno_colors))
    } else {
        stop('No differentially bound regions.')
    }
} # }}}

annotate = function(grange, annotations, output_path, target,# {{{
                    priority = c('2000-gene_start-500', 'first-exon', 'active_enhancers',
                                 'poised_enhancers', 'exons', 'first-intron', 'introns', 'genes')) {

    stopifnot(is(grange, "GRanges"))
    total = length(grange)

    per_region = data.frame(closest_feature = rep(NA, total))
    rownames(per_region) = names(grange)

    per_feature = setNames(vector("list", length(priority)), priority)
    bound = data.frame(Gene = character(0), Bound = character(0))
    # '5000-tss-2000' contains all possible TSS for a gene
    for (i in priority) {
        ov = findOverlaps(annotations[[i]], grange)
        # using names and not index as that will change when subsetting granges downstream
        hits = names(grange)[(1:length(grange)) %in% unique(subjectHits(ov))]
        per_region[hits, 'closest_feature'] = rep(i, length(hits))

        print(paste0('Calculating %peaks overlapping ', i, '.', sep=''))

        if(i %in% grep('enhancer', priority, invert=T, value=T)){
            if(grepl('exon', i) || grepl('intron', i)) {
                    x = 'ensembl_id'
            } else {
                    x = 'ensembl_gene_id'
            }
            # Even if queryHits are unique the Gene names might not be as grange has all possible TSS 
            # for a gene
            tmp = data.frame(Gene = unique(values(annotations[[i]][unique(queryHits(ov))])[[x]])) %>%
                mutate(Bound = i)
            # Find all genes that are not bound at a feature with higher priority
            # Avoiding duplicated values in bound$Gene
            tmp = anti_join(tmp, bound, by="Gene")
            bound = dplyr::bind_rows(bound, tmp)

            # For each region keep the ID of feature is binding to. For regions binding to more than one gene etc
            # will only keep the first occurrence in the ov table
            tmp = as.data.frame(ov) %>% mutate(qNames = values(annotations[[i]])[[x]][queryHits],
                                               sNames = names(grange)[subjectHits])
            per_region[tmp[['sNames']], 'feature_id'] = tmp[['qNames']]
        } else {
            x = 'name'
            if (! "name" %in% colnames(values(annotations[[i]])) | all(is.na(values(annotations[[i]])[['name']]))) {
                values(annotations[[i]])[['name']] = paste(seqnames(annotations[[i]]), 
                                                           start(annotations[[i]]), end(annotations[[i]]), sep=":")
            }
            tmp = as.data.frame(ov) %>% mutate(qNames = values(annotations[[i]])[[x]][queryHits],
                                               sNames = names(grange)[subjectHits])
            per_region[tmp[['sNames']], 'feature_id'] = tmp[['qNames']]
        }

        n_overlapping = length(unique(subjectHits(ov)))
        grange = grange[! names(grange) %in% hits ]
        per_feature[[i]] = n_overlapping
    }

    per_region[['closest_feature']] = gsub("^NA;", '',  per_region[['closest_feature']])

    bound_table = data.frame(Gene = values(annotations[['genes']])[['ensembl_gene_id']]) %>% left_join(bound, by="Gene")
    write.table(bound_table,
                file = file.path(output_path, paste0(target,'.binding-per-gene.tsv')),
                quote = F, row.names = F, sep="\t")
    write.table(bound,
                file = file.path(output_path, paste0(target,'.bound-genes-only.tsv')),
                quote = F, row.names = F, sep="\t")
    write.table(add_rownames(per_region, var = "region"),
                file = file.path(output_path, paste0(target,'.binding-per-region.tsv')),
                quote = F, row.names = F, sep="\t")

    if(per_feature[['genes']] != 0){
        stop('I find regions overlapping with genes even if tss, introns and exons have higher priority. Something went wrong')
    } else {
        per_feature = per_feature[which(names(per_feature)!='genes')]
    }

    per_feature[['intergenic']] = total - sum(unlist(per_feature))
    per_feature = unlist(per_feature)
    per_feature = per_feature / total
    return(list(per_feature = per_feature, per_region = per_region))
}
# }}}

# When wanting to retrieve a cluster call the function to get IDs in particular cluster# {{{
# df = tmp$plot_data, regions = bdf
get_cluster = function(results, k = 1, output_path, label, id = "ID"){
    # k must be equal to cutree_rows used in pheatmap for defining number of clusters
    # http://stackoverflow.com/questions/27820158/pheatmap-in-r-how-to-get-clusters
    clusters = data.frame(cluster_no = cutree(results$tree_row, k = k))
    heatmap_order = add_rownames(clusters[results$tree_row$order,,drop = F], var = id)
    write.table(heatmap_order, file.path(output_path, paste0(label, '-', k, 'clusters-heatmap_order.tsv')),
                sep = "\t", col.names = T, row.names = F, quote = F)
    return(heatmap_order)
}
# }}}



# annotate genes and plot TPM and binding affinity for all de genes# {{{
per_de_gene = function(comparisons){
    
    for(p in proteins){
    diffbind = read_diffbind(comparisons, p)
    db_regions = with(diffbind,
                      GRanges(seqnames = seqnames,
                              IRanges(start, end),
                              strand = strand))
    names(db_regions) = paste(seqnames(db_regions), start(db_regions), end(db_regions), sep = "_")
    values(db_regions) = diffbind %>% select(strand:ncol(.), -strand)
    annotated = annotatePeakInBatch(genes, AnnotationData = db_regions)
    annotated = as.data.frame(annotated)
    # since annotating genes with peaks (and not the other way round) 
    # gene identifiers will be under column peaks, select genes d.e.,
    # reorder based on clusters from gene expression heatmap and subset db_regions
    keep =  annotated %>% filter(peak %in% heatmap_order) %>%
            mutate(peak = factor(peak, levels = heatmap_order)) %>%
            arrange(peak) %>% select(feature) %>% .[[1]]
    plot_data_db = db_regions[keep]
    names(plot_data_db) = NULL
    keep = make.unique(keep, sep=".")
    plot_data_db = as.data.frame(plot_data_db)

    row_annotation = select(plot_data_db, contains("padj")) %>%
                     mutate_each(funs(ifelse(. > 0.05 | is.na(.), 'non significant', 'significant')))
    row_annotation = row_annotation %>% select(one_of(mixedsort(colnames(row_annotation))))
    rownames(row_annotation) = keep
    plot_data_db = plot_data_db %>% select(matches(p), -contains("log2"), -contains("padj")) %>%
                   log2() %>% as.matrix()
    rownames(plot_data_db) = keep

    column_annotation = data.frame( Time = gsub("_.*", "", colnames(plot_data_db)))
    rownames(column_annotation) = colnames(plot_data_db)

    anno_colors = lapply(1:ncol(row_annotation), function(x) c("non significant" = "black",
                                                               "significant" = "firebrick"))
    names(anno_colors) = colnames(row_annotation)
    tmp = cool_cols[1:length(unique(column_annotation[['Time']]))]
    names(tmp) = mixedsort(unique(column_annotation[['Time']]))
    anno_colors$Time = tmp

    pdf(file.path(plot_path, paste0('expression-', p, '.pdf')), paper = 'a4')
    pheatmap(plot_data_db,
             clustering_distance_rows = drows,
                       clustering_method = "ward.D2", scale = 'row', cutree_rows = 10,
                       cluster_cols = TRUE, cluster_rows = TRUE,
                       annotation_row = row_annotation, annotation_col = column_annotation,
                       annotation_legend = T, annotation_colors = anno_colors,
                       show_rownames = FALSE)
    qpcr_keep = annotated %>% filter(peak %in% qpcr[['gene_id']]) %>%
                mutate(peak = factor(peak, levels = qpcr_heatmap_order)) %>%
                arrange(peak) %>% select(feature) %>% .[[1]]
    pheatmap(plot_data_db[qpcr_keep,],
             clustering_distance_rows = qpcr_drows,
                       clustering_method = "ward.D2", scale = 'row', cutree_rows = 2,
                       cluster_cols = TRUE, cluster_rows = TRUE,
                       annotation_row = row_annotation[qpcr_keep,], annotation_col = column_annotation,
                       annotation_legend = T, annotation_colors = anno_colors,
                       show_rownames = TRUE)
    dev.off()
    }
} # }}}


# combine db regions with expression data for nearest gene, annotation = ensembl, jan2013 etc if not in_file # {{{
combine_expression = function(comparisons, annotation = "in_file"){
# anotate db sites and plot TPM for associated genes# {{{
per_db_region = function(db_regions, deseq_results, pattern = 'Expression'){
    genes = retrieve_annotation(deseq_results)
    # IMPORTANT: use read.delim instead of read.table as quote="\"" is set by default
    # if that is not set in read.table then the file is not read through he end
    db_regions = read.delim(db_regions, head = T, sep = "\t")
    deseq_results = read.delim(deseq_results, head = T, sep = "\t")
    
    # Keep TMM for replicates as I will need that downstream
    # It's not easy to keep it in granges object without knowing the names
    db_coverage = db_regions %>%
                    mutate(peak = paste(seqnames, start, end, sep = "_")) %>%
                    dplyr::select(-seqnames, -start, -end, -width, -strand, -Conc, -Fold, -p.value, -FDR)
    db_regions = with(db_regions,
                      GRanges(seqnames = seqnames,
                              IRanges(start, end),
                              strand = strand,
                              Fold = Fold,
                              FDR = FDR))
    names(db_regions) = paste(seqnames(db_regions), start(db_regions), end(db_regions), sep = "_")


    # data annotation and plotting# {{{
    annotated = annotatePeakInBatch(db_regions, AnnotationData = genes)
    # Filtering, keeping only db regions in or including genes or withn 5Kb from a gene TSS or TTS
    annotated = as.data.frame(annotated) %>%
                filter(insideFeature %in% c("inside", "includeFeature") |
                       (insideFeature %in% c("upstream", "downstream") & shortestDistance < 5000)) %>%
                # joining based on gene_id downstream
                rename(gene_id = feature) %>%
                inner_join(deseq_results, by = "gene_id") %>%
                # select fold, fdr for peak and log2FoldChange, padj and tpms for gene
                dplyr::select(Fold:gene_id, log2FoldChange, padj:ncol(.)) %>%
                # filtering based on significance of expression change
                filter(!is.na(padj)) %>% mutate(DE = ifelse(padj < 0.05, 'p < 0.05', 'p >= 0.05')) %>%
                inner_join(db_coverage, by = 'peak')

    plot_data = annotated %>% dplyr::select(contains(pattern),
                                     one_of(colnames(db_coverage)[3:ncol(db_coverage)]), gene_id, DE)
    row_annotation = data.frame(Expression = plot_data[['DE']],
                              row.names = paste(plot_data[['peak']], plot_data[['gene_id']], sep = "-"))

    # selecting only TPM columns and scaling based on the rowMeans for the respective experiment
    # RNA seq or ChiP-seq
    plot_data  = scale_rows(plot_data %>% dplyr::select(-peak, -gene_id, -DE))

    # Id rownames are not present (and dplyr removes them, so...) the annotation row
    # in the heatmap will not work
    rownames(plot_data) = rownames(row_annotation)

    pheatmap(plot_data,
             clustering_distance_rows = "correlation", scale = 'none',
             cluster_cols = FALSE, cluster_rows = TRUE,
             annotation_row = row_annotation, annotation_legend = T,
             show_rownames = FALSE)

    p = ggplot(annotated, aes(x = log2FoldChange, y = Fold, colour = DE)) + geom_point()
    p = p + theme_classic() + xlab("log2 Fold Change of Nearest Gene") + ylab("log2 Fold Change of DB regions") + ggtitle(label)
    p = p + scale_colour_manual(values = c("#AA3929", "#8E9CA3"))
    print(p)

    # Performing a GSEA on the genes associated with a db region
    gene_scores = deseq_results %>% filter(!is.na(padj)) %>% select(gene_id, padj)
    gene_scores = with(gene_scores, structure(padj, names = as.character(gene_id)))
    # }}}

    # Doing GSEA on de genes and de and db genes
#        ob = get_set_enrichment(gene_scores = gene_scores, label = label)
#        ob = get_set_enrichment(gene_scores = gene_scores, bound = (select(annotated_df, gene_id) %>% .[[1]]), label = label)
}
# }}}
    library(ChIPpeakAnno)

    scale_rows = function(tpms, pattern = regex("^Expression"), reverse_pattern = regex("^[^Expression]")){# {{{
        # IMPORTANT: if instead mutate_each_(funx(...), ~matches(...)) I use mutate_each(funs(...), matches(...))
        # then the selected columns are not mutated but new columns holding the correct values are created.
        # Those are names varX where X is the index of the column after subseting the table
        # Once done scaling remove the rowMeans columns
        tpms = tpms %>% log() %>%
               mutate(Expression_RowMean = rowMeans(.[grep(pattern, names(.))]),
                      Binding_RowMean = rowMeans(.[grep(pattern, names(.), invert = T)])) %>%
                      mutate_each_(funs(. - Expression_RowMean), vars = ~matches(pattern)) %>%
                      mutate_each_(funs(. - Binding_RowMean), vars = ~matches(reverse_pattern)) %>%
                      select(-matches('RowMean'))
       return(tpms)
    }# }}}


    # Match chip-seq with rna-seq# {{{
    for( x in 1:nrow(comparisons) ) {
        db_regions = as.character(comparisons[x, 'comparisons'])
        deseq_results = as.character(comparisons[x, 'expression_data'])
        
        filename = paste0(file_path_sans_ext(basename(db_regions)),
                          "_RNAseq-",
                          file_path_sans_ext(basename(deseq_results)),
                          ".pdf")
        label = paste0("binding affinity: ",
                       file_path_sans_ext(basename(db_regions)),
                       ", expression: ",
                       file_path_sans_ext(basename(deseq_results)))
    
        pdf(file.path(plot_path, filename), paper = 'a4')
        per_db_region(db_regions, deseq_results)
        dev.off()
    }# }}}
}
# }}}



# main # {{{
main = function () {
    setwd("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013")
    args = list()
    args$sheet = "mm10/config/TF-steady-state-diff-bind.config"

    cnf_df = read.table(args$sheet, header=T, stringsAsFactors=F)
    proteins = levels(as.factor(cnf_df[['Factor']]))

    for (p in proteins[2]) {
        print(p)
        cnf = cnf_df %>% filter(Factor %in% c(p, keep))
        differential_binding(cnf, p)
    }

    comparisons = list.files(pattern = "[^db].tsv", path = output_path, full.name = T)
    # Keeping only the comparisons of the time points. If design is different this won't work
    comparisons = comparisons[grep('!', comparisons, invert = T)]
    comparisons = as.data.frame(comparisons)

    tmp = comparisons %>% .[[1]] %>% as.character() %>% basename() %>% file_path_sans_ext
    tmp = gsub(".*_", '', tmp) %>% strsplit("-")
    tmp = as.data.frame(do.call(rbind, tmp))
    colnames(tmp) = c('c1', 'c2')
    comparisons = cbind(comparisons, tmp)

    # Retrieve annotation from gtf files# {{{
    pattern = paste0(file_path_sans_ext(args$gtf), '.rtracklayer-')
    rtracklayer_output = c('5000-tss-2000', '2000-tss-500',
                           '5000-gene_start-2000', '2000-gene_start-500',
                           'exons', 'first-exon',
                           'introns', 'first-intron',
                           'genes', 'canonical-transcripts')
    files = paste0(paste(pattern, rtracklayer_output, sep=''), '.gtf')

    library(rtracklayer)
    if (any(file.exists(files))){
        print('Importing files.')
        annotations = lapply(files, function(x) import.gff3(x))
        names(annotations) = rtracklayer_output
        annotations$active_enhancers = import.bed(args$active_enhancers)
        annotations$poised_enhancers = import.bed(args$poised_enhancers)
    } else {
        command = paste("Rscript ~/source/get-annotation-rscript.R -G", args$gtf, '-p', ncores, collapse = " ")
        stop(paste('Annotation GTFs do not exist. Run `', command, '` first.', collapse = ' '))
    }# }}}

    # Heatmap for histone mods across time points# {{{
    for (p in proteins){
        fs = comparisons %>% dplyr::select(comparisons) %>% 
                filter(grepl(p, comparisons)) %>%
                unique() %>% .[[1]] %>% as.character()
                
        names(fs) = mutate(comparisons, label = paste(c1, c2, sep= "-")) %>%
                        filter(grepl(p, comparisons)) %>%
                        dplyr::select(label) %>% .[[1]] %>% as.character()


        bdf_cnt = aggregate_counts(fs) %>% mutate(rowname = paste(seqnames, start, end, sep = "_"))
        # format_cdf returns error if no db regions so use tryCatch to save that
        bdf = tryCatch(format_cdf(bdf_cnt),
                       error=function(e) e,
                       warning=function(w) w)
        if(is(bdf, 'error')) next;

        db_regions = with(bdf_cnt,
                          GRanges(seqnames = seqnames,
                                  IRanges(start, end),
                                  strand = strand))
        names(db_regions) = paste(seqnames(db_regions), start(db_regions), end(db_regions), sep = "_")
        x = annotate(db_regions, annotations, output_path, p)

        library(ChIPpeakAnno)
        db_distance_from_genes = annotatePeakInBatch(db_regions, AnnotationData = annotations[['genes']])
        db_distance_from_genes = as.data.frame(db_distance_from_genes)
        dup = db_distance_from_genes %>% filter(duplicated(peak)) %>% select(peak) %>% .[[1]]
        tmp_ra = add_rownames(bdf$ra, var = 'peak')
        # filtering duplicated(peak) because ChIPpeakAnno will report 2 genes if peak has equal distance from them
        # here just selecting the first one reported
        bdf$ra = db_distance_from_genes %>% select(peak, insideFeature, shortestDistance) %>%
                    rename(relative_to_gene = insideFeature,
                           d_relative_to_gene = shortestDistance) %>%
                    mutate(d_relative_to_gene = d_relative_to_gene / 1000) %>%
                    right_join(tmp_ra, by = 'peak') %>% filter(!duplicated(peak))
        rownames(bdf$ra) = bdf$ra[['peak']]
        bdf$ra = bdf$ra[,-1]

        # discretize distance to look good on annotation tracks. Doing it manually cause not sure how to get 
        # what I want with a package
        tmp = c()
        bins = c(0, 0.5, 1, 5, 10, 50, 100, 300)
        for (i in 1:length(bins)) {
            if(i == length(bins)) {
                x = y >= bins[i]
                tmp[x] = paste0(bins[i], 'kb+')
            } else {
                x = y >= bins[i] & y < bins[i+1]
                tmp[x] = paste0(bins[i], '-', bins[i+1], 'kb')
            }
        }
        bdf$ra$d_relative_to_gene = tmp

        continuous_colors = colorRampPalette(brewer.pal(9,"RdYlGn"))(length(unique(tmp)))
        bdf$ac$d_relative_to_gene = continuous_colors
        names(bdf$ac$d_relative_to_gene) = mixedsort(unique(tmp))
        
        png(file.path(plot_path, paste0('differential_binding-', p, '.png')), height = 700, width = 700)

        results = pheatmap(bdf$plot_data,
                           clustering_distance_rows = "correlation",
                           clustering_method = "ward.D2", scale = 'row',
                           cluster_cols = FALSE, cluster_rows = TRUE,
                           annotation_row = bdf$ra, annotation_col = bdf$ca,
                           annotation_legend = T, annotation_colors = bdf$ac,
                           show_rownames = FALSE, 
                           main = paste0('Diff. binding: ', p))
        get_cluster(results, output_path = output_path, label = p)
        dev.off()
    }
    #}}}
    
    expression_data = list.files(pattern = "[^de].tsv", path = args$de, full.name = T)
    expression_data = as.data.frame(expression_data)

    tmp = expression_data %>% .[[1]] %>% as.character() %>%
          basename() %>% file_path_sans_ext %>% str_replace_all("_DE|_de", "") %>% strsplit("-")
    tmp = as.data.frame(do.call(rbind, lapply(tmp, sort)))
    colnames(tmp) = c('c1', 'c2')
    expression_data = cbind(expression_data, tmp)
    
    comparisons = inner_join(comparisons, expression_data, by = c("c1", "c2"))
    
    gene_scores = edf_cnt %>% dplyr::select(gene_id, matches('FDR'))
    rownames(gene_scores) = gene_scores[[1]]
    gene_scores = gene_scores[,-1]
    g = list()
    for (x in colnames(gene_scores)){
        name = gsub("FDR.", "", x)
        g[[name]] = gene_scores[,x]
        names(g[[name]]) = rownames(gene_scores)
        g[[name]] = g[[name]][!is.na(g[[name]])]
    }

    # Doing GSEA on de genes and de and db genes
    source("~/source/Rscripts/functions.R")
    x = lapply(names(g), function(name){
        pdf(file.path(plot_path, paste0(name, '.gsea.pdf')), paper = 'a4')
        ob = get_set_enrichment(gene_scores = g[[name]], label = name)
        dev.off()
                       })

#        ob = get_set_enrichment(gene_scores = gene_scores, bound = (select(annotated_df, gene_id) %>% .[[1]]), label = label)


}

# }}}





