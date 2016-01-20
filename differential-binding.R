#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library = function (...) suppressMessages(base::library(...))
select = dplyr::select
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")

parser =  ArgumentParser(description="Perform differential binding analysis")
parser$add_argument('-s', '--sheet', metavar= "file", required='True', type= "character", help= "Sample sheet must have SampleID,Condition,Treatment,Replicate,bamReads,bamControl,Peaks (last 3 are the path to the respective bed files)")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-k', '--keep', type= "character", default='', help= "Factors to ignore for differential binding but plot")
parser$add_argument('-d', '--de', metavar = "path", type= "character", default='', help= "Path to look for the list of de expressed genes")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args = parser$parse_args()

output_path = file.path(args$out, 'DiffBind')
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)
#}}}

# Load Packages
library(DiffBind)
library(dplyr)
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

# doing all the work # {{{
differential_binding = function (cnf, protein){
    tmp = import_data(cnf)
    config = tmp[['cnf_df']]
    peakset = tmp[['peakset']]
    
    # Using only this peak caller data, a correlation heatmap can be generated 
    # which gives an initial clustering of the samples
    pdf(file.path(plot_path, paste0(protein, '.pdf')))
    plot(peakset)
    peakset = dba.peakset(peakset, consensus = -DBA_REPLICATE, minOverlap = 2)

    counts = dba.count(peakset, score = DBA_SCORE_TMM_MINUS_FULL,
#    counts = dba.count(peakset, score = DBA_SCORE_READS_MINUS,
                       minOverlap = 2,
#                       peaks = peakset$masks$Consensus,
#                       bRemoveDuplicates = TRUE,
                       bScaleControl = TRUE, bParallel = TRUE,
                       bCorPlot = FALSE)
    plot(counts)
    dev.off()

    contrast = dba.contrast(counts, categories = DBA_CONDITION, minMembers = 2)
    
    diff_bound = dba.analyze(contrast, method = DBA_DESEQ2,
                             bSubControl= TRUE, bFullLibrarySize = TRUE,
                             bTagwise = TRUE,
                             bParallel = TRUE,
                             bCorPlot = FALSE)

    for (index in 1:length(contrast$contrasts)){
        comparison = paste(protein,
                           contrast$contrasts[[index]]$name1,
                           contrast$contrasts[[index]]$name2, sep='-')
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
        
        results =  dba.report(diff_bound, method = DBA_DESEQ2,
                              contrast = index,
                              # add count data for individual samples
                              bCounts = TRUE,
                              # only include sites with an absolute Fold value greater than equal
                              # fold= 2
                              # Threshold will be FDR and not pval
                              bUsePval = FALSE,
                              th = 0.05)

        write.table(as.data.frame(results),
                    file = file.path(output_path, paste0(comparison, '-db.tsv')), sep="\t", quote=FALSE, row.names=FALSE)

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
                    file = file.path(output_path, paste0(comparison, '.tsv')), sep="\t", quote=FALSE, row.names=FALSE)

    }
}
# }}}

# if working with meryem files just give one of the files (not the -de.tsv though)# {{{
retrieve_annotation = function(fs){

    # If de files from Meryem, gene coordinates are in deseq files. # {{{
    # No need to load get_annotation(assembly) to get gene coordinates
    # chromosome names and strands needs to be changed
    if(!missing(fs)){
        fs = read.delim(fs, head = T, sep = "\t", stringsAsFactors = F)
        source("~/source/Rscripts/granges-functions.R")
        genes = with((fs %>% dplyr::select(chromosome_name:end_position)),
                     GRanges(seqnames = chromosome_name,
                             IRanges(start_position, end_position),
                             strand = strand))
        values(genes) = (fs %>% dplyr::select(starts_with("gene")))
        names(genes) = genes$gene_id
        genes = change_seqnames(genes)
    } else {
        #ToDO: load genes from get_annotation(assembly)
    }
    # }}}
    return(genes)
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

# combine db regions with expression data for nearest gene, annotation = ensembl, jan2013 etc if not in_file # {{{
combine_expression = function(comparisons, annotation = "in_file"){
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


# Prepare a data.frame containing counts/TMM for multiple conditions/samples for plotting# {{{
# ca is the name for column annotation
format_cdf = function(cdf, statistic = "FDR", ca = "Time", threshold = 0.05){
    # Significant changes between any samples in the cdf
    significant = dplyr::select(cdf, contains(statistic)) %>% apply(1, function(x) min(x[!is.na(x)]) < threshold)
    plot_data = cdf[significant,]
    # Save because it is lost with dplyr
    row_names = plot_data[['rowname']]
    
    # Row annotation needs to be defined before selecting columns
    row_annotation = dplyr::select(plot_data, contains(statistic)) %>%
                     mutate_each(funs(ifelse(. > 0.05 | is.na(.), 'non significant', 'significant')))
    row_annotation = row_annotation %>% dplyr::select(one_of(rev(mixedsort(colnames(row_annotation)))))
    rownames(row_annotation) = row_names

    plot_data = plot_data %>% dplyr::select(-starts_with('log2FoldChange'),
                                            -starts_with(statistic),
                                            -seqnames, -strand, -start, -end,
                                            -rowname)
    plot_data = plot_data %>% dplyr::select(match(mixedsort(colnames(plot_data)), colnames(plot_data)))

    plot_data = as.matrix(plot_data)
    rownames(plot_data) = row_names

    # Column names should be sample_rep with no '_' in "sample"
    column_annotation = data.frame( x = gsub("_.*", '', colnames(plot_data)))
    colnames(column_annotation) = ca
    rownames(column_annotation) = colnames(plot_data)

    source("~/source/Rscripts/ggplot-functions.R")
    color = c('slategrey', 'violetred')
    anno_colors = lapply(1:ncol(row_annotation), function(x) c("non significant" = color[1],
                                                               "significant" = color[2]))
    names(anno_colors) = colnames(row_annotation)

    cool_colors = colorRampPalette(brewer.pal(9,"YlGnBu"))(8)
    tmp = cool_colors[1:length(unique(column_annotation[[ca]]))]
    names(tmp) = mixedsort(unique(column_annotation[[ca]]))
    anno_colors[[ca]] = tmp

    return(list(plot_data = plot_data,
                ra = row_annotation, ca = column_annotation,
                ac = anno_colors))
} # }}}

aggregate_count_data = function(fs, type = "binding"){ # {{{
    flag = FALSE
    for (x in 1:length(fs)) {
        tmp = read.delim(fs[x], head = T, sep = "\t", stringsAsFactors = F)

        rename_map = c(paste("log2FoldChange", names(fs)[x], sep = "."),
                       paste("FDR", names(fs)[x], sep = "."))
        combine_by = c('seqnames', 'strand', 'start', 'end')

        if(type == "expression") {# {{{
            # Don't try to read it in and select on one go as it won't work
            tmp = tmp %>% dplyr::select(gene_id, chromosome_name:end_position,
                                        log2FoldChange, padj:ncol(.)) %>%
                  rename_(.dots = setNames(list("log2FoldChange", "padj",
                                                "chromosome_name", "start_position", "end_position"),
                          c(rename_map, "seqnames", "start", "end")))

            combine_by = c('gene_id', combine_by)
        } else {
            tmp = tmp %>% dplyr::select(seqnames:end, strand, Fold, FDR:ncol(.)) %>%
                  rename_(.dots = setNames(list("Fold", "FDR"), rename_map))
        }
        # }}}

        if(flag) {# {{{
            ofs = inner_join(ofs, tmp, by = combine_by)
            # since samples at h0 will be in all comparisons expression columns for those will be
            # duplicated and dplyr add '.x' or '.y' to them depending from which df they came for
            # remove extension, find duplicates and remove. This would work nicer with rename() but I
            # can't get it to work on selected columns
            colnames(ofs) = gsub("\\.x|\\.y$", "", colnames(ofs))
            discard = ofs %>% names(.) %>% duplicated()
            ofs = ofs[,!discard]
        } else {
            ofs = tmp
            flag = TRUE
        }
    }# }}}
    return(ofs)
}
# }}}

# When wanting to retrieve a cluster call the function to get IDs in particular cluster# {{{
# df = tmp$plot_data, regions = bdf
get_cluster = function(df, regions, results, total_k, cluster, threshold = 0.05, label){
    # k must be equal to cutree_rows used in pheatmap for defining number of clusters
    # http://stackoverflow.com/questions/27820158/pheatmap-in-r-how-to-get-clusters
    drows = as.dist(1-cor(t(df), method = "pearson"))
    clusters = cutree(results$tree_row, k = total_k)
    heatmap_order = names(clusters[results$tree_row$order])

    keep = clusters[clusters == cluster]
    # Create df with just FDR scores for all regions in cluster
    filters = bdf %>% filter(rowname %in% names(keep)) %>% select(starts_with("FDR"))
    rownames(filters) = bdf %>% filter(rowname %in% names(keep)) %>% select(rowname) %>% .[[1]]
    keep = keep[apply(filters < threshold, 1, all)]

    cluster = regions %>% filter(rowname %in% names(keep))

    cluster_gr = with(cluster,
                       GRanges(seqnames = seqnames,
                               IRanges(start, end),
                               strand = strand))

    values(cluster_gr) = (cluster %>% dplyr::select(-seqnames, -start, -end, -strand))

    annotated = annotatePeakInBatch(cluster_gr, AnnotationData = genes)
    annotated = as.data.frame(annotated)
    minimal = as.data.frame(values(genes)) %>% rename(feature = gene_id)

    annotated = inner_join(annotated, minimal, by = "feature")
    write.table(annotated, file = file.path(output_path, paste0(label,".tsv")),
                sep = "\t", col.names = T, row.names = F, quote = F)
    pdf(file.path(plot_path, paste0('distances-', label, '.pdf')), paper = 'a4')
    hist(annotated$shortestDistance)
    dev.off()
}
# }}}

# main # {{{
main = function () {
#    setwd("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013")
#    args = list()
#    args$sheet = list.files(pattern="diff-bind.config", full.names=T)

    cnf_df = read.table(args$sheet, header=T, stringsAsFactors=F)
    proteins = levels(as.factor(cnf_df[['Factor']]))
    keep = unlist(strsplit(args$keep, ':'))
    for (x in keep){
        proteins = proteins[grep(x, proteins, invert = TRUE)]
    }

    for (p in proteins) {
        print(p)
        cnf = cnf_df %>% filter(Factor %in% c(p, keep))
        differential_binding(cnf, p)
    }

    comparisons = list.files(pattern = "[^db].tsv", path = output_path, full.name = T)
    # Keeping only the comparisons of the time points. If design is different this won't work
    # comparisons = comparisons[grep('h.-h.', comparisons)]
    comparisons = comparisons[grep('!', comparisons, invert = T)]
    comparisons = as.data.frame(comparisons)

    tmp = comparisons %>% .[[1]] %>% as.character() %>% basename() %>% file_path_sans_ext %>% strsplit("-")
    tmp = as.data.frame(do.call(rbind, tmp)) %>% select(V2, V3)
    tmp = lapply(1:nrow(tmp), function(x) sort(tmp[x,]))
    tmp = as.data.frame(do.call(rbind, tmp))
    colnames(tmp) = c('c1', 'c2')
    comparisons = cbind(comparisons, tmp)

    expression_data = list.files(pattern = "[^de].tsv", path = args$de, full.name = T)
    expression_data = as.data.frame(expression_data)

    tmp = expression_data %>% .[[1]] %>% as.character() %>%
          basename() %>% file_path_sans_ext %>% str_replace_all("_DE|_de", "") %>% strsplit("-")
    tmp = as.data.frame(do.call(rbind, lapply(tmp, sort)))
    colnames(tmp) = c('c1', 'c2')
    expression_data = cbind(expression_data, tmp)
    
    comparisons = inner_join(comparisons, expression_data, by = c("c1", "c2"))
    genes = retrieve_annotation(as.character(comparisons[1,'expression_data']))

    # Create matrix of expression for de genes only
    # Do not plot just yet. Use edf$ra (row annotation) on the heatmap of peaks associated with a gene
    efs = comparisons %>% dplyr::select(expression_data) %>% unique() %>% .[[1]] %>% as.character()
    names(efs) = mutate(comparisons, label = paste(c1, c2, sep= "-")) %>%
                 dplyr::select(expression_data, label) %>% unique() %>%
                 dplyr::select(label) %>% .[[1]] %>% as.character()
    edf_cnt = aggregate_count_data(efs, type = "expression")
    edf = edf_cnt %>% rename(rowname = gene_id)
    rename_columns = gsub('_2i', '',
                          gsub('h_', '_',
                               gsub(".*Rescue_", 'h', colnames(edf))))
    colnames(edf) = rename_columns
    edf = format_cdf(edf)

    # Heatmap for histone mods across time points# {{{
    for (p in proteins){
        fs = comparisons %>% dplyr::select(comparisons) %>% 
                filter(grepl(p, comparisons)) %>%
                unique() %>% .[[1]] %>% as.character()
                
        names(fs) = mutate(comparisons, label = paste(c1, c2, sep= "-")) %>%
                        filter(grepl(p, comparisons)) %>%
                        dplyr::select(label) %>% .[[1]] %>% as.character()

        bdf_cnt = aggregate_count_data(fs) %>% mutate(rowname = paste(seqnames, start, end, sep = "_"))
        bdf = format_cdf(bdf_cnt)

        db_regions = with(bdf_cnt,
                          GRanges(seqnames = seqnames,
                                  IRanges(start, end),
                                  strand = strand))
        names(db_regions) = paste(seqnames(db_regions), start(db_regions), end(db_regions), sep = "_")

        library(ChIPpeakAnno)
        library(rtracklayer)
#        db_annotated = annotatePeakInBatch(db_regions, AnnotationData = genes)
#        db_annotated = as.data.frame(db_annotated) %>% select(peak, feature, insideFeature, shortestDistance)
        gene_promoters = promoters(genes, upstream = 2000, downstream = 500)
        ov = findOverlaps(db_regions, gene_promoters, select = "arbitrary")
        df = data.frame(region = names(db_regions))
        df = left_join(df,
                        data.frame(region = names(db_regions)[!is.na(ov)],
                                   promoter = names(gene_promoters)[ov[!is.na(ov)]]),
                        by = "region")
        active_enhancers = import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-active_2015-08-10.bed")
        names(active_enhancers) = paste(seqnames(active_enhancers), start(active_enhancers), end(active_enhancers), sep="_")
        ov = findOverlaps(db_regions, active_enhancers, select = "arbitrary")
        df = left_join(df, 
                       data.frame(region = names(db_regions)[!is.na(ov)],
                                  active_enhancers = names(active_enhancers)[ov[!is.na(ov)]]),
                       by = "region")
        poised_enhancers = import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-poised_2015-08-10.bed")
        poised_enhancers = poised_enhancers[queryHits(findOverlaps(poised_enhancers, active_enhancers))]
        names(poised_enhancers) = paste(seqnames(poised_enhancers), start(poised_enhancers), end(poised_enhancers), sep="_")
        ov = findOverlaps(db_regions, poised_enhancers, select = "arbitrary")
        df = left_join(df, 
                       data.frame(region = names(db_regions)[!is.na(ov)],
                                  poised_enhancers = names(poised_enhancers)[ov[!is.na(ov)]]),
                       by = "region")
        ov = findOverlaps(db_regions, genes, select = "arbitrary")
        df = left_join(df, 
                       data.frame(region = names(db_regions)[!is.na(ov)],
                                  genes = names(genes)[ov[!is.na(ov)]]),
                       by = "region")
        df = df %>% mutate(status = ifelse((is.na(promoter) & is.na(genes) & is.na(active_enhancers) & is.na(poised_enhancers))
                                            , 'intergenic',
                                            ifelse(!is.na(active_enhancers),
                                                   'active_enhancers',
                                                   ifelse(!is.na(poised_enhancers),
                                                          'poised_enhancers',
                                                          ifelse(!is.na(promoter),
                                                          'promoter','genes')
                                                          )
                                                   )
                                            ))

        # when I want to find the percentage of overlap
#        x = db_regions[queryHits(ov)]
#        y = active_enhancers[subjectHits(ov)]
#        table(width(pintersect(x, y))/width(x) < 0.01)
#        table(width(pintersect(x, y))/width(y) < 0.)

        if(FALSE){# {{{
            annotated = annotatePeakInBatch(genes, AnnotationData = db_regions)
            annotated = as.data.frame(annotated)
            # since annotating genes with peaks (and not the other way round) 
            # gene identifiers will be under column peaks, select genes d.e.,
            # reorder based on clusters from gene expression heatmap and subset db_regions
            # downstream of gene, overlapStart of gene i.e. TSS, includeFeature = include peak
            a = annotated %>% dplyr::select(one_of('gene_id', 'feature', 'insideFeature', 'shortestDistance'))

            b = add_rownames(bdf$ra, var = "feature")
            ra = left_join(b, a, by='feature')
            row_annotation = add_rownames(edf$ra, var = 'gene_id')
            colnames(row_annotation) = gsub('FDR.', 'Expression.FDR', colnames(row_annotation))
            ra = left_join(ra, row_annotation, by="gene_id")

    #        a = right_join(a,row_annotation, by='gene_id') %>% dplyr::select(-gene_id)
    #        b = add_rownames(bdf$ra, var = "feature")
    #        ra = left_join(b, a, by='feature')

            if(any(duplicated(ra$feature))) {
                dups = ra %>% filter(duplicated(feature)) %>% dplyr::select(feature) %>% .[[1]]
                keep = lapply((ra %>% filter(feature %in% dups) %>% group_by(feature) %>% do(x = as_data_frame(.)) %>% .$x),
                       function(x){
                           k = x %>% dplyr::select(ncol(.):2) %>% mutate(x= min_rank(shortestDistance)) %>%
                                 filter(x ==1) %>% dplyr::select(gene_id) %>% .[[1]]
                           k[1]
                       })
                ra = ra %>% filter(!((feature %in% dups) & !(gene_id %in% unlist(keep))))
            }
            levels(ra[["insideFeature"]]) = c(levels(ra[["insideFeature"]]), 'non significant')
            ra = ra %>% select(-shortestDistance) %>% mutate_each(funs(replace(., which(is.na(.)), "non significant"))) %>%
                    dplyr::select(-gene_id)
            
            rownames(ra) = ra[['feature']]
            ra = ra[,-1]
            color = gg_color_hue(2)
            ac = lapply(1:length(grep('Expression', colnames(ra))),
                        function(x) c("non significant" = color[1],
                                      "significant" = color[2]))
            names(ac) = colnames(ra)[grep('Expression', colnames(ra))]

            ac$insideFeature = gg_color_hue(length(levels(ra$insideFeature)))
            names(ac$insideFeature) = levels(ra$insideFeature)

            ac = c(ac, bdf$ac)
        }

        ra = left_join(add_rownames(bdf$ra, var = "region"), select(df, region, status), by = "region")
        rownames(ra) = ra[['region']]
        ra = ra[,-1]
        ac = list()
        set1_col = brewer.pal(9, 'Set1')
        ac$status = c(set1_col[2], set1_col[5], set1_col[9], set1_col[4], set1_col[1])
        names(ac$status) = levels(factor(df$status))
        ac = c(ac, bdf$ac)


        png(file.path(plot_path, paste0('differential_binding-', p, '.png')), height = 700, width = 700)

        results = pheatmap(bdf$plot_data,
                           clustering_distance_rows = "correlation",
                           clustering_method = "ward.D2", scale = 'row', cutree_rows = 5,
                           cluster_cols = FALSE, cluster_rows = TRUE,
                           annotation_row = ra, annotation_col = bdf$ca,
                           annotation_legend = T, annotation_colors = ac,
                           show_rownames = FALSE, 
                           main = paste0('Diff. binding: ', p))
        dev.off()
        rm(ac)
    }# }}}# }}}

    # Heatmap for expression across time points # {{{
    png(file.path(plot_path, 'expression.png'), height=700, width =700)

    results = pheatmap(edf$plot_data,
                       clustering_distance_rows = "correlation",
                       clustering_method = "ward.D2", scale = 'row', cutree_rows = 5,
                       cluster_cols = FALSE, cluster_rows = TRUE,
                       annotation_row = edf$ra, annotation_col = edf$ca,
                       annotation_legend = T, annotation_colors = edf$ac,
                       show_rownames = FALSE, 
                       main = "Gene expression across time points")
    dev.off()
    drows = as.dist(1-cor(t(edf$plot_data), method = "pearson"))
    clusters = cutree(results$tree_row, k = 5)
    heatmap_order = names(clusters[results$tree_row$order])
    heatmap_order_df = left_join((as.data.frame(heatmap_order) %>% rename(gene_id = heatmap_order)),
                                 add_rownames(as.data.frame(clusters), var = 'gene_id'), by = 'gene_id')
    a = read.delim("/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/deseq//h1-h0.tsv")
    a = a %>% dplyr::select(gene_id) %>% mutate(clusters = 0) %>% filter(!gene_id %in% heatmap_order)
    heatmap_order_df = rbind(a, heatmap_order_df)
    write.table(heatmap_order_df, file.path(output_path, 'expression_heatmap_order.tsv'),
                sep = "\t", col.names = T, row.names = F, quote = F)

    steady_state_de = structure("/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/deseq/2i_wt-ko.tsv", names = 'steady_state')
    steady_state_counts = structure("/nfs/research2/bertone/user/ralser/RNAseq/Brian/expression/20150811_corr_Expression_SL_and_2i", names = 'steady_state')
    steady_state_counts = read.delim(steady_state_counts, head = T, sep = "\t", stringsAsFactors = F)

    steady_state_counts = steady_state_counts %>% dplyr::select(gene_id, matches("2i"))
    rownames(steady_state_counts) = steady_state_counts[['gene_id']]
    steady_state_counts = steady_state_counts[,-1]
    
    pdf(file.path(plot_path, 'steady_state_expression.pdf'))
    pheatmap(steady_state_counts[rownames(edf$plot_data),],
             clustering_distance_rows = drows,
             clustering_method = "ward.D2", scale = 'row', cutree_rows = 5,
             cluster_cols = TRUE, cluster_rows = TRUE,
#             annotation_row = row_annotation[qpcr_keep,], annotation_col = column_annotation,
#             annotation_legend = T, annotation_colors = anno_colors,
             show_rownames = FALSE)
    dev.off()
    gtf = import.gff("/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.70.gtf")
    gene_names = as.data.frame(values(gtf))[['group']]
    a = parallel::mclapply(gene_names, function(x)
                           gsub(".* ", "", gsub("\"", "", unlist(strsplit(gsub(';', ':', x), ":"))[1])))
    a = unlist(a)
    names(gtf) = a
    # 3 genes were not in my annotation
    keep = heatmap_order[heatmap_order %in% names(gtf)]
    write.table((as.data.frame(keep) %>% rename(gene_id = keep)), file.path('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/results/profiles', 'expression_heatmap_order_list.tsv'),
                sep = "\t", col.names = T, row.names = F, quote = F)
    export(gtf[keep],"/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/results/profiles/Mus_musculus.GRCm38.70.DExpression_set.gtf", "GFF")

    # make a test.gtf to see that ordering in python works
    i = length(keep)-5
    j = length(keep)
    write.table((as.data.frame(keep) %>% rename(gene_id = keep)), file.path('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/results/profiles', 'test_list.tsv'),
                sep = "\t", col.names = T, row.names = F, quote = F)
    export(gtf[keep[c(1:5,i:j)]],"/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/results/profiles/test.gtf", "GFF")
    # }}}
    
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





