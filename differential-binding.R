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
cool_cols = colorRampPalette(c('steelblue','midnightblue'))(8)
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

    # anotate db sites and plot TPM for associated genes# {{{
    per_db_region = function(db_regions, deseq_results, pattern = 'Expression'){
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

        # If de files from Meryem, gene coordinates are in deseq files. # {{{
        # No need to load get_annotation(assembly) to get gene coordinates
        # chromosome names and strands needs to be changed
        if(annotation == "in_file"){
            source("~/source/Rscripts/granges-functions.R")
            genes = deseq_results %>% select(gene_id, chromosome_name, strand, start_position, end_position, gene_name)
            genes = with(genes, 
                         GRanges(seqnames = chromosome_name,
                                 IRanges(start_position, end_position),
                                 strand = strand,
                                 gene_id = gene_id,
                                 gene_name = gene_name))
            genes = change_seqnames(genes)
            names(genes) = genes$gene_id
        } else {
            #ToDO: load genes from get_annotation(assembly)
        }
        # }}}

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
                    select(Fold:gene_id, log2FoldChange, padj:ncol(.)) %>%
                    # filtering based on significance of expression change
                    filter(!is.na(padj)) %>% mutate(DE = ifelse(padj < 0.05, 'p < 0.05', 'p >= 0.05')) %>%
                    inner_join(db_coverage, by = 'peak')

        plot_data = annotated %>% select(contains(pattern),
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

    # annotate genes and plot TPM and binding affinity for all de genes# {{{
    per_de_gene = function(comparisons){

        # scaled (log() - rowMeans) expression for all genes (where NaN priot to log gene # {{{
        # expression was zero for all samples
        flag = FALSE
        for (x in (comparisons %>% select(expression_data) %>% unique() %>% .[[1]] %>% as.character()) ) {
            print(x)
            tmp = read.delim(x, head = T, sep = "\t", stringsAsFactors = F)
            if(!flag) qpcr = tmp %>% filter(gene_name %in% c('Htra1', 'Ppp2r2c', 'Cbr3', 'Bmp4')) %>%
                                     select(gene_id, gene_name)
            # Don't try to read it in and select on one go as it won't work
            tmp = tmp %>% dplyr::select(gene_id, chromosome_name:end_position,
                                        log2FoldChange, padj:ncol(.))
            colnames(tmp) = gsub('log2FoldChange',
                                 paste('log2FoldChange', file_path_sans_ext(basename(x)), sep = "."),
                                 gsub('padj',
                                      paste('padj', file_path_sans_ext(basename(x)), sep = "."),
                                      colnames(tmp)))

            if(flag) {
                print('in if')
                deseq = inner_join(deseq, tmp,
                                   by = c('gene_id', 'chromosome_name', 'strand', 'start_position', 'end_position'))
                # since samples at h0 will be in all comparisons expression columns for those will be
                # duplicated and dplyr add '.x' or '.y' to them depending from which df they came for
                # remove extension, find duplicates and remove. This would work nicer with rename() but I
                # can't get it to work on selected columns
                colnames(deseq) = gsub("\\.x|\\.y$", "", colnames(deseq))
                discard = deseq %>% names(.) %>% duplicated()
                deseq = deseq[,!discard]
            } else {
                print('in else')
                deseq = tmp
                flag = TRUE
            }
        }# }}}

        de_in_any = select(deseq, contains("padj")) %>% apply(1, function(x) min(x[!is.na(x)]) < 0.05)
        plot_data = deseq[de_in_any,]
        row_names = plot_data[['gene_id']]
        
        row_annotation = select(plot_data, contains("padj")) %>%
                         mutate_each(funs(ifelse(. > 0.05 | is.na(.), 'non significant', 'significant')))
        row_annotation = row_annotation %>% select(one_of(mixedsort(colnames(row_annotation))))
        rownames(row_annotation) = plot_data[['gene_id']]

        plot_data = plot_data %>% dplyr::select(starts_with('Expression'))
        plot_data = as.matrix(plot_data)
        rownames(plot_data) = row_names

        column_annotation = data.frame( Time = gsub(".*_", '', (gsub("h.*", 'h', colnames(plot_data)))))
        rownames(column_annotation) = colnames(plot_data)

        anno_colors = lapply(1:ncol(row_annotation), function(x) c("non significant" = "black",
                                                                   "significant" = "firebrick"))
        names(anno_colors) = colnames(row_annotation)
        tmp = cool_cols[1:length(unique(column_annotation[['Time']]))]
        names(tmp) = mixedsort(unique(column_annotation[['Time']]))
        anno_colors$Time = tmp

        pdf(file.path(plot_path, 'expression.pdf'), paper = 'a4')# {{{
        results = pheatmap(plot_data,
                           clustering_distance_rows = "correlation", clustering_method = "ward.D2", scale = 'row',
                           cluster_cols = TRUE, cluster_rows = TRUE, cutree_rows = 10,
                           annotation_row = row_annotation, annotation_col = column_annotation,
                           annotation_legend = T, annotation_colors = anno_colors,
                           show_rownames = FALSE)

        # k must be equal to cutree_rows used in pheatmap for defining number of clusters
        # http://stackoverflow.com/questions/27820158/pheatmap-in-r-how-to-get-clusters
        drows = as.dist(1-cor(t(plot_data), method = "pearson"))
        clusters = cutree(results$tree_row, k = 10)
        heatmap_order = names(clusters[results$tree_row$order])
        
        qpcr_genes_expression = plot_data[rownames(plot_data) %in% qpcr[["gene_id"]],]
        qpcr_results = pheatmap(qpcr_genes_expression,
                                clustering_distance_rows = "correlation", clustering_method = "ward.D2", scale = 'row',
                                cluster_cols = TRUE, cluster_rows = TRUE, cutree_rows = 2,
                                annotation_row = row_annotation[rownames(row_annotation) %in% qpcr[["gene_id"]],],
                                annotation_col = column_annotation,
                                annotation_legend = T, annotation_colors = anno_colors,
                                show_rownames = TRUE, labels_row = qpcr[['gene_name']])
        qpcr_drows = as.dist(1-cor(t(qpcr_genes_expression), method = "pearson"))
        qpcr_clusters = cutree(qpcr_results$tree_row, k = 2)
        qpcr_heatmap_order = names(qpcr_clusters[qpcr_results$tree_row$order])

        dev.off()# }}}

        # If de files from Meryem, gene coordinates are in deseq files. # {{{
        # No need to load get_annotation(assembly) to get gene coordinates
        # chromosome names and strands needs to be changed
        if(annotation == "in_file"){
            source("~/source/Rscripts/granges-functions.R")
            genes = deseq %>% select(gene_id:end_position)
            genes = with(genes,
                         GRanges(seqnames = chromosome_name,
                                 IRanges(start_position, end_position),
                                 strand = strand,
                                 gene_id = gene_id))
            genes = change_seqnames(genes)
            names(genes) = genes$gene_id
        } else {
            #ToDO: load genes from get_annotation(assembly)
        }
        # }}}

        # read_diffbind # {{{
        read_diffbind = function(comparisons, protein){
            flag = FALSE
            files = comparisons %>% filter(grepl(protein, comparisons)) %>% select(comparisons) %>%
                    unique() %>% .[[1]] %>% as.character()
            for (x in files) {
                print(x)
                tmp = read.delim(x, head = T, sep = "\t", stringsAsFactors = F)
                # Don't try to read it in and select on one go as it won't work
                tmp = tmp %>% dplyr::select(seqnames:strand,
                                            Fold, FDR:ncol(.))
                colnames(tmp) = gsub('Fold',
                                     paste('log2FoldChange', file_path_sans_ext(basename(x)), sep = "."),
                                     gsub('FDR',
                                          paste('padj', file_path_sans_ext(basename(x)), sep = "."),
                                          colnames(tmp)))

                if(flag) {
                    print('in if')
                    diffbind = inner_join(diffbind, tmp,
                                       by = c('seqnames', 'start', 'end', 'width', 'strand'))
                    # since samples at h0 will be in all comparisons expression columns for those will be
                    # duplicated and dplyr add '.x' or '.y' to them depending from which df they came for
                    # remove extension, find duplicates and remove. This would work nicer with rename() but I
                    # can't get it to work on selected columns
                    colnames(diffbind) = gsub("\\.x|\\.y$", "", colnames(diffbind))
                    discard = diffbind %>% names(.) %>% duplicated()
                    diffbind = diffbind[,!discard]
                } else {
                    print('in else')
                    diffbind = tmp
                    flag = TRUE
                }
            }
            return(diffbind)
        }# }}}
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
    }
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

}# }}}
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
    test = function(){
        db_regions = as.character(comparisons[2,1])
        deseq_results = as.character(comparisons[2,4])
        combine_expression(db_regions, deseq_results)
    }


}

# }}}

efs = comparisons %>% select(expression_data) %>%
      unique() %>% .[[1]] %>% as.character()
names(efs) = mutate(comparisons, label = paste(c1, c2, sep= "-")) %>%
             select(expression_data, label) %>% unique() %>%
             select(label) %>% .[[1]] %>% as.character()

fs = comparisons %>% select(comparisons) %>% 
     filter(grepl(proteins[1], comparisons)) %>%
     unique() %>% .[[1]] %>% as.character()
names(fs) = mutate(comparisons, label = paste(c1, c2, sep= "-")) %>%
            filter(grepl(proteins[1], comparisons)) %>%
            select(label) %>% .[[1]] %>% as.character()
#fs = efs
#type = "expression"

bdf = aggregate_count_data(fs)
pdf(file.path(plot_path, paste0('test-', proteins[1], '.pdf')), paper = 'a4')
    pheatmap(plot_data_db,
             clustering_distance_rows = drows,
             clustering_method = "ward.D2", scale = 'row', cutree_rows = 10,
             cluster_cols = FALSE, cluster_rows = TRUE)
             #annotation_row = row_annotation, annotation_col = column_annotation,
             #annotation_legend = T, annotation_colors = anno_colors,
             #show_rownames = FALSE)
dev.off()

aggregate_count_data = function(fs, type = "binding"){# {{{
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

# {{{
#setwd("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/macs/")

#SampleSheet must have SampleID,Condition,Treatment,Replicate,bamReads,bamControl,Peaks (last 3 are the path to the respective bed files)
DB<-function(sheet, fdr=0.05){
   library(DiffBind)

   output <- gsub(".csv", 'pdf', gsub(".*/", '', sheet))

   #Cause it will otherwise complain and terminate when you run it on a cluster
   pdf(output)
   #working object is of dba class

   #Check how many Conditions are there
   line <- unlist(strsplit(tolower(readLines(sheet[1], n=1)), ','))
   if(!"condition" %in% line) stop("Sample Sheet does not contain a Condition column")

   protein <- dba(sampleSheet= sheet)
   #Will print the correlation heatmap using occupancy(peak caller score) data
   plot(protein)

   sheet_table <-read.table(sheet, sep=',')
   conditions <- sort(levels(sheet_table$condition), decreasing=T)

   nCond <- match("condition", line)
   tmp <- rep('NULL', length(line))
   tmp[nCond] <- NA
   nCond <- sort(levels(read.table(sheet, colClasses=tmp, sep=",", head=T)[[1]]), decreasing=T)

   #calculates a binding matrix, scores based on read counts (affinity scores) rather than the previous one that plotted confidence scores.
   #DBA_SCORE_READS_MINUS => read count for interval-read count for interval in control..might be a better way of doing this
   #bScaleControl=TRUE => scale the control reads based on relative library size
   wobj <- dba.count(wobj, score=DBA_SCORE_READS_MINUS, bScaleControl=TRUE, bParallel=TRUE)
   btag <- TRUE

   if(sum(wobj$masks[[nCond[1]]]) < 2 | sum(wobj$masks[[nCond[2]]]) < 2){
      wobj <- dba.contrast(wobj, wobj$masks[[nCond[1]]], wobj$masks[[nCond[2]]], nCond[1], nCond[2])
      btag <- FALSE
   }
   else{
      #step is actually optional cause it's set up
      wobj <- dba.contrast(wobj, categories=DBA_CONDITION)
   }


   wobj <- dba.analyze(
                       wobj,
                       method= DBA_DESEQ,
                       bSubControl= TRUE,
                       #bFullLibrarySize= TRUE,
                       #SOS if no replicates this option should be set to FALSE
                       bTagwise= btag,
                       bParallel= TRUE)

   wobj.DB <- dba.report(wobj,
                         method=DBA_DESEQ,
                         #only sites with an absolute Fold value greater than equal to this will be included in the report
                         #fold= 2,
                         th= fdr,
                         initString= fdr,
                         file= paste(gsub(".csv", '', gsub(".*/", '', sheet)), '.diffbind', sep=''))
   wobj.df <- dba.report(wobj, method=DBA_DESEQ, th=fdr, DataType = DBA_DATA_FRAME)


   write.table(wobj.df, file=gsub(".csv", 'pdf', gsub(".*/", '', sheet)), sep="\t", quote=FALSE, row.names=FALSE)
   #plot(wobj)
   #dba.plotMA(wobj)
   #PCA based on the afinity data for all sites
   #dba.plotPCA(wobj, DBA_CONDITION)
   #PCA based on the affiity data for the DB sites
   #dba.plotPCA(wobj, contrast=1, th=.05)
   dev.off()
   return(list("db"=wobj.DB))
}

#myPeakFiles <- list.files(path=getwd(), full.name=T )
#samplesheets must point to modified macs peaks.bed with 4 columns instead of five
mySampleSheets <- list.files("..", full.name=T, pattern="csv")

#system.time(DB(mySampleSheets[1]), .01)


#Comment in and out depending what you are running
#Mi2b_2i <- DB(mySampleSheets[1], 0.05)
#Mi2b_Epi <- DB(mySampleSheets[2], 0.05)
#Mi2b_serum <- DB(mySampleSheets[3], 0.05)



#source("filewith get.annotation")
#mm9.annot<-get.annotation('may2012.archive.ensembl.org', 'mmusculus_gene_ensembl')




# }}}
