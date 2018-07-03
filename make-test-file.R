#!/usr/bin/env Rscript

# annotation plots # {{{
a = read.delim('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/Condition_7E12_2i_Proteins_M2_meryem_filtered_AND_Chd4_meryem_filtered_Ov1bp_complete_annotated.annotated.tsv', head=T)

plot_data = as.data.frame(a) %>% select(-Regions, -File, -exons, -first.exon, -first.intron, -introns) %>% rowwise() %>% mutate(gene_body = 1 - rowSums(.)) %>% t()
plot_data = as.data.frame(plot_data) %>% add_rownames(var = 'Feature')
colnames(plot_data) = c('Feature', c('Chd4 only', 'NuRD', 'Mbd3 only'))
plot_data = reshape2::melt(plot_data) %>% mutate(Feature = ifelse(Feature == 'X2000.tss.500', '2000-tss-500',Feature))


ordering = c('2000-tss-500', 'gene_body', 'intergenic', 'poised_enhancers', 'active_enhancers')
plot_data$Feature = factor(plot_data$Feature, levels = ordering)
plot_data = plot_data %>% arrange(Feature)

pdf('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/plots/NuRD_nice_col.pdf')
p = ggplot(plot_data, aes(x = variable, y = value, fill = Feature)) + geom_bar(stat="identity")
p = p + scale_fill_brewer(palette = "Set2")
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + ylab('% of peaks') + xlab('')
p
dev.off()

b = read.delim('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/Condition_2i_AND_EpiSC_Proteins_NuRD_Ov1bp_complete_annotated.annotated.tsv', head=T)
plot_data = as.data.frame(b) %>% select(-Regions, -File, -exons, -first.exon, -first.intron, -introns) %>% rowwise() %>% mutate(gene_body = 1 - rowSums(.)) %>% t()
plot_data = as.data.frame(plot_data) %>% add_rownames(var = 'Feature')
colnames(plot_data) = c('Feature', c('ESCs', 'ESCs-EpiSCs', 'EpiSCs'))
plot_data = reshape2::melt(plot_data) %>% mutate(Feature = ifelse(Feature == 'X2000.tss.500', '2000-tss-500',Feature))

plot_data = mutate(plot_data, Feature = ifelse(Feature %in% c('poised_enhancers', 'active_enhancers'), paste0('ESCs_', Feature), Feature))
ordering = c('2000-tss-500', 'gene_body', 'intergenic', 'ESCs_poised_enhancers', 'ESCs_active_enhancers')
plot_data$Feature = factor(plot_data$Feature, levels = ordering)
plot_data = plot_data %>% arrange(Feature)

pdf('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/plots/NuRD_nice_col-ESCvsEpiSC.pdf')
p = ggplot(plot_data, aes(x = variable, y = value, fill = Feature)) + geom_bar(stat="identity")
p = p + scale_fill_brewer(palette = "Set2")
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(angle = 25, hjust = 1)) + ylab('% of NuRD peaks') + xlab('')
p
dev.off()# }}}

esc_targets = read.delim('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/Condition_7E12_2i_Proteins_M2_meryem_filtered_AND_Chd4_meryem_filtered_Ov1bp.binding-per-gene.tsv')# {{{
episc_targets = read.delim('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/Condition_7E12_EpiSC_Proteins_M2_AND_Chd4_Ov1bp.binding-per-gene.tsv')

esc_regulated = read.delim("/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10//deseq/0+2i_wt-ko.deseq-de.tsv") %>% select(Gene) %>% mutate(ESCs_regulated = Gene)
episc_regulated = read.delim("/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10//deseq/episc_wt-ko.deseq-de.tsv") %>% select(Gene) %>% mutate(EpiSCs_regulated = Gene)
tmp_reg = full_join(esc_regulated, episc_regulated, by = "Gene")

esc_targets = esc_targets %>% filter(Bound == '2000-tss-500') %>% select(Gene) %>% mutate(ESCs_bound = Gene)
episc_targets = episc_targets %>% filter(Bound == '2000-tss-500') %>% select(Gene) %>% mutate(EpiSCs_bound = Gene)

tmp_targets = full_join(esc_targets, episc_targets, by = "Gene")
NuRD_targeted_genes = full_join(tmp_reg, tmp_targets, by = "Gene") %>% select(ESCs_bound, ESCs_regulated, EpiSCs_bound, EpiSCs_regulated, -Gene)
sets = colnames(NuRD_targeted_genes)
x = lapply(1:length(sets), function(x) c(0,1))
names(x) = sets
mat = expand.grid(x)
mat$Counts = 0
for(i in 2:nrow(mat)){
    current_columns = colnames(mat)[mat[i,] == 1]
    if(length(current_columns) == 1){
        tmp = as.character(NuRD_targeted_genes[, current_columns, drop = T])
        mat[i, 'Counts'] = length(tmp[!is.na(tmp)])
    } else {
        mat[i, 'Counts'] = NuRD_targeted_genes[, current_columns] %>% na.omit() %>% nrow()
    }
}
library(UpSetR)
upsetR_cnt = do.call(rbind, apply(mat, 1, function(x) {
                                      if(x['Counts'] != 0){
                                          y = as.data.frame(t(x[-length(x)]))
                                          cnt = x['Counts']
                                          y[1:cnt,] = y[1,]
                                          return(y)
                                      }
}))
colnames(upsetR_cnt) = c('NuRD_bound_ESC', 'DE_Mbd3KO_ESC', 'NuRD_bound_EpiSC', 'DE_Mbd3KO_EpiSC')

#pdf('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/plots/NuRD_bound_promoters-ESCvsEpiSC.pdf')
pdf('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/plots/NuRD_bound_promoters-ESCvsEpiSC_better_labs.pdf')
upset(upsetR_cnt, order.by=c("freq", "degree"), mb.ratio=c(0.7, 0.3),
      queries = list(
                     list( query = intersects, params = list('NuRD_bound_ESC', 'DE_Mbd3KO_ESC'), color = default_colors["lightpinkish"], active = T),
                     list( query = intersects, params = list('NuRD_bound_EpiSC', 'DE_Mbd3KO_EpiSC'), color = default_colors["lightpinkish"], active = T),
                     list( query = intersects, params = list('DE_Mbd3KO_ESC', 'DE_Mbd3KO_EpiSC'), color = default_colors["hotpink"], active = T),
                     list( query = intersects, params = list('NuRD_bound_ESC', 'NuRD_bound_EpiSC', 'DE_Mbd3KO_ESC', 'DE_Mbd3KO_EpiSC'), color = default_colors["hotpink"], active = T),
                     list( query = intersects, params = list('NuRD_bound_ESC', 'NuRD_bound_EpiSC'), color = default_colors["blue"], active = T)
                     ),
      empty.intersections = "on", 
      sets.bar.color = default_colors["grey"],
      mainbar.y.label = 'Gene Sets Intersections', sets.x.label = 'Genes Per Set')
      
grid.newpage()
source('~/source/Rscripts/plotting-functions.R')
a = nrow(esc_targets)
b = nrow(episc_targets)
n12  = table(esc_targets[['Gene']] %in% episc_targets[['Gene']])[['TRUE']]
plot_venn(a, b, n12, c('ESCs targets', 'EpiSCs targets'), c('green', 'purple'))
dev.off()# }}}

comparisons = list.files(pattern = "[^db].tsv", path = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/before_Jan_14/DiffBind/",
                         full.name = T)
comparisons = comparisons[grep('*cluster*', comparisons, invert = T)]
comparisons = comparisons[grep('*binding-per*', comparisons, invert = T)]
comparisons = comparisons[grep('*bound*', comparisons, invert = T)]
comparisons = comparisons[grep('*-2i-*', comparisons, invert = T)]
comparisons = comparisons[grep('*h0-*', comparisons)]
# Keeping only the comparisons of the time points. If design is different this won't work
comparisons = comparisons[grep('!', comparisons, invert = T)]
comparisons = as.data.frame(comparisons)

library(dplyr)
library(tools)
tmp = comparisons %>% .[[1]] %>% as.character() %>% basename() %>% file_path_sans_ext
tmp = gsub(".*_", '', tmp) %>% strsplit("-")
tmp = as.data.frame(do.call(rbind, tmp))
colnames(tmp) = c('protein', 'c1', 'c2')
comparisons = cbind(comparisons, tmp)
args = list()
args$gtf = "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.70.gtf"

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
    annotations$active_enhancers = import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-active_2015-08-10.bed")
    annotations$poised_enhancers = import.bed("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-poised_2015-08-10.bed")
} else {
    command = paste("Rscript ~/source/get-annotation-rscript.R -G", args$gtf, '-p', ncores, collapse = " ")
    stop(paste('Annotation GTFs do not exist. Run `', command, '` first.', collapse = ' '))
    
}# }}}

aggregate_counts = function(fs){ # {{{
    first = TRUE
    for (x in 1:length(fs)) {
        tmp = read.delim(fs[x], head = T, sep = "\t", stringsAsFactors = F)
        rename_map = c(paste("FC", names(fs)[x], sep = "."),
                       paste("FDR", names(fs)[x], sep = "."))
        combine_by = c('seqnames', 'strand', 'start', 'end')
        tmp = tmp %>% select(seqnames, start, end, strand, Fold, FDR, starts_with("Conc")) %>%
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

        plot_data = plot_data %>% dplyr::select(starts_with('Conc_'))

        # Formatting time in row_annotation and plot_data and consequenctly column_annotation and ac# {{{
        #if(ca == "Time") {
        if(FALSE) {
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

            timepoints = data.frame(x = colnames(plot_data)) %>% tidyr::separate(x, into = c('tmp', "timepoint"), sep="_")
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
        column_annotation = data.frame( x = gsub(".*_", '', colnames(plot_data)))
        colnames(column_annotation) = "Time"
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
                    priority = c('2000-gene_start-500', 'active_enhancers',
                                 'poised_enhancers', 'first-exon',  'exons', 'first-intron', 'introns', 'genes')) {

    stopifnot(is(grange, "GRanges"))
    total = length(grange)

    per_region = data.frame(overlapping_feature = rep(NA, total))
    rownames(per_region) = names(grange)

    per_feature = setNames(vector("list", length(priority)), priority)
    bound = data.frame(Gene = character(0), Bound = character(0))
    # '5000-tss-2000' contains all possible TSS for a gene
    for (i in priority) {
        ov = findOverlaps(annotations[[i]], grange)
        # using names and not index as that will change when subsetting granges downstream
        hits = names(grange)[(1:length(grange)) %in% unique(subjectHits(ov))]
        per_region[hits, 'overlapping_feature'] = rep(i, length(hits))

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
    bound_table = data.frame(Gene = values(annotations[['genes']])[['ensembl_gene_id']]) %>% left_join(bound, by="Gene")
    per_region[is.na(per_region$overlapping_feature), 'overlapping_feature'] = 'intergenic'
    row_names = rownames(per_region)
    per_region = per_region %>% mutate(overlapping_feature = ifelse(overlapping_feature %in% c('exons', 'first-exon', 'first-intron', 'intros'), 'gene_body', overlapping_feature))
    rownames(per_region) = row_names

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

proteins = levels(comparisons$protein)
# Heatmap for histone mods across time points# {{{
for (p in proteins){
    print(paste0('Plotting heatmap for ', p))
    fs = comparisons %>% dplyr::select(comparisons) %>% 
            filter(grepl(p, comparisons)) %>%
            unique() %>% .[[1]] %>% as.character()

    names(fs) = mutate(comparisons, label = paste(c1, c2, sep= "-")) %>%
                    filter(grepl(p, comparisons)) %>%
                    dplyr::select(label) %>% .[[1]] %>% as.character()
    fs = mixedsort(fs)

    bdf_cnt = aggregate_counts(fs) %>% mutate(rowname = paste(seqnames, start, end, sep = "_"))
    library('RColorBrewer')
        bdf = tryCatch(format_cdf(bdf_cnt),
                       error=function(e) e,
                       warning=function(w) w)
        if(is(bdf, 'error')) next;

        # Rsgions considered in DB analysis. Not differentially bound ones!
        db_regions = with(bdf_cnt,
                          GRanges(seqnames = seqnames,
                                  IRanges(start, end),
                                  strand = strand))
        names(db_regions) = paste(seqnames(db_regions), start(db_regions), end(db_regions), sep = "_")
        annotated_regions = annotate(db_regions, annotations, output_path, p)

        library(ChIPpeakAnno)
        db_distance_from_genes = annotatePeakInBatch(db_regions, AnnotationData = annotations[['genes']])
        db_distance_from_genes = as.data.frame(db_distance_from_genes)
        dup = db_distance_from_genes %>% filter(duplicated(peak)) %>% select(peak) %>% .[[1]]
        tmp_ra = add_rownames(bdf$ra, var = 'peak')
        # filtering duplicated(peak) because ChIPpeakAnno will report 2 genes if peak has equal distance from them
        # here just selecting the first one reported
        tmp_ra = right_join((add_rownames(annotated_regions$per_region, var='peak') %>% select(-feature_id)), tmp_ra, by ='peak')
        bdf$ra = db_distance_from_genes %>% select(peak, insideFeature, shortestDistance) %>%
                    rename(relative_to_gene = insideFeature,
                           d_relative_to_gene = shortestDistance) %>%
                    mutate(d_relative_to_gene = d_relative_to_gene / 1000) %>%
                    right_join(tmp_ra, by = 'peak') %>% filter(!duplicated(peak))
        rownames(bdf$ra) = bdf$ra[['peak']]
        bdf$ra = bdf$ra[,-1]

        # discretize distance to look good on annotation tracks. Doing it manually cause not sure how to get # {{{
        # what I want with a package
        tmp = c()
        bins = c(0, 0.5, 1, 5, 10, 50, 100, 300)
        for (i in 1:length(bins)) {
            if(i == length(bins)) {
                x = bdf$ra$d_relative_to_gene >= bins[i]
                tmp[x] = paste0(bins[i], 'kb+')
            } else {
                x = bdf$ra$d_relative_to_gene >= bins[i] & bdf$ra$d_relative_to_gene < bins[i+1]
                tmp[x] = paste0(bins[i], '-', bins[i+1], 'kb')
            }
        }
        bdf$ra$d_relative_to_gene = tmp 
        # }}}

        continuous_colors = colorRampPalette(brewer.pal(9,"RdYlGn"))(length(unique(tmp)))
        bdf$ac$d_relative_to_gene = continuous_colors
        names(bdf$ac$d_relative_to_gene) = mixedsort(unique(tmp))
        
        ordering = c('gene_body', '2000-tss-500','intergenic', 'active_enhancers', 'poised_enhancers')
        # Setting levels for factor cause otherwise  pheatmap complains
        bdf$ra$overlapping_feature = factor(bdf$ra$overlapping_feature, levels = ordering)
        bdf$ac$overlapping_feature = rev(colorRampPalette(brewer.pal(8,"Set1"))(length(ordering)))
        names(bdf$ac$overlapping_feature) = ordering
        #not keeping chippeakanno column
        bdf$ra = bdf$ra[-1]
        # not keeping d_relative
        bdf$ra = bdf$ra[-1]
        bdf$ra = bdf$ra[c(length(bdf$ra):1)]
        
        png(file.path(plot_path, paste0('differential_binding-', p, '.png')), height = 700, width = 700)

        results = pheatmap(as.data.frame(bdf$plot_data),
                           clustering_distance_rows = "correlation",
                           clustering_method = "ward.D2", scale = 'row',
                           cluster_cols = FALSE, cluster_rows = TRUE, cutree_rows = 2,
                           annotation_row = bdf$ra, annotation_col = bdf$ca,
                           annotation_legend = T, annotation_colors = bdf$ac,
                           show_rownames = FALSE, 
                           main = paste0('Diff. binding: ', p))
#        cluster = get_cluster(results, output_path = output_path, label = p)
        dev.off()
    }
    #}}}
#}}}

episc = read.delim('/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/deseq/episc_wt-ko.deseq-de.tsv')
gene_names = read.delim(paste0(file_path_sans_ext(args$gtf), '.ensembl2gene_name.tsv'), head = T, sep = "\t")
gene_names = rename(gene_names, Gene = ensembl_gene_id)
episc = episc %>% left_join(gene_names, by = "Gene")
nurd_targets = read.delim('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/NuRD-EpiSCs_unique_sites.binding-per-gene.tsv', head=T)
episc = episc %>% left_join(nurd_targets, by = "Gene")
episc = rename(episc, EpiSC_specific = Bound)
nurd_targets = read.delim('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/NuRD-ESCs_unique_sites.binding-per-gene.tsv', head=T)
episc = episc %>% left_join(nurd_targets, by = "Gene") %>% rename(ESC_specific = Bound)
nurd_targets = read.delim('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/Condition_2i_AND_EpiSC_Proteins_NuRD_Ov1bp.binding-per-gene.tsv', head=T)
episc = episc %>% left_join(nurd_targets, by = "Gene") %>% rename(ESC_EpiSC = Bound)

write.table(episc, "/nfs/research2/bertone/user/mxenoph/hendrich/rna/episc-DEG-with-binding-info.tsv", quote=FALSE, sep="\t", row.names=FALSE)

library(rtracklayer)
gtf = '/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.70.metaseq.genes.filtered.gtf'
deseq = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/deseq/2i_wt-ko.tsv'

deseq = read.table(deseq, header=T, sep="\t")

library(dplyr)
library(dplyrExtras)
library(gtools)

fout = deseq %>% select(gene_id, baseMean:Expression_spl2_2i) %>%
    mutate(groups = ifelse(is.na(padj) | padj >= 0.05, 'unchanged', ifelse(log2FoldChange > 0, 'up', 'down'))) %>%
    mutate(expression = log2(1 + baseMean), quartiles_all = ntile(expression, 4))

tmp = fout %>% filter(is.finite(baseMean) & baseMean != 0) %>% mutate(quartiles_valid = ntile(expression, 4))
fout = left_join(fout, tmp, by = c("gene_id", "baseMean", "expression", "quartiles_all", "padj", "log2FoldChange", "groups")) %>%
    mutate(quartiles_valid = ifelse(is.na(quartiles_valid), 0, quartiles_valid)) %>%
    arrange(log2FoldChange, padj) %>% select(gene_id, groups, expression)

write.table(fout, "/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/results/report/metagene_groups-2i_wt-ko.tsv",
            sep="\t", row.names=F, quote=F)

inducible_path = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/deseq/'
ind_deseq = read.table(paste0(inducible_path, "h48-h0.tsv"), header=T, sep="\t")
ind_fout = ind_deseq %>% select(gene_id, baseMean:Expression_Mbd3_Rescue_48h_3_2i) %>%
    mutate(groups = ifelse(is.na(padj) | padj >= 0.05, 'unchanged', ifelse(log2FoldChange > 0, 'up', 'down'))) %>%
    mutate(expression = log2(1 + baseMean), quartiles_all = ntile(expression, 4))

tmp = ind_fout %>% filter(is.finite(baseMean) & baseMean != 0) %>% mutate(quartiles_valid = ntile(expression, 4))
ind_fout = left_join(ind_fout, tmp, by = c("gene_id", "baseMean", "expression", "quartiles_all", "padj", "log2FoldChange", "groups")) %>%
    mutate(quartiles_valid = ifelse(is.na(quartiles_valid), 0, quartiles_valid)) %>%
    arrange(log2FoldChange, padj) %>% select(gene_id, groups, expression)

write.table(ind_fout, "/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/results/report/metagene_groups-h0-h48.tsv",
            sep="\t", row.names=F, quote=F)

test = '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/DiffBind/TF-steady-state-diff-bind/Nanog.binding-per-gene.tsv'
test = read.delim(test, head=T)
nanog_all = test %>% filter(!is.na(Bound)) %>% select(Gene) %>% .[[1] 

tf_path = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/DiffBind/TF-steady-state-diff-bind"
annotation_path = '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/feature-annotation/'

klf4 = paste0(annotation_path, c('2lox_2i-Klf4_peaks.bound-genes-only.tsv', 'Spl_2i-Klf4_peaks.bound-genes-only.tsv'))


# all peaks in wt and ko and annotated
klf = read.delim(file.path(tf_path, "Klf4.binding-per-region.tsv"), head = T, sep = "\t")
nanog = read.delim(file.path(tf_path, "Nanog.binding-per-region.tsv"), head = T, sep = "\t")

klf = klf %>% filter(grepl("ENS", feature_id))


comparisons = list.files(pattern = "db.tsv", path = tf_path, full.name = T)
# Keeping only the comparisons of the time points. If design is different this won't work
comparisons = comparisons[grep('!', comparisons, invert = T)]

cmp_klf = grep('h0', comparisons[grep('Klf4_', comparisons)], value=T)
x = lapply(cmp_klf, read.delim)

cmp_nanog = grep('h0', comparisons[grep('Nanog_', comparisons)], value=T)

