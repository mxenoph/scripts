#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Get exons, introns, intergenic for given annotation")
parser$add_argument('-b', '--bed', metavar= "file", required='True', type= "character",
                    help= "bed file for PolII ChIP")
parser$add_argument('-g', '--gtf', metavar= "file", required='True', type= "character",
                    default = '/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.70.rtracklayer-genes.gtf',
                    help= "GTF file for genes created with get-annotation-rscript.R")
parser$add_argument('-a', '--annotation', metavar= "file", required='True', type= "character",
                    default = '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/Ens75_genes+transscripts.txt',
                    help= "GTF file for genes provided from Meryem")
parser$add_argument('-s', '--subsets', metavar= "file", type= "character",
                    default = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/deseq/*_de.tsv',
                    help= "TSV files containing subsets")
parser$add_argument('-t', '--targets', metavar= "file", type= "character", default = NULL,
                    default = '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/Condition_7E12_2i_Proteins_M2_meryem_filtered_AND_Chd4_meryem_filtered_Ov1bp.binding-per-gene.tsv',
                    help= "TSV files containing targets for annotating fold changes")
parser$add_argument('p', '--pattern', type = "character", default = NULL, help = "Pattern used to exclude samples from calculating fold change")
parser$add_argument('-u', '--upstream', metavar= "file", required='True',
                    default = 250, help= "Upstream of promoter")
parser$add_argument('-d', '--downstream', metavar= "file", required='True',
                    default = 250, help= "Downstream of promoter")
parser$add_argument('-e', '--downstream_gene_end', metavar= "file", required='True',
                    default = 500, help= "Downstream of gene end")

args = parser$parse_args()

# For testing interactively# {{{
if(FALSE){
    args = list()
    args$annotation = '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/Ens75_genes+transscripts.txt'
    args$bed = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/bowtie/coverage/ddup/"
    args$out = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10"
    args$gtf = "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.70.rtracklayer-genes.gtf"
    args$subsets = "/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/deseq/*_de.tsv"
    args$targets = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/from_meryem/feature-annotation/Condition_7E12_2i_Proteins_M2_meryem_filtered_AND_Chd4_meryem_filtered_Ov1bp.binding-per-gene.tsv"
    args$pattern = "h24-S5P_1"
}# }}}

output_path = file.path(args$out, 'pausing')
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)
# }}}

# helper function to plot in preview.pdf when debugging# {{{
preview = function(f) {
    pdf('~/public_html/preview.pdf')
    eval(f)
    dev.off()
}# }}}

# Return a list with indexes of times in minutes and the respective values in h# {{{
convert_time = function(x) {
    index = grepl('m', x)
    values = x[index]
    # R introduces X in front of the string if starting with number
    values = gsub("^X", '', values)
    values = as.numeric(gsub("m", '', values))
    values = gsub("^", "h", as.character(values / 60))
    return(list(index = index, values = values))
}
# }}}

# calculating travelling ratio based on Remco's/Meryem analysis# {{{
calculate_travelling_ratio_meryem = function(files, annotation = args$annotation){
    library(ChIPpeakAnno)
    library(data.table)
    library(Vennerable)
    library(GenomicRanges)
    library(gridExtra)
    library(gplots)
    require(VennDiagram)
    library(RColorBrewer)
    rename = dplyr::rename

    # Define promoter and gene body positions# {{{
    get_annotation_from_file = function(annotation){
        annotation = read.delim(file = annotation, header = T)
        genes = annotation %>% select(one_of(c("Ensembl.Gene.ID", "Associated.Gene.Name", "Strand", "Description", "Gene.Biotype"))) %>%
                # why are there replicates in the annotation file, basically where does it come from?
                filter(!duplicated(Ensembl.Gene.ID))

        rename_target = c("gene_id","gene_name",
                          "chr","gene_start","gene_end",
                          "transcript_start","transcript_end","strand")
        rename_map = list("Ensembl.Gene.ID", "Associated.Gene.Name", "Chromosome.Name", 
                       "Gene.Start..bp.", "Gene.End..bp.", 
                       "Transcript.Start..bp.", "Transcript.End..bp.", "Strand")
        # for sense genes: promoter_start is upstream 2000 TSS and promoter end is downstream 500 TSS
        # for antisense genes: promoter_start is downstream 500 TSS and promoter end is upstream 2000 TSS
        # promoter_start and promoter_end are not named based on orientation but on what comes first 5->3
        annotation = annotation %>% select(one_of(unlist(rename_map))) %>%
                    dplyr::rename_(.dots = setNames(rename_map, rename_target)) %>%
                    filter(!duplicated(gene_id)) %>%
                    filter(grepl("^\\d{1,2}$", chr, perl = TRUE) | grepl("^[X,Y]{0,1}$", chr, perl = TRUE)) %>%
                    mutate(chr = paste0('chr', chr)) %>%
                    ## bed file format 0-based
                    mutate(promoter_start = ifelse(strand == 1, transcript_start - 2001, transcript_end - 501)) %>%
                    mutate(promoter_end = ifelse(strand == 1, transcript_start + 500, transcript_end + 2000)) %>%
                    mutate(gene_body_start = ifelse(strand == 1, gene_start + 501, gene_start),
                           gene_body_end = ifelse(strand == 1, gene_end, gene_end - 502)) %>%
                    mutate(gene_length = gene_end - gene_start + 1, 
                           long = ifelse(gene_length >= 800, TRUE, FALSE))

        coordinates = rbind((annotation %>% select(gene_id, chr, promoter_start:promoter_end) %>%
                             mutate(id = paste0(gene_id, ':promoter')) %>%
                             rename(start = promoter_start, end = promoter_end)),
                            (annotation %>% select(gene_id, chr, gene_body_start:gene_body_end) %>%
                             mutate(id = paste0(gene_id, ':gene_body')) %>%
                             rename(start = gene_body_start, end = gene_body_end)) %>% filter(start < end)) %>%
                            arrange(chr, gene_id)
        # Filtering out genes for which the gene_body_start is after the body_gene_end
        # i.e. genes that are shorter than 500bp. 
        # IMPORTANT: Meryem's way of filtering out short genes is flawed cause she removes them from the bodycoord
        # but not promcoord which results in genes having promoter coordinates but not gene_body coordinates in coords
        # She potentially removes them further down in her analysis but this may result in logical errors!
        coordinates = coordinates %>% group_by(gene_id) %>%
            mutate(keep = ifelse(start >= end, FALSE, TRUE)) %>% mutate(keep = ifelse(all(keep), TRUE, FALSE)) %>%
            ungroup() %>% filter(keep) %>% select(-keep)
        # strand is set to '+' for all genes because start and end are 5->3 appropriate (seennotation_from_file )
        coordinates_gr = with(coordinates, GRanges(chr, IRanges(start = start, end = end), strand = '+'))
        names(coordinates_gr) = coordinates[['id']]

        return(list(genes = genes, annotation = annotation, coordinates = coordinates))
    }# }}}

    collect_counts = function(files){# {{{
        counts = lapply(files, read.table, header = T)
        names(counts) = basename(gsub("_gene_signal.txt","",files))
        counts = lapply(names(counts), function(x){ tmp = counts[[x]] %>%
                        rename_(.dots = setNames(c('gene.id', 'gene.name', 'prom.max', 'body.max'),
                                                 c('gene_id', 'gene_name',
                                                   paste0(x, '.promoter_max'),
                                                   paste0(x, '.gene_body_max'))))
                        rownames(tmp) = tmp[['gene_id']]
                        tmp
                            })
        # merging list of dataframes into one dataframe (cool use of Reduce, ey?)
        counts = counts %>% Reduce(function(df1, df2) full_join(df1, df2, by = c('gene_id', 'gene_name')), .)
        return(counts)
    }# }}}
    
    travelling_ratio = function(cnt, for_sample){# {{{
        tr = cnt %>% select(gene_id, gene_name, starts_with(for_sample)) %>%
            # making sure the order is promoter, gene_body counts
            select(gene_id, gene_name, ends_with('promoter_max'), ends_with('gene_body_max')) %>%
            # what is the filtering based on, why 20 reads in gene_body? Might have to remove this if normalising based on input
            mutate(total = .[[3]] + .[[4]], 
                   travelling_ratio = round(.[[3]] / (total + 0.001), digits = 2), 
                   travelling_ratio_filtered = ifelse(total < 20, NA, travelling_ratio)) %>%
            rename_(.dots = setNames(c('total', 'travelling_ratio', 'travelling_ratio_filtered'),
                                     c(paste0(for_sample, '.total'),
                                       paste0(for_sample, '.travelling_ratio'),
                                       paste0(for_sample, '.travelling_ratio_filtered'))))
        return(tr)
    }# }}}

    # g() is defined in Rprofile and it assigns the list returned from get_annotation_from_meryem() 
    g(genes, annotation, coordinates) %=% get_annotation_from_file(annotation)
    counts = collect_counts(files)
    travelling_ratio = lapply(basename(gsub("_gene_signal.txt","",files)), function(x) travelling_ratio(counts, x))
    travelling_ratio = travelling_ratio %>% Reduce(function(df1, df2) full_join(df1, df2, by = c('gene_id', 'gene_name')), .)
    return(travelling_ratio)
}
# }}}

calculate_fc = function(travelling_ratio, pattern = "travelling_ratio_filtered"){ #{{{
    per_condition_tr = travelling_ratio %>% select(gene_id, gene_name, ends_with(pattern))
    per_condition_tr =  melt(per_condition_tr) %>%
        tidyr::separate(variable, c('Condition','Replicate'), sep='-') %>%
        group_by(Condition, gene_id) %>% summarise(mean = mean(value), standard_deviation = sd(value)) %>%
        ungroup()

    # id mean is NA after this, it means that at least one replicate was NA
    per_condition_tr = lapply(split(per_condition_tr, per_condition_tr$Condition), function(x){ 
        label = unique(x$Condition)
        x %>% dplyr::rename_(.dots = setNames(c('mean', 'standard_deviation'), 
                                              paste(label, c('mean', 'standard_deviation'), sep='-'))) %>%
        select(-Condition)}) %>%
        Reduce(function(df1, df2) full_join(df1, df2, by = 'gene_id'), .)

    # relies on the fact that naming is consistent: sample_condition-factor_rep
    conditions = strsplit(gsub('-mean', '', colnames(per_condition_tr %>% select(matches('mean')))), '_')
    tmp = sapply(conditions, function(x) x[2])
    map = convert_time(tmp)
    tmp[map$index] = map$values
    conditions = paste(sapply(conditions, function(x) x[1]), tmp, sep='_')
    # if any changes happened to the names then I need to change them on the dataframe too and need to do that before using mixedsort
    tmp = colnames(per_condition_tr %>% select(matches('mean')))
    # because if the labels start with numbers then parse doesn't work
    tmp = paste0("`", tmp, "`")
    per_condition_tr = per_condition_tr %>% dplyr::rename_(.dots = setNames(tmp, paste0(conditions, '-mean')))
    per_condition_tr = per_condition_tr[, c('gene_id', paste0(mixedsort(conditions), '-mean'))]

    combination = combn(mixedsort(conditions), 2)
    fold_change = apply(combination, 2, function(x){
                            nominator = per_condition_tr %>% select(matches(paste0(x[2], '-mean')))
                            denominator = per_condition_tr %>% select(matches(paste0(x[1], '-mean')))
                            fc = (nominator / denominator)
                            fc = data.frame(gene_id = per_condition_tr[['gene_id']], fc = fc)
                            colnames(fc) = c('gene_id', paste(x[1], x[2], sep = '-'))
                            return(fc)
        }) %>%
        Reduce(function(df1, df2) full_join(df1, df2, by = 'gene_id'), .)

    fold_change[, sapply(fold_change, is.numeric)] = log2(fold_change[, sapply(fold_change, is.numeric)])
    return(list(per_condition_tr = per_condition_tr, fc = fold_change, conditions = conditions))
} #}}}

format_melted = function(melted_df){# {{{
    contains_times = any(gregexpr("h[[:digit:]]+", melted_df$variable, ignore.case = T) != -1) | any(gregexpr("[[:digit:]]+h", melted_df$variable, ignore.case = T) != -1) | any(gregexpr("[[:digit:]]+m", melted_df$variable, ignore.case = T) != -1)
    if(contains_times){
        melted_df = melted_df %>% separate(variable, c("Sample", "Replicate"), sep = "-") %>% 
            separate(Replicate, c("Replicate", "tmp"), sep = "\\.") %>% 
            mutate(Time = gsub('_', '', str_extract(Sample, "_.*")),
                   Sample = gsub('_', '', str_extract(Sample, ".*_")),
                   Mark = gsub('_', '', str_extract(Replicate, "[^_]*_")),
                   Replicate = gsub('_', '', str_extract(Replicate, "_\\d_"))) %>%
            select(-tmp)
        
        tmp = convert_time(melted_df$Time)
        melted_df[tmp$index, 'Time'] = tmp$values
        melted_df$Time = factor(melted_df$Time, levels = mixedsort(unique(melted_df$Time)))

        tmp = sapply(melted_df %>% select(-gene_id, -value), function(x) !any(is.na(x)))[c('Sample', 'Time', 'Mark', 'Replicate')]
        if(all(tmp)){
            melted_df = melted_df %>% mutate(variable = paste0(Sample,'_', Time, '-', Mark, '_', Replicate)) %>%
                select(-Sample, -Mark) %>% arrange(Time, variable)
        } else{
            # If not all(tmp) != TRUE working on an fc dataframe where Mark and Replicate columns are missing
            melted_df = melted_df %>% mutate(variable = paste0(Sample,'_', Time)) %>% 
                select(-Sample, -Mark, -Replicate) %>% arrange(Time, variable)
        }
        melted_df$variable = factor(melted_df$variable, levels = mixedsort(unique(melted_df$variable)))
    } else{
        melted_df$variable = factor(melted_df$variable, levels = mixedsort(unique(melted_df$variable)))
        melted_df = melted_df %>% arrange(variable)
    }
    return(melted_df)
}# }}}

annotate_df = function(df, metadata = args$targets){# {{{
    metadata = read.delim(metadata)
    # Ensure that genes in metadata object are as many as in df => annotation is the same
    if(nrow(df) != nrow(metadata)) warning(paste0("Number of genes in fold change table does not match number of genes in metadata provided (genes in fold change table:",
                                                nrow(df), ", metadata:", nrow(metadata), ", missing:", nrow(df) - nrow(metadata), ").\n"))
    x = grep('Gene', colnames(metadata), ignore.case = TRUE)
    if(length(x) >= 1){
        # there's a column for gene ids
        keep = sapply(x, function(y) grep('ENS', metadata[1, y]))
        to_rename = colnames(metadata)[keep]
        metadata = metadata %>% dplyr::rename_(.dots = setNames(to_rename, 'gene_id'))

        # join whatever the columns are, don't check or filter based on Bound column
        # This means I need to be sure of the annotation used for metadata matches the one used in this script
        # metadata = metadata %>% filter(!is.na(Bound))
        annotated = df %>% full_join(metadata, by = 'gene_id') %>% 
            mutate(annotation_files = ifelse((gene_id %in% df[['gene_id']]) & (gene_id %in% metadata[['gene_id']]),
                                             'both',
                                             ifelse((gene_id %in% df[['gene_id']]), 'args_annotation', 'target_file')))
        return(annotated)
    } else {
        stop("Not annotating genes with regard to metadata provided. No gene column found")
    }
}# }}}

# Plotting density# {{{
plot_density = function(pausing, subsets = NULL, subset_name = '', label = '', metric = 'Travelling Ratio', combination = NULL){
    get_valid_n = function(pausing){# {{{
        return_list = list()
        if('Group' %in% colnames(pausing)){
            n_valid_subset = pausing %>% group_by(Group, variable) %>% dplyr::summarise(number = table(is.na(value))['FALSE'])
            return_list$n_valid_subset = n_valid_subset
        }

        n_valid = pausing %>% group_by(variable) %>%
            dplyr::summarise(number = table(is.na(value))['FALSE'])
        # Remove NA values
        pausing = pausing %>% filter(!is.na(value))
        
        return_list$pausing = pausing
        return_list$n_valid = n_valid
        return(return_list)
    }# }}}

    test_significance = function(df, combination){# {{{
        stats = apply(combination, 2, function(x){
                          distributions = list((df %>% filter(str_detect(variable,
                                                                         paste0("^", x[1], "\\-{0,1}[^\\.]*$"))) %>% .[['value']]),
                                               (df %>% filter(str_detect(variable,
                                                                         paste0("^", x[2], "\\-{0,1}[^\\.]*$"))) %>% .[['value']]))
                          names(distributions) = x
                          ks = ks.test(distributions[[1]], distributions[[2]])
                          ks = data.frame(comparison = paste0(x[1], '-VS-', x[2]),
                                          pvalue = ks$p.value, alternative = ks$alternative, method = ks$method)
                   }) %>% Reduce(function(df1, df2) bind_rows(df1, df2), .)
        return(stats)
    }# }}}

    # ggplot extrapolated when computing the ecdf (stat_ecdf(pad = F) does not work).# {{{
    # Temporary fix from https://github.com/hadley/ggplot2/issues/1467
    stat_myecdf = function(mapping = NULL, data = NULL, geom = "step",
                            position = "identity", n = NULL, na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE, direction="vh", ...) {
        layer(data = data,
              mapping = mapping,
              stat = StatMyecdf,
              geom = geom,
              position = position,
              show.legend = show.legend,
              inherit.aes = inherit.aes,
              params = list(n = n,
                            na.rm = na.rm,
                            direction=direction, ...))
    }

    StatMyecdf = ggproto("StatMyecdf", Stat,
                         compute_group = function(data, scales, n = NULL) {
                             # If n is NULL, use raw values; otherwise interpolate
                             if (is.null(n)) {
                                 # Dont understand why but this version needs to sort the values
                                 xvals = sort(unique(data$x))
                             } else {
                                 xvals = seq(min(data$x), max(data$x), length.out = n)
                             }
                             y = ecdf(data$x)(xvals)
                             x1 = max(xvals)
                             y0 = 0
                             data.frame(x = c(xvals, x1), y = c(y0, y))
                         },
                         default_aes = aes(y = ..y..),
                         required_aes = c("x"))# }}}

    get_ks_dist = function(first_dist, second_dist) {# {{{
        first_cdf = ecdf(first_dist)
        second_cdf = ecdf(second_dist)
        # find min and max statistics to draw line between points of greatest distance
        min_max_stats = seq(min(first_dist, second_dist), max(first_dist, second_dist),
                            length.out = min(length(first_dist), length(second_dist)))
        x0 = min_max_stats[which( abs(first_cdf(min_max_stats) - second_cdf(min_max_stats)) == max(abs(first_cdf(min_max_stats) - second_cdf(min_max_stats))) )]
        y0 = first_cdf(x0)
        y1 = second_cdf(x0)
        return(list('x0' = x0, 'y0' = y0, 'y1' = y1))
    }# }}}

    return_list = list()
    if(! 'Gene' %in% colnames(pausing)){
        pausing = add_rownames(as.data.frame(pausing), var = "Gene")
    }
    g(pausing, n_valid) %=% get_valid_n(pausing)
    
    # sample column to be used in violinplot so I can tweak the geom_text positions
    annotations = pausing %>% 
        mutate(Sample = plyr::mapvalues(variable, from = levels(variable), to = 1:length(levels(variable)))) %>% 
        select(-one_of('Gene', 'gene_id', 'gene_name', 'value')) %>% unique() %>%
        inner_join(n_valid, by = 'variable')
    
    # Plot basics# {{{
    p = ggplot(data = pausing)
    p = p + theme_bw() + theme(legend.position= "right", aspect.ratio = 1)
    p = p + ggtitle(label) + xlab(metric) + labs(colour = 'Sample')
    
    aesthetics_linetype = NULL
    if('Replicate' %in% colnames(pausing)){
        aesthetics_linetype = 'Replicate'
        tmp = annotations$Replicate
    } else {
        tmp = annotations$variable
    }

    if('Time' %in% colnames(pausing)){# {{{
        divergent_colors_seq = setNames(colorRampPalette(brewer.pal(11,"PRGn")[c(1:4,8:11)])(length(levels(pausing$Time))),
                                  levels(pausing$Time))
        names(divergent_colors_seq) = levels(pausing$Time)
        fill = 'Time'
        cols = divergent_colors_seq
        
#        pviolin = p + geom_violin(aes_string(x = 'variable', y = 'value', fill = 'Time'))
#        pviolin = pviolin + scale_fill_manual(values = divergent_colors_seq) + labs(fill = 'Time')
        
#        pcumulative = p + scale_colour_manual(values = divergent_colors_seq) + labs(colour = 'Time') + ylab('CDF')
#        pcumulative = pcumulative + scale_shape_discrete(solid = F)
        #pcumulative = p + stat_ecdf(aes_string('value', color = 'Time', shape = aesthetics_linetype), geom = "point", pad = F)
#        pcumulative_final = pcumulative + geom_point(aes_string('value', color = 'Time', linetype = aesthetics_linetype), stat ="myecdf", size = 2)

        # density = counts / sum(counts * bar width)
#        p = p + geom_line(stat = "density", aes_string(x = 'value', color = 'Time', linetype = aesthetics_linetype))
#        p = p + scale_colour_manual(values = divergent_colors_seq) + labs(colour = 'Time')
#        tmp = paste0(annotations[['Time']], '_', tmp)
    } else{
        pviolin = p + geom_violin(aes_string(x = 'variable', y ='value', fill = 'variable')) + labs(fill = 'Sample')
        
        pcumulative = p + scale_shape_discrete(solid = T) + ylab('CDF')
        #pcumulative = pcumulative + stat_ecdf(aes_string('value', color = 'variable', shape = aesthetics_linetype), geom = "point")
        
        p = p + geom_line(stat = "density", aes_string(x = "value", color = 'variable', linetype = aesthetics_linetype))

        if(any(duplicated(annotations$Replicate))) tmp = as.character(annotations$variable)

        if(length(levels(pausing$variable)) > 9){
            default_colors= setNames(colorRampPalette(brewer.pal(8,"Set2"))(length(levels(pausing$variable))),
                                      levels(pausing$variable))
            names(default_colors) = levels(pausing$variable)
            cols = default_colors

#            pviolin = pviolin + scale_fill_manual(values = default_colors)
#            pcumulative = pcumulative + scale_colour_manual(values = default_colors)
#            p = p + scale_colour_manual(values = default_colors)
        } else{
            cols = "Set2"
#            pviolin = pviolin + scale_fill_brewer(palette = "Set2")
#            pcumulative = pcumulative + scale_colour_brewer(palette = "Set2")
#            p = p + scale_colour_brewer(palette = "Set2")
        }
        pcumulative_final = pcumulative + geom_point(aes_string('value', color = 'variable', shape = aesthetics_linetype), stat ="myecdf")
    }# }}}
    pviolin = p + geom_violin(aes_string(x = 'variable', y = 'value', fill = fill))
    pviolin = pviolin + scale_fill_manual(values = cols) + labs(fill = fill)
        
    pcumulative = p + scale_colour_manual(values = divergent_colors_seq) + labs(colour = fill) + ylab('CDF')
    pcumulative = pcumulative + scale_shape_discrete(solid = F)
    #pcumulative = p + stat_ecdf(aes_string('value', color = 'Time', shape = aesthetics_linetype), geom = "point", pad = F)
    pcumulative_final = pcumulative + geom_point(aes_string('value', color = 'Time', linetype = aesthetics_linetype), stat ="myecdf", size = 2)

        # density = counts / sum(counts * bar width)
        p = p + geom_line(stat = "density", aes_string(x = 'value', color = 'Time', linetype = aesthetics_linetype))
        p = p + scale_colour_manual(values = divergent_colors_seq) + labs(colour = 'Time')

    # adding annotations
    pviolin_final = pviolin + annotate("text",
                                       x = annotations$Sample %>% as.numeric() + 0.35,
                                       y = 0,
                                       # left justified
                                       hjust = 0,
                                       label = paste0("n = ", annotations$number))
    pviolin_final = pviolin_final + coord_flip() + ylab(metric) + xlab('')
    
    # Getting the ggplot range for y axis
    from_y_limit = ggplot_build(p)$panel$ranges[[1]]$y.range[2]
    to_y_limit = from_y_limit - (nrow(annotations) * 0.15)

    # because if the labels start with numbers then parse doesn't work
    tmp = paste0("`", tmp, "`")
    p_final = p + annotate("text",
                     x = 0,
                     y = seq(from = from_y_limit, to = to_y_limit, length.out = nrow(annotations)),
                     # left justified
                     hjust = 0, parse = TRUE,
                     label = paste("n[", eval(tmp), "] == ", annotations$number))

    plot(p_final)
    plot(pviolin_final)
    plot(pcumulative_final)
    # }}}

    if(!is.null(combination)){# {{{
        apply(combination, 2, function(x){
                  pairwise_comp = pausing %>% filter(str_detect(variable, paste0("^", x[1], "\\-{0,1}[^\\.]*$")) | str_detect(variable,  paste0("^", x[2], "\\-{0,1}[^\\.]*$"))) %>% droplevels()
                  
                  if(length(levels(pairwise_comp$variable)) == 2){
                      # finding the greatest distance in the ECDFs
                      tmp1 = (pairwise_comp %>% filter(variable == levels(variable)[1]) %>% .[['value']])
                      tmp2 = (pairwise_comp %>% filter(variable == levels(variable)[2]) %>% .[['value']])
                      g(x0, y0, y1) %=% get_ks_dist(tmp1, tmp2)
                      df = data.frame(x0 = x0[1], y0 = y0[1], y1 = y1[1])
                      pcumulative = pcumulative %+% pairwise_comp
                      pcumulative = pcumulative + geom_line(aes_string('value', color = 'Time', linetype = aesthetics_linetype),
                                                            stat ="myecdf", size = 2)

                      pcumulative = pcumulative + geom_segment(data = df, aes(x = x0, y = y0, xend = x0, yend = y1),
                                                               linetype = "dotted", color = "#A94777", size = 1)
                      pcumulative = pcumulative + geom_point(data = df, aes(x = x0 , y = y0), color = "#A94777", size = 3)
                      pcumulative = pcumulative + geom_point(data = df, aes(x = x0 , y = y1), color = "#A94777", size = 3)

                      # Testing if the distributions are the same with Kolmogorov smirnoff 
                      stats = test_significance(pairwise_comp, as.matrix(x))
                      pcumulative = pcumulative + annotate("text",
                                                           x = 0,
                                                           y = 0.95,
                                                           # left justified
                                                           hjust = 0,
                                                           label = paste0(stats$method, " (", stats$alternative, ")\n",
                                                                          "p-value = ", signif(stats$pvalue, 3)))
                      return_list$stats = stats
                  } else {
                      print('in if')
                      print(aesthetics_linetype)
                      preview(plot(pcumulative_final))
                      pcumulative = pcumulative_final %+% pairwise_comp
                  }
                  plot(pcumulative)})
    }# }}}

    # Subsets # {{{
    if (!is.null(subsets)){
        print('in subsets')
        pausing_subset = pausing %>% inner_join(subsets, by = "gene_id") %>% droplevels()
        g(n_valid_subset, pausing_subset, n_valid) %=% get_valid_n(pausing_subset)

        annotations = pausing %>% 
            mutate(Sample = plyr::mapvalues(variable, from = levels(variable), to = 1:length(levels(variable)))) %>% 
            select(-one_of('Gene', 'gene_id', 'gene_name', 'value')) %>% unique() %>%
            inner_join(n_valid_subset, by = 'variable')

        q = p %+% pausing_subset + facet_grid(.~ Group) + ggtitle(paste(label, subset_name, sep=''))
        
        n_groups = length(levels(factor(pausing_subset$Group)))
        # Getting the ggplot range for y axis
        from_y_limit = ggplot_build(q)$panel$ranges[[1]]$y.range[2]
        to_y_limit = from_y_limit - ((nrow(annotations)/n_groups) * 0.5)
        y_range4labels = seq(from = from_y_limit, to = to_y_limit, length.out = (nrow(annotations)/n_groups))
        y_range4labels = rep(y_range4labels, n_groups)
        q = q + geom_text(data = annotations,
                         aes(x = 0, y = y_range4labels),
                         # left justified
                         hjust = 0, parse = TRUE,
                         label = paste("n[", rep(tmp, n_groups), "] == ", annotations$number))

        # not adding annotation to this plot as it gets too crowded
        qviolin = pviolin %+% pausing_subset + facet_grid(.~ Group)
        qviolin = qviolin + ggtitle(paste(label, subset_name, sep='')) + coord_flip() + ylab(metric) + xlab('')

        plot(q)
        plot(qviolin)

        if(!is.null(combination)){
            preview(apply(combination, 2, function(x){
                      pairwise_comp = pausing_subset %>% filter(grepl(glob2rx(paste0("*", x[1], "-*")), variable) | grepl(glob2rx(paste0("*", x[2], "-*")), variable))
                      pcumulative = pcumulative %+% pairwise_comp
            plot(pcumulative)
            return(x)}))
        }
    }# }}}

}# }}}

#PRR = promoter release ratios from Fei Xavier Chen 2015
calculate_PRR = function(rpms, tss, peaks){

}

calculate_pausing_index = function(){
    promoter_avg / gene_body_avg
}

main(){
    library(dplyr)
    library(tidyr)
    # for extracting pattern from string with str_extract
    library(stringr)
    library(rtracklayer)

    files = list.files(pattern=glob2rx(paste("*_gene_signal*",".txt",sep="")), args$bed, full.name = T)
    travelling_ratio = calculate_travelling_ratio_meryem(files)
    write.table(travelling_ratio,
                file = file.path(output_path, 'travelling-ratios-unfiltered.tsv'), quote = F, row.names = F, sep ="\t")
    write.table((fc_annotated %>% select(-annotation_files)),
                file = file.path(output_path, 'travelling-ratios-fc-annotated.tsv'), quote = F, row.names = F, sep ="\t")

    g(per_condition_tr, fc, conditions) %=% calculate_fc(travelling_ratio)
    write.table(per_condition_tr,
                file = file.path(output_path, 'travelling-ratios-replicate-summary.tsv'), quote = F, row.names = F, sep ="\t")
    write.table(fc,
                file = file.path(output_path, 'travelling-ratios-fc.tsv'), quote = F, row.names = F, sep ="\t")

    dlab = ''
    # removing samples from calculating TR per condition and consequently TR fold change# {{{
    # all subsequent plots and calculations are done based on the fc computed without those samples
    if(!is.null(args$pattern)){
        g(per_condition_tr, fc, conditions) %=% calculate_fc(travelling_ratio %>% select(-matches(args$pattern)))
        dlab =  paste0('-excluding-', gsub("\\*", "ANY", args$pattern))
        write.table(per_condition_tr,
                    file = file.path(output_path, paste0('travelling-ratios-replicate-summary',
                                                         dlab, '.tsv')),
                    quote = F, row.names = F, sep ="\t")
        write.table(fc,
                    file = file.path(output_path, paste0('travelling-ratios-fc',
                                                         dlab, '.tsv')),
                                     quote = F, row.names = F, sep ="\t")
    }# }}}

    if(!is.null(args$targets)){
        fc_annotated = tryCatch(annotate_df(fc),
                       error=function(e) e)
        if(is(fc_annotated, 'error')) {
            warning('Could not annotate the data frame (no file written).')
        } else {
            write.table((fc_annotated %>% select(-annotation_files)),
                        file = file.path(output_path, paste0('travelling-ratios-fc-annotated', dlab, '.tsv')),
                        quote = F, row.names = F, sep ="\t")
        }
    }

    # Melt travelling ratios and order according to times if samples from timecourse
    TR = melt(travelling_ratio %>% select(gene_id, gene_name, ends_with("travelling_ratio_filtered")))
    TR = format_melted(TR)

    pdf(file.path(plot_path, 'density-plots.pdf'))
    plot_density(TR, combination = combn(mixedsort(conditions),2), label = 'Full dataset')
    dev.off()
    
    # per_condition_tr is of class "tbl_df" and melt doesn't play nice with that
    tmp = melt(as.data.frame(per_condition_tr))
    tmp = tmp %>% separate(variable, c("variable", 'x'), sep = '-') %>% select(-x)
    tmp = format_melted(tmp)
    
    pdf(file.path(plot_path, paste0('density-plots-mean', dlab, '.pdf')))
    plot_density(tmp, combination = combn(mixedsort(conditions),2), label = paste0('Mean per condition', dlab))
    dev.off()

    subsets = travelling_ratio[1:100, 1, drop=F]
    subsets = subsets %>% mutate(Group = c(rep('A', 50), rep('B', 50)))
    preview(plot_density(TR, subsets = subsets, subset_name = 'test'))

}

calculating_travelling_ratio_maria = function(){# {{{
    args$bed = '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/bowtie/coverage/ddup/3KO_30m-S5P_1_ddup.bedgraph'
    coverage = import.bed(args$bed)
    genes = import.gff3(args$gtf)

# Ensuring seqlengths are there for using trim() afterwards
    if(any(is.na(seqlengths(genes)))){
        chromosomes = file.path(dirname(args$gtf), paste0(basename(dirname(args$gtf)), '.genome'))
        chromosomes = read.delim(chromosomes, head = FALSE, sep = '\t')
        chromosomes = setNames(chromosomes[[2]], chromosomes[[1]])
        seqlengths(genes) = chromosomes[names(seqlengths(genes))]
    }

    extended_genes = trim(resize(genes, fix = 'start', width = width(genes) + 500))
    gene_start = trim(promoters(genes,
                              downstream = args$downstream,
                              upstream = args$upstream))

    gene_body_antisense = extended_genes[strand(extended_genes) == "-"]

    gene_body_sense = extended_genes[strand(extended_genes) == "+"]
    gene_body_sense = GRanges(seqnames(gene_body_sense),
                  IRanges(end(gene_start[strand(gene_start) == "+"]) + 1, end(gene_body_sense)),
                  strand(gene_body_sense))
    values(gene_body_sense) = values(gene_body_sense)


    gene_body = GRanges(seqnames(genes), IRanges(end(gene_start) + 1, end(flank(genes, start=FALSE, width=500))), strand(genes))
}# }}}
