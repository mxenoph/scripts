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

# reset_par# {{{
reset_par = function(){
    # Reset par to defaults in the parent environment
    eval(quote(par(mfrow = c(1,1))), parent.frame())
    eval(quote(par(oma = c(0,0,0,0))), parent.frame())
    eval(quote(par(mar = c(5, 4, 4, 2) + 0.1)), parent.frame())}
# }}}

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
    per_condition_tr = travelling_ratio %>% select(gene_id, gene_name, ends_with(pattern)) %>%
                        gather(variable, value, -gene_id, -gene_name) %>%
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
        return(list(annotated = annotated, targets = metadata))
    } else {
        stop("Not annotating genes with regard to metadata provided. No gene column found")
    }
}# }}}

test_significance_ks = function(df, combination){# {{{
    # unlike the t-statistic, the value of the D statistic (and hence the P value) 
    # is not affected by scale changes like using log. The KS-test is a robust test that 
    # cares only about the relative distribution of the data.
    # Unequal dataset size is not a problem for the KS-test.
    # distribution of the test statistic is based on the assumption that the distributions are continuous 
    # => ties are impossible. When there are ties the distribution is affected in such a way that depends on the 
    # pattern of ties and approximations of the p-values are reported instead. p-value is HEAVILY affected when 
    # large number of ties
    # If so, then that will affect the results, the assumption is that the 2 samples are independent of each 
    # other and if you have some patients values in both data sets then that makes them not independent.
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

test_significance_t = function(df, combination){# {{{
    # unless the deviation from normality is really obvious uou shouldn't worry
    # about using the t-test
   stats =  apply(combination, 2, function(x){
                      distributions = df %>% filter(str_detect(variable, paste0("^", x[1], "\\-{0,1}[^\\.]*$")) | str_detect(variable,  paste0("^", x[2], "\\-{0,1}[^\\.]*$"))) %>% droplevels()

                      if(length(levels(distributions$variable)) == 2){
                          effect_size = seq(0.1, 1, 0.2)
                          power_estimation = pwr.t.test(n = nrow(distributions), d = effect_size,
                                                        sig.level = c(0.05), type = "paired")
                          power_estimation = data.frame(effect = effect_size,
                                                       power = power_estimation$power,
                                                       significance_level = power_estimation$sig.level,
                                                       method = power_estimation$method)
                          write.table(power_estimation, file.path(output_path,
                                                                  paste0(x[1], '-VS-', x[2],'.t-test-detection-power.tsv')),
                                      quote = FALSE, row.names = FALSE, sep = "\t")
                          ttest = t.test((distributions %>% filter(variable == levels(variable)[1]) %>% .[['value']]),
                                         (distributions %>% filter(variable == levels(variable)[2]) %>% .[['value']]),
                                         conf.level = 0.95)
                                         #conf.level = 0.95, alternative = "two.sided")
                          ttest = data.frame(comparison = paste0(x[1], '-VS-', x[2]),
                                             pvalue = ttest$p.value, df = ttest$parameter, null_value = ttest$null.value,
                                             t_statistic = ttest$statistic,
                                             alternative = ttest$alternative, method = ttest$method)
                      } else{
                          ttest = NA
                      }
               })
    if(all(is.na(stats))) return(NA)

    stats = stats[!is.na(stats)] %>% Reduce(function(df1, df2) bind_rows(df1, df2), .) %>% 
            mutate(FDR = p.adjust(pvalue, method = "BH", n = nrow(.)))
    return(stats)
}# }}}

# Plotting density# {{{
plot_density = function(pausing, subsets = NULL, subset_name = '', label = '', metric = 'Travelling Ratio', combination = NULL){
    library(ggplot2)
    library(RColorBrewer)

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

    return_list = list()
    if(! 'gene_id' %in% colnames(pausing)){
        pausing = add_rownames(as.data.frame(pausing), var = "gene_id")
    }
    g(pausing, n_valid) %=% get_valid_n(pausing)
    
    # sample column to be used in violinplot so I can tweak the geom_text positions
    annotations = pausing %>% 
        mutate(Sample = plyr::mapvalues(variable, from = levels(variable), to = 1:length(levels(variable)))) %>% 
        select(-one_of('gene_id', 'gene_name', 'value')) %>% unique() %>%
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
        # set what function to be called for colours and with what parameters
        fcol = 'scale_colour_manual'
        fparam = 'values'
        tmp = paste0(annotations[['Time']], '_', tmp)
    } else{
        fill = 'variable'
        if(any(duplicated(annotations$Replicate))) tmp = as.character(annotations$variable)

        if(length(levels(pausing$variable)) > 9){
            default_colors= setNames(colorRampPalette(brewer.pal(8,"Set2"))(length(levels(pausing$variable))),
                                      levels(pausing$variable))
            names(default_colors) = levels(pausing$variable)
            cols = default_colors
            fcol = 'scale_colour_manual'
            fparam = 'values'

        } else{
            cols = "Set2"
            fcol = 'scale_colour_brewer'
            fparam = 'palette'
        }
    }# }}}

    # Violin plots# {{{
    #pviolin = p + geom_violin(aes_string(x = 'variable', y = 'value', colour = fill))
    pviolin = p + geom_boxplot(aes_string(x = 'variable', y = 'value', colour = fill), notch = FALSE, varwidth = TRUE)
    #pviolin = pviolin + eval(parse(text = paste(gsub('colour', 'fill', fcol), '(', fparam, '=cols', ')'))) + labs(colour = fill)
    pviolin = pviolin + eval(parse(text = paste(fcol, '(', fparam, '=cols', ')'))) + labs(colour = fill)
    
    # adding annotations
    pviolin = pviolin + annotate("text",
                                       x = annotations$Sample %>% as.numeric() + 0.35,
                                       y = 0,
                                       # left justified
                                       hjust = 0,
                                       label = paste0("n = ", annotations$number))
    pviolin = pviolin + coord_flip() + ylab(metric) + xlab('')
    # }}}
    
    # ECDF plots# {{{
    #pcumulative = p + stat_ecdf(aes_string('value', color = fill, shape = aesthetics_linetype), geom = "point", pad = F)
    pcumulative = p + geom_point(aes_string('value', color = fill, shape = aesthetics_linetype), stat ="myecdf", size = 2)
    pcumulative = pcumulative + scale_shape_discrete(solid = F)
    pcumulative = pcumulative + eval(parse(text = paste(fcol, '(', fparam, '=cols', ')'))) + labs(colour = fill) + ylab('CDF')
    # }}}

    # Density plots# {{{
    # density = counts / sum(counts * bar width
    pdens = p + geom_line(stat = "density", aes_string(x = 'value', color = fill, linetype = aesthetics_linetype))
    pdens = pdens + eval(parse(text = paste(fcol, '(', fparam, '=cols', ')'))) + labs(colour = fill)

    # Getting the ggplot range for y axis
    from_y_limit = ggplot_build(pdens)$panel$ranges[[1]]$y.range[2]
    to_y_limit = from_y_limit - (nrow(annotations) * 0.15)

    # because if the labels start with numbers then parse doesn't work
    tmp = paste0("`", tmp, "`")
    pdens = pdens + annotate("text",
                                  x = 0,
                                  y = seq(from = from_y_limit, to = to_y_limit, length.out = nrow(annotations)),
                                  # left justified
                                  hjust = 0, parse = TRUE,
                                  label = paste("n[", eval(tmp), "] == ", annotations$number))
    # }}}

    plot(pdens)
    plot(pviolin)
    plot(pcumulative)
    # }}}

    if(!is.null(combination)){# {{{
        # Testing if the difference of the distributions is zero 
        stats = test_significance_t(pausing, combination)
        if(is(stats, 'data.frame')) {
            filename = gsub(' ', '_', label)
            if(filename != '') filename = paste0(filename, '.')

            write.table(stats, file.path(output_path, paste0(filename, 'paired-t-test.tsv')),
                        quote = FALSE, row.names = FALSE, sep = "\t")
        }
        apply(combination, 2, function(x){
                  pairwise_comp = pausing %>% filter(str_detect(variable, paste0("^", x[1], "\\-{0,1}[^\\.]*$")) | str_detect(variable,  paste0("^", x[2], "\\-{0,1}[^\\.]*$"))) %>% droplevels()
                  
                  if(length(levels(pairwise_comp$variable)) == 2){
                      pcumulative = p %+% pairwise_comp
                      pcumulative = pcumulative + geom_line(aes_string('value', color = fill, linetype = aesthetics_linetype),
                                                            stat ="myecdf", size = 2)
                      pcumulative = pcumulative + eval(parse(text = paste(fcol, '(', fparam, '=cols', ')'))) 
                      pcumulative = pcumulative + labs(colour = fill) + ylab('CDF')

                      tmp = stats %>% filter(comparison == paste0(x[1], '-VS-', x[2]))
                      pcumulative = pcumulative + annotate("text",
                                                           x = 0,
                                                           y = 0.95,
                                                           # left justified
                                                           hjust = 0,
                                                           label = paste0(tmp$method, " (", tmp$alternative, ")\n",
                                                                          "FDR = ", round(tmp$FDR, 3), "\n"))#,
                                                                         # "p-value = ", round(tmp$pvalue, 3), "\n"))
                      return_list$stats = stats
                  } else {
                      pcumulative = pcumulative %+% pairwise_comp
                  }
                  plot(pcumulative)})
    }# }}}

    # Subsets # {{{
    if (!is.null(subsets)){
        pausing_subset = pausing %>% inner_join(subsets, by = "gene_id") %>% droplevels()
        g(n_valid_subset, pausing_subset, n_valid) %=% get_valid_n(pausing_subset)

        annotations_subset = pausing_subset %>% 
            mutate(Sample = plyr::mapvalues(variable, from = levels(variable), to = 1:length(levels(variable)))) %>% 
            select(-one_of('gene_id', 'gene_name', 'value')) %>% unique() %>%
            inner_join(n_valid_subset, by = c('variable', 'Group')) %>%
            mutate(tmp = paste0("n[", variable, "] == ", number), 
                   x = as.numeric(Sample) + 0.35)
        
        pausing_subset_ann = pausing_subset %>% 
            left_join(annotations_subset, by = c(unique(c(fill, 'variable')), 'Group'))

        a = pausing_subset_ann %>% group_by(Group) %>%
            do(a = ggplot(data = .) +
               geom_boxplot(aes_string(x = 'variable', y = 'value', colour = fill), notch = FALSE, varwidth = TRUE) +
               facet_grid(~ Group) +
               eval(parse(text = paste(fcol, '(', fparam, '=cols', ')'))) + labs(colour = fill) +
               annotate("text", y = 0, x = unique(.$x),
                        label = paste0("n = ",
                                       data.frame(a = .$variable, b = .$number) %>% unique() %>% .[['b']]),
                        hjust = 0) + 
               ggtitle(paste(label, subset_name, sep='')) + coord_flip() + ylab(metric) + xlab('') +
               theme_bw() + theme(legend.position= "right", aspect.ratio = 1)
                )
        b = pausing_subset_ann %>% group_by(Group) %>%
            do(pviolin %+% . +
               geom_boxplot(aes_string(x = 'variable', y = 'value', colour = fill), notch = FALSE, varwidth = TRUE) +
               facet_grid(~ Group) +
               eval(parse(text = paste(fcol, '(', fparam, '=cols', ')'))) + labs(colour = fill) +
               annotate("text", y = 0, x = unique(.$x),
                        label = paste0("n = ",
                                       data.frame(a = .$variable, b = .$number) %>% unique() %>% .[['b']]),
                        hjust = 0) + 
               ggtitle(paste(label, subset_name, sep='')) + coord_flip() + ylab(metric) + xlab('') +
               theme_bw() + theme(legend.position= "right", aspect.ratio = 1)
                )
        preview(lapply(a$a, plot))

        if(!is.null(combination)){
            apply(combination, 2, function(x){
                      pairwise_comp = pausing_subset %>% 
                          filter(str_detect(variable,
                                            paste0("^", x[1], "\\-{0,1}[^\\.]*$")) | str_detect(variable,
                                            paste0("^", x[2], "\\-{0,1}[^\\.]*$"))) %>% droplevels()
                      pcumulative = pcumulative %+% pairwise_comp
            plot(pcumulative)
            return(x)})
        }
    }# }}}

}# }}}

# Calculating qqline # {{{
get_qqline = function(v){
    # from http://stackoverflow.com/questions/4357031/qqnorm-and-qqline-in-ggplot2
    y = quantile(v[!is.na(v)], c(0.25, 0.75))
    x = qnorm(c(0.25, 0.75))
    slope = diff(y)/diff(x)
    int = y[1L] - slope * x[1L]
    return(list(slope = slope, int = int))
}# }}}

# Check for normality via plots and gofstats -- interpret that with caution# {{{
check_normality = function(x, id = ''){
    library(MASS)
    library(e1071)
    library(raster)
    # for finding distribution of the data
    library(fitdistrplus)
    library(gridExtra)
    # Remove NA values as descdist complains - affects the qq-lot but how exactly
    if(any(is.na(x))) warning(sprintf("NA values not allowed, removing %s. Interpret results with caution",
                                      table(is.na(x))[['TRUE']]))

    x = x[!is.na(x)]
    gstats = apply(t(bind_cols(as.data.frame(t(as.matrix(test))),
                       data.frame(sd = sd(x),
                                  # type = 2 refers to wquation used to compute kurtosis, only type 2 is 
                                  # unbiased under normality
                                  kurtosis = kurtosis(x, type = 2),
                                  skewness = skewness(x, type = 2),
                                  cv = cv(x)))), 
                   2, round, 2)

    descdist(x, discrete = FALSE)
    # Plotting histogram, qqplot, pp-plot and ecdf
    check_against = c("norm", "lnorm", "pois", "exp", "gamma",
                      "nbinom", "geom", "beta", "unif", "logis")
    fits = lapply(check_against, function(d) tryCatch(fitdist(x, d, discrete = F), error = function(e) e))
    names(fits) = check_against
    fits = fits[which(unlist(lapply(fits, function(x) !is(x, 'error'))))]
    lapply(names(fits), function(x) {
                       par(oma = c(0, 0, 3, 0))
                       plot(fits[[x]])
                       title(paste(id, x), outer = TRUE)
                       reset_par()})

    gofstats = gofstat(fits)
    tmp_1 = t(apply(do.call(rbind, gofstats[c(1,3:4, 6, 8, 10, 12:13)]), 1, function(x) round(as.numeric(x), 2)))
    colnames(tmp_1) = names(gofstats$chisq)
    tmp_2 = do.call(rbind, gofstats[c(7,9,11)])
    gofstats = rbind(as.data.frame(tmp_1), as.data.frame(tmp_2))

    return(list(gstats = gstats, gofstats = gofstats, fits = fits))
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
    library(gridExtra)
    library(ggplot2)

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
        g(fc_annotated, targets) %=% tryCatch(annotate_df(fc),
                       error=function(e) e)
        if(is(fc_annotated, 'error')) {
            warning('Could not annotate the data frame (no file written).')
        } else {
            write.table((fc_annotated %>% select(-annotation_files)),
                        file = file.path(output_path, paste0('travelling-ratios-fc-annotated', dlab, '.tsv')),
                        quote = F, row.names = F, sep ="\t")
        }
    }

    # Melt travelling ratios and order according to times if samples from timecourse - gather() == melt()
    TR = travelling_ratio %>% select(gene_id, gene_name, ends_with("travelling_ratio_filtered")) %>%
        gather(variable, value, -gene_id, -gene_name)
    TR = format_melted(TR)

    pdf(file.path(plot_path, paste0('QQ-plots-normality.pdf')))# {{{
    ablines = lapply(levels(TR$variable), function(x) {
                         as.data.frame(get_qqline((TR %>% filter(variable == x) %>% .[['value']]))) %>%
                             mutate(variable = x) }) %>%
                Reduce(function(df1, df2) bind_rows(df1, df2), .) %>%
                full_join((TR %>% dplyr::select(Replicate, Time, variable) %>% group_by(variable) %>% dplyr::slice(1:1)), by = "variable")

    pp = ggplot(TR, aes(sample = value)) + stat_qq(geom = "point", size = 1) 
    pp = pp + geom_abline(data = ablines, aes(slope = slope, intercept = int), colour = "red")
    pp = pp + facet_grid(Time ~ Replicate)
    pp = pp + theme_bw() + theme(legend.position= "right", aspect.ratio = 1)
    plot(pp)

    gstats = lapply(levels(TR$variable), function(x){
                        tmp = TR %>% filter(variable == x) %>% .[['value']]
                        g(gstats, gofstats, fits) %=% check_normality(tmp, id = x)
                        plot.new()
                        grid.table(gofstats)
                        return(as.data.frame(t(gstats)))}) %>% Reduce(function(df1, df2) bind_rows(df1, df2), .)
    gstats = gstats %>% mutate(sample = levels(TR$variable))
    write.table(gstats,
                file = file.path(output_path, paste0('per_replicate_stats.tsv')),
                quote = F, row.names = F, sep ="\t")
    dev.off()# }}}

    pdf(file.path(plot_path, 'density-plots.pdf'))
    plot_density(TR, combination = combn(mixedsort(conditions),2), label = 'Full dataset')
    dev.off()
    
    # per_condition_tr is of class "tbl_df" and melt doesn't play nice with that
    tmp = per_condition_tr %>%
          gather(variable, value, -gene_id) %>%
          separate(variable, c("variable", 'x'), sep = '-') %>% select(-x)
    tmp = format_melted(tmp)
    
    pdf(file.path(plot_path, paste0('density-plots-mean', dlab, '.pdf')))
    plot_density(tmp, combination = combn(mixedsort(conditions),2), label = paste0('Mean per condition', dlab))
    dev.off()

    valid = per_condition_tr %>% na.omit() %>% gather(variable, value, -gene_id)
    valid = valid %>% separate(variable, c("variable", 'x'), sep = '-') %>% select(-x)
    valid = format_melted(valid)
    
    pdf(file.path(plot_path, paste0('density-plots-mean-common', dlab, '.pdf')))
    plot_density(valid, combination = combn(mixedsort(conditions),2), label = paste0('Mean per condition-common in all', dlab))
    dev.off()
    
    subsets = targets %>% mutate(Group = ifelse(is.na(Bound), 'Not Bound', as.character(Bound))) %>% select(-Bound)
    preview(plot_density(valid, subsets = subsets, subset_name = 'test'))

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
