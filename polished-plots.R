#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library = function (...) suppressMessages(base::library(...))
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
select = dplyr::select
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")
source("~/source/Rscripts/granges-functions.R")

parser =  ArgumentParser(description="Polished thesis plots")
parser$add_argument('--f2.1', type= "character", help= "Plot figure 1 in chapter 2")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")
args = parser$parse_args()

# For testing interactively# {{{
if(FALSE){
    args = list()
    args$f2.1 = 'true'
    args$out = "/nfs/research2/bertone/user/mxenoph/hendrich/thesis"
}# }}}

plot_path = file.path(args$out, 'plots')
dir.create(plot_path, recursive= TRUE)
chip_path = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/"
rna_path = "/nfs/research2/bertone/user/mxenoph/hendrich/rna/"

gg_param = modules::import('ggplots')

f2.1 = function(){# {{{
    tsv_files = data.frame(file = list.files(pattern = ".*annotated.tsv",
                                             file.path(chip_path,
                                                       "hendrich_2013/mm10/feature-annotation"),
                                             full.names = TRUE)) %>%
                mutate( tmp = gsub('^.*/', '', file)) %>% separate(tmp, c('desc', 'tmp'), '_peaks') %>%
                separate(desc, c('genotype', 'condition', 'rep'), '_') %>%
                separate(condition, c('condition', 'protein'), '-') %>%
                filter(protein %in% c('M2', 'M12b', 'Chd4'), condition %in% c('2i', 'EpiSC')) %>%
                group_by(condition, protein) %>% do(filter(., rep == 'pooled' | is.na(rep) )) %>%
                ungroup()

    plot_data = lapply(tsv_files %>% .$file %>% as.character(), read.delim) %>%
        Reduce(function(df1, df2) bind_rows(df1, df2), .) %>% 
        mutate( tmp = gsub('^.*/', '', File)) %>% separate(tmp, c('desc', 'tmp'), '_peaks') %>%
        separate(desc, c('genotype', 'condition', 'rep'), '_') %>%
        separate(condition, c('condition', 'protein'), '-') %>%
        mutate(condition = ifelse(condition == '2i', 'ESC', 'EpiSC')) %>%
        mutate(protein = ifelse(protein == 'M2', 'Mbd3', 'Chd4')) %>%
        select(-Regions, -File, -tmp, -rep)

    plot_data_reduced = plot_data %>% 
        filter(genotype != '7G9') %>%
        group_by(genotype, condition, protein) %>%
        gather('feature', 'value', 1:8) %>%
        filter(!feature %in% c('exons', 'first.exon', 'introns', 'first.intron')) %>%
        do(bind_rows(., data.frame(feature = 'gene_body', value = (1 - sum(.$value)), 
                                   genotype= unique(.$genotype),
                                   condition= unique(.$condition),
                                   protein= unique(.$protein)))) %>%
        ungroup()
            
    ordering = c('X2000.tss.500','gene_body', 'active_enhancers', 'poised_enhancers', 'intergenic')
    plot_data_reduced$feature = factor(plot_data_reduced$feature, levels = ordering)
    plot_data_reduced = plot_data_reduced %>% arrange(feature) %>% mutate(value = value * 100)
    
    p = ggplot(plot_data_reduced, aes(x = protein, y = value,
                                      fill = feature, group = genotype)) + facet_grid(.~condition)
    p = p + geom_bar(aes(width = 0.5), stat ="identity")
    p = p + scale_fill_brewer(palette = "Set2",
                              labels = setNames(c('Proximal promoter','Gene body',
                                                  'Active enhancers', 'Inactive enhancers', 'Intergenic'), ordering))
    p = p + gg_param$theme_publication + labs(x = '', y = "% of peaks", fill = 'Feature')

    pdf(file.path(plot_path, paste0("f2.1a.pdf")))
    plot(p)
    dev.off()

}
# }}}

f2.2 = function(){# {{{
    merged_files = data.frame(file = list.files(pattern = ".*Ov1bp_complete_annotated.bed",
                                 file.path(chip_path, 'hendrich_2013/mm10/merged-ranges/'),
                                 full.names = TRUE)) %>%
                mutate( tmp = gsub('^.*/', '', file)) %>% separate(tmp, c('condition', 'protein'), '_Proteins_') %>%
                mutate(condition = gsub('Condition_', '', condition),
                       condition =  gsub('_AND_', ':', condition),
                       protein = gsub('_Ov1bp.*', '', protein),
                       protein = gsub('_AND_', ':', protein)) %>%
                filter(grepl(".*_filtered.*", protein)) %>%
                mutate(protein = gsub('_filtered', '', protein),
                       reps = protein,
                       protein = gsub("_[^:]*", '', protein),
                       reps = ifelse(grepl("_", reps), gsub("[^:]*_", '', reps), NA)) %>%
                group_by(condition, protein) %>%
                # reps are numeric except the 'pooled' one which is the first value if desc() used
                arrange(desc(protein)) %>%
                filter(grepl("M2|Chd4", protein))

    merged_ranges = lapply(setNames(merged_files %>% .$file %>% as.character,
                                    merged_files %>% .$condition %>% as.character),
                           function(x){
                               as.data.frame(rtracklayer::import(x)) %>%
                                   rename(score = 'fdr', seqnames = 'chr')
                           })
    
    tab_files = data.frame(file = list.files(pattern = ".*.tab",
                                             file.path(chip_path,
                                                       "hendrich_2013/mm10/macs2/sharp"),
                                             full.names = TRUE)) %>%
                mutate( tmp = gsub('^.*/', '', file)) %>%
                separate(tmp, c('condition', 'protein', 'x', 'y'), '-') %>%
                mutate(id = condition) %>%
                select(-c(x, y)) %>%
                separate(condition, c('genotype', 'condition'), '_') %>%
                mutate(protein = gsub('_logFE.*', '', protein))
    
    tab_df = lapply(setNames(tab_files %>% .$id %>% unique,
                             tab_files %>% .$id %>% unique),
                    function(x){
                        y = setNames(tab_files %>% filter(id == x) %>% .$file %>% as.character,
                                     tab_files %>% filter(id == x) %>% .$protein %>% as.character)
                        lapply(names(y),
                               function(z){
                                   read.table(y[z], sep = "\t", skip=1,
                                              col.names = c('chr', 'start', 'end', 'score')) %>%
                                   # tab file has different index than imported bed files
                                   mutate(start = start + 1) %>%
                                   mutate(condition = x,
                                          protein = z) %>%
                                   inner_join(merged_ranges[[x]], by=c('chr', 'start', 'end'))
                               }) %>% Reduce(function(df1, df2) bind_rows(df1, df2), .)
                    }) %>%
                Reduce(function(df1, df2) bind_rows(df1, df2), .)

    tab_df = tab_df %>% mutate(tmp = gsub("M2", "Mbd3",
                                                         gsub("_filtered", '',
                                                              gsub("^[^-]*-", '', gsub("_AND[^-]*", '', name))))) %>%
    mutate(tmp = ifelse(grepl('-', tmp), tmp, paste(tmp, 'unique sites', sep = " ")), 
           condition = gsub(".*_", '', condition), 
           protein = gsub('M2', 'Mbd3', protein))

    p = ggplot(tab_df, aes(x = tmp, y = score, fill = protein)) + geom_violin()
    p = p + facet_grid(.~condition)
    p = p + gg_param$theme_publication + labs(y = "Log Fold Enrichment", fill = 'ChIP-seq', x = '')
    preview(print(p))
    
    pdf(file.path(plot_path, paste0("f2.2.pdf")))
    plot(p)
    dev.off()

}
# }}}

#}}}
