#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)
source("/homes/mxenoph/source/Rscripts/granges-functions.R")

parser =  ArgumentParser(description="Annotate regions (overlap with genomic features and get genes bound by protein)")
parser$add_argument('-b', '--bed', metavar= "file", required='True', type= "character",
                    help= "BED file of regions to be annotated. If --computed is set then this file is expected to be of the same format as the *annotated.tsv files created from the script.")
parser$add_argument('-G', '--gtf', metavar= "file", type= "character",
                    default = "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.75.gtf",
                    help= "GTF file. Looks for the rtracklayer generated gtfs for that annotation. If not found it runs get-annotation-rscript.R")
parser$add_argument('-a', '--active_enhancers', metavar= "file", type= "character",
                    default = "/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-active_2015-08-10.bed",
                    help= "BED file of active enhancers")
parser$add_argument('-p', '--poised_enhancers', metavar= "file", type= "character",
                    default = "/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers-poised_2015-08-10.bed",
                    help= "BED file of active enhancers")
parser$add_argument('-o', '--out', metavar= "path", type= "character",
                    default= getwd(), help= "Output directory")
parser$add_argument('-c', '--computed', action ="store_true", default = FALSE,
                    help = "If set it just plots the values computed in the file. Default = FALSE")

args = parser$parse_args()

output_path = file.path(args$out, 'feature-annotation', basename(file_path_sans_ext(args$gtf)))
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

target = file_path_sans_ext(basename(args$bed))
# }}}

# Functions# {{{
library(rtracklayer)
library(dplyr)
library(ggplot2)

get_subsets = function(grange){# {{{
    sets = levels(factor(grange$name))
    if (length(sets) < 2) {
        print('According to the names provided, regions are not subsetted. Returning full grange.')
        return(grange)
    } else {
        subsets = lapply(sets, function(x){
                         pattern = glob2rx(paste(c("^", x, "$"), collapse=""))
                         grange[grepl(pattern, values(grange)[['name']])] })
        names(subsets) = sets
        return(subsets)
    }
}
# }}}

annotate = function(grange, annotations, output_path, target,# {{{
                    priority = c('2000-tss-500', 'first-exon', 'active_enhancers',
                                 'poised_enhancers', 'exons', 'first-intron', 'introns', 'genes')) {

    stopifnot(is(grange, "GRanges"))
    total = length(grange)
    per_feature = setNames(vector("list", length(priority)), priority)
    bound = data.frame(Gene = character(0), Bound = character(0))
    # '5000-tss-2000' contains all possible TSS for a gene
    for (i in priority) {
        ov = findOverlaps(annotations[[i]], grange)
        n_overlapping = length(unique(subjectHits(ov)))
        grange = grange[! (1:length(grange)) %in% unique(subjectHits(ov))]
        per_feature[[i]] = n_overlapping

        print(paste0('Calculating %peaks overlapping ', i, '.', sep=''))

        if(i %in% c('2000-tss-500', 'first-exon', 'exons', 'first-intron', 'introns', 'genes')){
            if (grepl('intron', i)) {
                fov = findOverlaps(annotations[[i]], annotations[['genes']])
                intron_in_genes = values(annotations[['genes']][unique(subjectHits(fov[queryHits(fov) %in% queryHits(ov)]))])[['ensembl_gene_id']]
                tmp = data.frame(Gene = intron_in_genes) %>% mutate(Bound = i)
            } else {
                if(grepl('exon', i)) {
                    x = 'ensembl_id'
                } else {
                    x = 'ensembl_gene_id'
                }
                # Even if queryHits are unique the Gene names might not be as grange has all possible TSS 
                # for a gene
                tmp = data.frame(Gene = unique(values(annotations[[i]][unique(queryHits(ov))])[[x]])) %>%
                        mutate(Bound = i)
            }
            # Find all genes that are not bound at a feature with higher priority
            # Avoiding duplicated values in bound$Gene
            tmp = anti_join(tmp, bound, by="Gene")
            bound = dplyr::bind_rows(bound, tmp)
        }
    }
    bound_table = data.frame(Gene = values(annotations[['genes']])[['ensembl_gene_id']]) %>% left_join(bound, by="Gene")
    write.table(bound_table,
                file = file.path(output_path, paste0(target,'.binding-per-gene.tsv')),
                quote = F, row.names = F, sep="\t")
    write.table(bound,
                file = file.path(output_path, paste0(target,'.bound-genes-only.tsv')),
                quote = F, row.names = F, sep="\t")

    if(per_feature[['genes']] != 0){
        stop('I find regions overlapping with genes even if tss, introns and exons have higher priority. Something went wrong')
    } else {
        per_feature = per_feature[which(names(per_feature)!='genes')]
    }

    per_feature[['intergenic']] = total - sum(unlist(per_feature))
    per_feature = unlist(per_feature)
    per_feature = per_feature / total
    return(per_feature)
}
# }}}

# plotting# {{{
plot_annotation  = function(anno, target){
    # if this is a matrix/data.frame it means it comes from the subset annotation
    if(is(anno, 'matrix') | is(anno, 'data.frame')) {
    #if((is(anno, 'matrix') | is(anno, 'data.frame')) & nrow(anno) > 1 ) {
        # if it has a column called regions that means that the function is called after reading an -annotations.tsv file
        if('Regions' %in% colnames(anno)){
            # keep Regions column and melt columns in two columns called Feature and value
            plot_data = anno %>% tidyr::gather(Feature, value, -Regions) %>%
            # because in the saved files values are in percentages
                mutate(value = value / 100)
        } else {
            plot_data = as.data.frame(anno) %>% add_rownames(var = 'Feature') %>%
                tidyr::gather(Regions, value, -Feature)
        }
    } else {
        plot_data = add_rownames(as.data.frame(anno), var = 'Feature') %>% rename(value = anno) %>%
            mutate(Regions = 1)
    }

    ordering = c('2000-tss-500','intergenic', 'first-exon', 'exons', 'first-intron', 'introns','active_enhancers', 'poised_enhancers')
    plot_data$Feature = factor(plot_data$Feature, levels = ordering)
    plot_data = plot_data %>% arrange(Feature) %>% mutate(value = value * 100)

    pdf(file.path(plot_path, paste0(target,'.annotated.pdf')))
    p = ggplot(plot_data, aes(x = Regions, y = value, fill = Feature)) + geom_bar(aes(width = 0.5), stat ="identity")
    p = p + scale_fill_brewer(palette = "Set2")
    p = p + theme_bw() + theme(legend.position= "right", aspect.ratio = 1, 
                                    axis.text.x = element_text(angle = 25, hjust = 1))
    p = p + labs(title = target, 
                 y = '% of peaks')

    if(length(unique(plot_data[['Regions']])) == 1) p = p + theme(axis.text.x = element_blank(),
                                                          axis.ticks = element_blank())
    plot(p)

    plot_data_reduced = plot_data %>% group_by(Regions) %>%
        filter(!Feature %in% c('exons', 'first-exon', 'introns', 'first-intron')) %>%
        do(bind_rows(., data.frame(Feature = 'gene_body', value = (100 - sum(.$value)), Regions = unique(.$Regions)))) %>%
        ungroup()

    ordering = c('2000-tss-500', 'gene_body', 'intergenic', 'poised_enhancers', 'active_enhancers')
    plot_data_reduced$Feature = factor(plot_data_reduced$Feature, levels = ordering)
    plot_data_reduced = plot_data_reduced %>% arrange(Feature)

    plot(p %+% plot_data_reduced)
    dev.off()

    return(plot_data)
}# }}}
# }}}

if(args$computed){
    anno = read.delim(args$bed)
    anno = anno %>% 
        select(-matches('File')) %>%
        rename_(.dots = setNames(colnames(.), 
                                 gsub("^X", '', gsub('\\.', '-', colnames(.)))))
    formatted_anno = plot_annotation(anno, target)
} else {
    pattern = paste0(file_path_sans_ext(args$gtf), '.rtracklayer-')
    rtracklayer_output = c('5000-tss-2000', '2000-tss-500',
                           'exons', 'first-exon',
                           'introns', 'first-intron',
                           'genes', 'canonical-transcripts')
    files = paste0(paste(pattern, rtracklayer_output, sep=''), '.gtf')

    if (any(file.exists(files))){
        print('Importing files.')
        annotations = lapply(files, function(x) import.gff3(x))
        names(annotations) = rtracklayer_output
        annotations$active_enhancers = import.bed(args$active_enhancers)
        annotations$poised_enhancers = import.bed(args$poised_enhancers)
    } else {
        command = paste("Rscript ~/source/get-annotation-rscript.R -G", args$gtf, '-p', ncores, collapse = " ")
        stop(paste('Annotation GTFs do not exist. Run `', command, '` first.', collapse = ' '))
    }

    if(file_ext(args$bed) != "bed") {
        if(file_ext(args$bed) == "narrowPeak") {
            regions = import_narrowPeak(args$bed)
        } else if (file_ext(args$bed) == "broadPeak"){
            target = paste0(target, '.broad-')
            regions = import_broadPeak(args$bed)
        } else if (file_ext(args$bed) == "gappedPeak"){
            target = paste0(target, '.gapped-')
            g(regions, blocks) %=% import_gappedPeak(args$bed)
        }
    } else {
        regions = import.bed(args$bed)
    }

    subsets = get_subsets(regions)

    if(is(subsets, "list")) {
        anno = sapply(subsets, function(x){
                      r = annotate(x, annotations, output_path, target)
                      return(r)})
    } else {
        anno = annotate(subsets, annotations, output_path, target)
    }

    formatted_anno = plot_annotation(anno, target)
    write.table((formatted_anno %>% tidyr::spread(Feature, value) %>%
                 dplyr::mutate(File = args$bed,
                        Active_enhancers_File = as.character(args$active_enhancers),
                        Poised_enhancers_File = as.character(args$poised_enhancers))),
                file.path(output_path, paste0(target,'.annotated.tsv')),
                quote = FALSE, row.names = F, sep ="\t")

    if(file_ext(args$bed) == "gappedPeak"){
        # if the file provided is gapped then no subsets provided and the name is the block name
        target = paste0(target, 'blocks-')
        blocks_anno = annotate(blocks, annotations, output_path, target)
        formatted_blocks = plot_annotation(blocks_anno, target)

        write.table((formatted_blocks %>% tidyr::spread(Feature, value) %>%
                     dplyr::mutate(File = args$bed,
                            Active_enhancers_File = as.character(args$active_enhancers),
                            Poised_enhancers_File = as.character(args$poised_enhancers)))
                    ,
                    file.path(output_path, paste0(target,'.annotated.tsv')),
                    quote = FALSE, row.names=F, sep="\t")
    }

# one I did for the TAC
# write.table(plot_data, '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/plots/NuRD-annotations.tsv', quote=FALSE, row.names=F, sep="\t")
}
