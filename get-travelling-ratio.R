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
                    help= "GTF file for genes created with get-annotation-rscript.R")
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


    annotate_peaks = function(rangedData,name) {# {{{
            annotatedSample = annotatePeakInBatch(rangedData, AnnotationData = coord_Ranged,output = "overlapping")
            annotatedSample = as.data.frame(annotatedSample)
            annotatedSample = annotatedSample[,c("space","start","end","width","names", "peak","strand","feature","start_position","end_position","insideFeature","distancetoFeature","shortestDistance","fromOverlappingOrNearest")]
            annotatedSample = subset(annotatedSample,!is.na(insideFeature))
            annotatedSample = cbind(annotatedSample, read.table(text = as.character(annotatedSample$feature), sep = ":"))
            colnames(annotatedSample)[c(15,16)] = c("gene_id","status")

            ## find interval
            annotatedSampleinterval_start = 0
            annotatedSample$interval_start = ifelse(annotatedSample$start>annotatedSample$start_position,annotatedSample$start,annotatedSample$start_position)
            annotatedSample$interval_end = 0
            annotatedSample$interval_end = ifelse(annotatedSample$end<annotatedSample$end_position,annotatedSample$end,annotatedSample$end_position)

            ## if a peak overlaps with the prom and gene body of a gene, select the prom
            annotatedSample$peak_gene_id = paste(annotatedSample$peak,annotatedSample$gene_id,sep = " ")
            index = annotatedSample[duplicated(annotatedSample[,c("peak_gene_id")]),"peak_gene_id"]
            annotatedSampleDuplicated = annotatedSample[annotatedSample$peak_gene_id %in% index,]
            annotatedSamplenotDuplicated = annotatedSample[!annotatedSample$peak_gene_id %in% index,]
            annotatedSample = rbind(annotatedSamplenotDuplicated,subset(annotatedSampleDuplicated,status=="prom"))

            ## a peak could be associated to two genes
            annotatedSample = annotatedSample[,-c(19)]
            assign(paste0(name,"_annotated"),annotatedSample,envir = .GlobalEnv)
            assign(paste0(name,"_gene_body"),subset(annotatedSample,status=="gene_body"),envir = .GlobalEnv)
            assign(paste0(name,"_prom"),subset(annotatedSample,status=="prom"),envir = .GlobalEnv)
    }# }}}

#PRR = promoter release ratios from Fei Xavier Chen 2015
calculate_PRR = function(rpms, tss, peaks){

}

calculate_pausing_index = function(){
    promoter_avg / gene_body_avg
}

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

# Plotting density# {{{
plot_density = function(pausing, subsets = NULL, subset_name = NULL, label = '', metric = 'Travelling Ratio'){
    get_valid_n = function(pausing){
        n_valid = pausing %>% group_by(variable) %>%
            dplyr::summarise(number = table(is.na(value))['FALSE'])
        # Remove NA values
        pausing = pausing %>% filter(!is.na(value))
        return(list(pausing = pausing , n_valid = n_valid))
    }

    if(! 'Gene' %in% colnames(pausing)){
        pausing = add_rownames(as.data.frame(pausing), var = "Gene")
    }
    g(pausing, n_valid) %=% get_valid_n(pausing)

    # Plot basics# {{{
    p = ggplot(data = pausing)
    p = p + theme_bw() + theme(legend.position= "right", aspect.ratio = 1)
    p = p + ggtitle(label) + xlab(metric) + labs(colour = 'Sample')
    
    aesthetics_linetype = NULL
    if('Replicate' %in% colnames(pausing)){
        aesthetics_linetype = 'Replicate'
    }

    if('Time' %in% colnames(pausing)){
        divergent_colors_seq = setNames(colorRampPalette(brewer.pal(11,"PRGn")[c(1:4,8:11)])(length(levels(pausing$Time))),
                                  levels(pausing$Time))
        names(divergent_colors_seq) = levels(pausing$Time)
        
        pviolin = p + geom_violin(aes_string(x = 'variable', y = 'value', fill = 'Time'))
        pviolin = pviolin + scale_fill_manual(values = divergent_colors_seq) + labs(fill = 'Time')

        p = p + geom_line(stat = "density", aes_string(x = 'value', color = 'Time', linetype = aesthetics_linetype))
        p = p + scale_colour_manual(values = divergent_colors_seq) + labs(colour = 'Time')

    } else{
        pviolin = p + geom_violin(aes_string(x = 'variable', y ='value', fill = 'variable')) + labs(fill = 'Sample')
        p = p + geom_line(stat = "density", aes_string(color = 'variable', linetype = aesthetics_linetype))

        if(length(levels(pausing$variable)) > 9){
            default_colors= setNames(colorRampPalette(brewer.pal(8,"Set2"))(length(levels(pausing$variable))),
                                      levels(pausing$variable))
            names(default_colors) = levels(pausing$variable)

            pviolin = pviolin + scale_fill_manual(values = default_colors)
            p = p + scale_colour_manual(values = default_colors)
        } else{
            pviolin = pviolin + scale_fill_brewer(palette = "Set2")
            p = p + scale_colour_brewer(palette = "Set2")
        }
    }

    # sample column to be used in violinplot so I can tweak the geom_text positions
    annotations = pausing %>% 
        mutate(Sample = mapvalues(variable, from = levels(variable), to = 1:length(levels(variable)))) %>% 
        select(variable, Sample, Time) %>% unique() %>%
        inner_join(n_valid, by = 'variable')
    pviolin = pviolin + annotate("text",
                                 x = annotations$Sample %>% as.numeric() + 0.35,
                                 y = 0,
                                # left justified
                                 hjust = 0,
                                 label = paste0("n = ", annotations$number))
    pviolin = pviolin + coord_flip() + ylab(metric) + xlab('')

    plot(p)
    plot(pviolin)# }}}

    # density = counts / sum(counts * bar width)
    if (!is.null(subsets)){
        pausing_subset = pausing %>% inner_join(subsets, by = "gene_id") %>% droplevels()
        q = p %+% pausing_subset + facet_grid(.~ Group) + ggtitle(paste(label, subset_name, sep=''))
        qviolin = pviolin %+% pausing_subset + facet_grid(.~ Group) + ggtitle(paste(label, subset_name, sep=''))

        plot(q)
        plot(qviolin)
    }
    
}# }}}

main(){
    library(dplyr)
    library(tidyr)
    # for extracting pattern from string with str_extract
    library(stringr)
    library(rtracklayer)

    files = list.files(pattern=glob2rx(paste("*_gene_signal*",".txt",sep="")), args$bed, full.name = T)
    travelling_ratio = calculate_travelling_ratio_meryem(files)
    TR = melt(travelling_ratio %>% select(gene_id, gene_name, ends_with("travelling_ratio_filtered")))

    # order according to times if samples from timecourse# {{{
    contains_times = any(gregexpr("h[[:digit:]]+", TR$variable, ignore.case = T) != -1) | any(gregexpr("[[:digit:]]+h", TR$variable, ignore.case = T) != -1) | any(gregexpr("[[:digit:]]+m", TR$variable, ignore.case = T) != -1)
    if(contains_times){
        TR = TR %>% separate(variable, c("Sample", "Replicate"), sep = "-") %>% 
            separate(Replicate, c("Replicate", "tmp"), sep = "\\.") %>% 
            mutate(Time = gsub('_', '', str_extract(Sample, "_.*")),
                   Sample = gsub('_', '', str_extract(Sample, ".*_")),
                   Mark = gsub('_', '', str_extract(Replicate, "[^_]*_")),
                   Replicate = gsub('_', '', str_extract(Replicate, "_\\d_"))) %>%
            select(-tmp)
        
        tmp = convert_time(TR$Time)
        TR[tmp$index, 'Time'] = tmp$values
        TR$Time = factor(TR$Time, levels = mixedsort(unique(TR$Time)))
        TR = TR %>% mutate(variable = paste0(Sample,'_', Time, '-', Mark, '_', Replicate)) %>%
            select(-Sample, -Mark)
        TR$variable = factor(TR$variable, levels = mixedsort(unique(TR$variable)))
        TR = TR %>% arrange(Time, variable)
    } else{
        TR$variable = factor(TR$variable, levels = mixedsort(unique(TR$variable)))
        TR = TR %>% arrange(variable)
    }# }}}

    preview(plot_density(TR))
    subsets = travelling_ratio[1:100, 1, drop=F]
    subsets = subsets %>% mutate(Group = c(rep('A', 50), rep('B', 50)))
    preview(plot_density(TR, subsets = subsets))

    pdf(file.path(plot_path, 'cumulative-distribution.pdf'))
    p = ggplot(pausing, aes(value, colour = variable)) + stat_ecdf(geom = "point") + theme_bw() + theme(legend.position= "bottom") + ggtitle('all genes')
    plot(p)

    combination = combn(c('h0', '30m', 'h4', 'h24'), 2)

    apply(combination, 2, function(x){
    subsets = pausing %>% filter(grepl(glob2rx(paste0("*", x[1], "*")), variable) | grepl(glob2rx(paste0("*", x[2], "*")), variable))
    p = ggplot(subsets, aes(value, colour = variable)) + stat_ecdf(geom = "point") + theme_bw() + theme(legend.position= "bottom") + ggtitle('all genes')
    plot(p)
                                               return(x)})
    dev.off()

    pausing_long = pausing %>% filter(gene_id %in% genes_long$gene_id)
    cdf = lapply((pi.table %>% select(-gene_id, -gene_name)), function(x) {x = x[!is.na(x)]; return(list(cdf = ecdf(x), l = length(x)))})

    pdf(file.path(plot_path, 'cumulative-distribution-long.pdf'))
    p = ggplot(pausing_long, aes(value, colour = variable)) + stat_ecdf(geom = "point") + theme_bw() + theme(legend.position= "bottom") + ggtitle('all genes')
    plot(p)

    combination = combn(c('h0', '30m', 'h4', 'h24'), 2)

    apply(combination, 2, function(x){
    subsets = pausing_long %>% filter(grepl(glob2rx(paste0("*", x[1], "*")), variable) | grepl(glob2rx(paste0("*", x[2], "*")), variable))
    p = ggplot(subsets, aes(value, colour = variable)) + stat_ecdf(geom = "point") + theme_bw() + theme(legend.position= "bottom") + ggtitle('all genes')
    plot(p)
                                               return(x)})
    dev.off()






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
}
