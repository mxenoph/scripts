#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Define enhancer elements")
parser$add_argument('-c', '--config', metavar= "file", required='True', type= "character", help= "config file")
parser$add_argument('-a', '--assembly', type= "character", default='mm10', help= "Give preferred assembly e.g. mm10. Default: mm10")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args = parser$parse_args()
output_path = file.path(args$out, 'enhancers')
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

#}}}

# Load packages# {{{
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPpeakAnno)
library(biovizBase)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(rtracklayer)# }}}

#Functions
#Caution: different from previous script: take intersection for data coming from different scripts
parse_peaks = function(config){# {{{
    source("~/source/Rscripts/granges-functions.R")
    
    peaks = lapply(config[['file']], function (fs){
           x = match.fun(paste('import', file_ext(fs), sep = '_'))(fs)
           if(file_ext(fs) == 'gappedPeak') x$blocks else x
        })

    names(peaks) = with(config, paste(system, factor, rownames(config), sep=':'))
    metadata_has_score = all(sapply(peaks, function(x) 'fe' %in% gtools::mixedsort(colnames(values(x)))))
    
    descriptors = unique(gsub(":[0-9]{1,2}$", '', names(peaks)))
    pooled = sapply(descriptors, function(x){
        i = grep(x, names(peaks))
        if(length(i) == 1) {
            return(peaks[[i]])
        }else {
            gr = intersect_multi(peaks[i])
            for(j in colnames(values(peaks[[1]])) ) values(gr)[j] = rep(NA, length(gr))

            if(metadata_has_score){
                scores = mclapply(peaks[i], function(x) get_hits_scores(gr, x, 'FE'))
                names(scores) = names(peaks[i])
                scores = do.call(cbind.data.frame, scores)
                scores[,'mean'] = apply(scores,1,mean)
                values(gr)['FE'] = scores[['mean']]
                #apply(scores,1,sd)
            }
            return(gr)
        }
    })
    return(list(raw_peaks=peaks, pooled=pooled))
}# }}}

config = read.table(args$config, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# list of granges grL[[1]] it's all the data while grL[[2]] it's pooled experiments for same mark
peaks_list = parse_peaks(config)

# .$factor selects only column named factor
active = config %>% filter(element == 'enhancer' & type =='active') %>% .$factor %>% unique
poised = config %>% filter(element == 'enhancer' & type =='poised') %>% .$factor %>% unique
promoter = config %>% filter(element == 'promoter' & type =='active') %>% .$factor %>% unique

# get index of datasets for each type of feature
get_index = function(type) unlist(mclapply(type, grep, names(peaks_list[['pooled']])))

promoter = intersect_multi(peaks_list[['pooled']][get_index(promoter)])
enhancers = mclapply(list('active' = active, 'poised' = poised),
                     function(x){
                         intersection = intersect_multi(peaks_list[['pooled']][get_index(x)])
                         ov = findOverlaps(intersection, promoter)
                         keep = (1:length(intersection))[! 1:length(intersection) %in% unique(queryHits(ov))]
                         intersection[keep]
                     })
# introduced in 27/9/2016
enhancers$poised = enhancers$poised[width(enhancers$poised) > 49]
enhancers$active = enhancers$active[width(enhancers$active) > 49]

# remove active enhancers from poised ones
ov = findOverlaps(enhancers$poised, enhancers$active)
enhancers$poised = enhancers$poised[!names(enhancers$poised) %in% names(enhancers$poised[unique(queryHits(ov))])]

# dealing with cases where metadata name and score are not set
enhancers = lapply(enhancers,
           function(x){ 
               if(all(is.na(values(x)[['name']]))) values(x)[['name']] = names(x)
               if(all(is.na(values(x)[['score']]))) values(x)[['score']] = rep(0, length(x))
               x
           })

validate = function (feature) {# {{{
    # from papers: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0114485
    # http://mcb.asm.org/content/25/14/6031.full
    known_enhancers = GRanges(c('chr3','chr17'),
                              IRanges(start=c(34752460, 35504494), end=c(34765460, 35505988)),
                              strand = rep('*',2),
                              name = c('Sox2-distal', 'Oct4-proximal')
                              )
    ov = findOverlaps(feature, known_enhancers)
    report = feature[queryHits(ov)]
    values(report)['name'] = c(known_enhancers[subjectHits(ov)]$name)
    report
}# }}}

export.bed(enhancers[['active']],
           file.path(output_path, paste0('enhancers-active_', Sys.Date(),'.bed')))

export.bed(enhancers[['poised']],
           file.path(output_path, paste0('enhancers-poised_', Sys.Date(),'.bed')))

format_data = function(feature) {
    cbind(as.data.frame(values(feature)[c('qvalue','fe')]),
          data.frame('width' = width(feature)/1000)
          )}

plot_data = lapply(names(enhancers), function(x) data.frame(region = names(enhancers[[x]]),
                                                    width = width(enhancers[[x]])/1000,
                                                    type = rep(x, length(enhancers[[x]])))) %>%
                    Reduce(function(df1, df2) bind_rows(df1, df2), .)
pd_vlines = plot_data %>% group_by(type) %>% do(summary = summary(.$width)) %>%
    tidy(summary) %>% full_join((plot_data %>% group_by(type) %>% summarize(n = n())), by = "type")

p = ggplot(plot_data, aes(x = width)) + geom_histogram(binwidth=0.1, pad=FALSE) + facet_grid(~type)
p = p + geom_vline(data = pd_vlines, aes(xintercept = median), colour = "red", linetype = "dotted")

annotation_y_position = sum(ggplot_build(p)$panel$ranges[[1]]$y.range)
annotation_x_position = max(ggplot_build(p)$panel$ranges[[1]]$x.range)

p = p + geom_text(data = pd_vlines, aes(median, annotation_y_position, label = paste0("median ==", median)),
                  angle = 270, vjust = -1,
                  hjust = -0.5,
                  colour = "red", parse = T)
p = p + geom_text(data = pd_vlines, aes(annotation_x_position, annotation_y_position, label = paste0("n ==", n)),
                  hjust = 1,
                  parse = T)

p = p + theme_bw() + xlab('Width (Kb)') + theme(aspect.ratio = 1)

pdf(file.path(plot_path, 'enhancer-widths.pdf'))
plot(p)
dev.off()
#peak_metrics = do.call(rbind,
#                       mclapply(peaks_list[['pooled']],
#                                format_data))

