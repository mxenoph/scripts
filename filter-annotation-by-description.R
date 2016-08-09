#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)
source("/homes/mxenoph/source/Rscripts/granges-functions.R")

parser =  ArgumentParser(description="Get exons, introns, intergenic for given annotation")
parser$add_argument('-G', '--gtf', metavar= "file", required='True', type= "character", help= "GTF file")
parser$add_argument('-t', '--etable', metavar= "file",
                    default = 'http://Sep2015.archive.ensembl.org/info/website/archives/assembly.html',
                    type= "character", help= "ENSEMBL table with month of version release")
parser$add_argument('-p', '--threads', required='False', default = 4, help= "Number of processors to use")
parser$add_argument('-f', '--filter', required='True', help= "keyword to filter by, should be included in the scripts dictionary")

args = parser$parse_args()
ncores = args$threads
# }}}

get_filtered_features = function(filters, values, archive = NULL, organism = NULL, output = NULL){ # {{{

    # Load packages# {{{
    library(GenomicFeatures)
    library(GenomicRanges)
    library(rtracklayer)
    library(parallel)
    library(BiocParallel)
    library(stringr)
    # }}}

    ensembl = useMart(host = archive,
                      biomart = 'ENSEMBL_MART_ENSEMBL',
                      dataset = paste0(organism, "_gene_ensembl"))
    
    filtered_features = getBM(attributes = c("ensembl_gene_id", "external_gene_id",
                                               "description"),
                                filters = filters,
                                # from PFAM "PF00125" are the core histones together with 
                                #  some other DNA binding proteins form a superfamily defined 
                                # by a common fold and distant sequence similarities 
                                # PF00538is for the linker histones
                                values = values,
                                mart = ensembl)

    # Removing the DNA binding proteins coming from PF00125
    # I have to use dplyrr::filter because filtered is masked by one of the bioc packages
    # stringr::coll() compares strings respecting standard collation rules
    filtered_features = filtered_features %>% dplyr::filter(str_detect(description, coll(args$filter, ignore_case = T)))
    
    write.table(filtered_features, output, sep="\t", quote=FALSE, col.names=T, row.names=F)
    print('Returning requested genes')
    return(filtered_features)
}
# }}}

main = function(){# {{{
    library(rvest)
    library(dplyr)
    library(tidyr)
    
    doc_text = read_html(args$etable) %>% html_nodes('th') %>% html_text()
    ensembl_versions = do.call(rbind, (strsplit(doc_text[grepl("\\d+", doc_text, perl=TRUE)], ' ' )))
    ensembl_versions = as.data.frame(ensembl_versions) %>% dplyr::rename(month = V1) %>%
                        separate(V2, c('year', 'version'), sep = 'v') %>% 
                        mutate(archive = paste0(tolower(month), year,".archive.ensembl.org")) %>%
                        distinct(version)

    # Build host name e.g for ensembl70 host='jan2013.archive.ensembl.org', ds= 'mmusculus_gene_ensembl'
    gtf_name = basename(args$gtf)

    # Converting e.g. Mus_musculus to mmusculus
    organism = unlist(strsplit(gsub("\\..*", "", gtf_name), '_'))
    organism[1] = substr(organism[1], 1, 1)
    organism = tolower(paste(organism, collapse=''))

    pattern = "\\.[[:digit:]]*\\."
    m = gregexpr(pattern, gtf_name)
    gtf_eversion = gsub("\\.", "", unlist(regmatches(gtf_name, m)))
    archive = ensembl_versions %>% filter(version == gtf_eversion) %>% dplyr::select(archive) %>% .[[1]]

    datasets_versions = listDatasets(useMart(host = archive,
                                             biomart = 'ENSEMBL_MART_ENSEMBL')) %>%
                        separate(dataset, into = c("Organism", "suffix_1", "suffix_2"), sep = "_") %>%
                        dplyr::select(-suffix_1, -suffix_2)
    
    target = gsub(".gtf", "", args$gtf)

    dictionary = list('histone' = list(filters = c("pfam"), values  = c("PF00125", "PF00538")))

    if(any(grepl(args$filter, names(dictionary)))){
        index = grep(args$filter, names(dictionary))
        target = paste0(target, '.filtered-', names(dictionary)[index], '-genes.tsv')

        filtered_features = get_filtered_features(archive = archive,
                                     organism = organism,
                                     output = target,
                                     filters = dictionary[[index]]$filters,
                                     values = dictionary[[index]]$values)
    } else {
        warning("The filter provided is not included in the dictionary. Update script and run again\n")
    }

}
# }}}

main()
