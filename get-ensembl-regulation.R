#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Get exons, introns, intergenic for given annotation")
parser$add_argument('-G', '--gtf', metavar= "file", required='True', type= "character", help= "GTF file")
parser$add_argument('-t', '--etable', metavar= "file",
                    default = 'http://Sep2015.archive.ensembl.org/info/website/archives/assembly.html',
                    type= "character", help= "ENSEMBL table with month of version release")
parser$add_argument('-p', '--threads', required='False', default = 4, help= "Number of processors to use")
parser$add_argument('-f', '--filter', required='True', help= "keyword to filter by, should be included in the scripts dictionary")

args = parser$parse_args()
ncores = args$threads

if(FALSE){
    args = list()
    args$gtf = "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.75.gtf"
    args$etable = "http://Sep2015.archive.ensembl.org/info/website/archives/assembly.html"
}# }}}
# }}}

get_ensembl_regulatory = function(filters, values, archive = NULL, organism = NULL, output = NULL, gtf_eversion = NULL){ # {{{

    # Load packages# {{{
    library(stringr)
    # }}}

    # Ensembl regulatory build was introduced in e73 from September 2013. 
    # If archive version requested is older stop the script with a message
    if(as.numeric(gtf_eversion) < 73) stop('Archived ensembl version does not have regulatory build')

    attr_to_get = list('regulatory_feature' =  c("regulatory_stable_id",
                                                 "bound_seq_region_start", "bound_seq_region_end",
                                                 "chromosome_name", "chromosome_start", "chromosome_end",
                                                 "feature_type_name", "so_accession"),
                                                 #"cell_type_name", "cell_type_description")
                       'feature_set' = c("regulatory_stable_id",
                                         "bound_seq_region_start_1051", "bound_seq_region_end_1051",
                                         "seq_region_name_1057", "seq_region_start_1057", "seq_region_end_1057",
                                         "feature_type_name_1057", "display_label_1057",
                                         "cell_type_name_1051", "cell_type_description_1051")
                       )

    dataset_prefix = 'regulatory_feature'
    # Dataset names change from '_feature_set' in Dec 2014 to '_regulatory_feature' in Mar 2015
    if(as.numeric(gtf_eversion) < 79) dataset_prefix = 'feature_set'

    so_accession = data.frame(so_accession = c('SO:0001952', 'SO:0001974',
                                               'SO:0000235', 'SO:0000167',
                                               'SO:0001747', 'SO:0000165'),
                              feature_type_name = c('Promoter Flanking Region', 'CTCF Binding Site',
                                                    'TF binding site', 'Promoter',
                                                    'Open chromatin', 'Enhancer'))

    x = lapply(paste(c('sep2013', 'dec2013', 'feb2014', 'aug2014', 'oct2014', 'dec2014', 'mar2015', 'may2015', 'jul2015', 'sep2015'),
                     'archive.ensembl.org', sep='.'),
               function(x) {
                   dataset_prefix = 'feature_set'
                   if(grepl('2015', x))  dataset_prefix = 'regulatory_feature'
                   ensembl = useMart(host = x, biomart = 'ENSEMBL_MART_FUNCGEN', dataset = paste(organism, dataset_prefix, sep = '_'))
#                   regulatory_features = getBM(attributes = attr_to_get[[dataset_prefix]],
#                                               mart = ensembl)
               })
    y = lapply(x, function(i){
                   regulatory_features = getbm(attributes = attr_to_get[[dataset_prefix]],
                                               mart = ensembl)
               })

    ensembl = useMart(host = archive, biomart = 'ENSEMBL_MART_FUNCGEN', dataset = paste0(organism, dataset_prefix, sep = '_'))
    regulatory_features = getBM(attributes = attr_to_get[[dataset_prefix]],
                                mart = ensembl)

    # Removing the DNA binding proteins coming from PF00125
    # I have to use dplyrr::filter because filtered is masked by one of the bioc packages
    # stringr::coll() compares strings respecting standard collation rules
    filtered_features = filtered_features %>% dplyr::filter(str_detect(description, coll(args$filter, ignore_case = T)))
    
    write.table(filtered_features, output, sep="\t", quote=FALSE, col.names=T, row.names=F)
    cat("Returning requested genes \n")
    return(filtered_features)
}# }}}

main = function(){# {{{
    library(rvest)
    library(dplyr)
    library(tidyr)
    library(biomaRt)
    
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
    dictionary = list('histone' = list(filters = c("?"), values  = c("embryonic stem cells")))

    if(any(grepl(args$filter, names(dictionary)))){
        index = grep(args$filter, names(dictionary))
        target = paste0(target, '.ensembl-regulation', names(dictionary)[index], '-genes.tsv')

        filtered_features = get_ensembl_regulatory(archive = archive,
                                     organism = organism,
                                     output = target,
                                     filters = dictionary[[index]]$filters,
                                     values = dictionary[[index]]$values,
                                     gtf_eversion = gtf_eversion)
    } else {
        warning("The filter provided is not included in the dictionary. Update script and run again\n")
    }

}
# }}}

main()
