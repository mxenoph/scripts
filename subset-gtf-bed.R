#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
# Make library loading silent
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Get exons, introns, intergenic for given annotation")
parser$add_argument('-b', '--bed', metavar = "file", type = "character", help = "Annotation file (BED)")
parser$add_argument('-g', '--gtf', metavar = "file", type = "character", help = "Annotation file (GTF)")
parser$add_argument('-o', '--out', required ='True', type = "character", help = "Output directory")
parser$add_argument('-s', '--subsets', metavar = "file", type= "character",
                    help = "TSV with columns gene_id and Group")

args = parser$parse_args()

# For testing interactively# {{{
if(FALSE){
    args = list()
    args$bed = "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.75.rtracklayer-5000-gene_start-2000.bed6"
    args$out = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10"
    args$gtf = "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.75.rtracklayer-5000-gene_start-2000.gtf"
    args$subsets = "/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/pausing/Mean_per_condition-common_in_all-excluding-h24-S5P_1.valid-in-subset.tsv"
}# }}}

dir.create(args$out, recursive = TRUE)

# }}}

subsets = read.delim(args$subsets, header = T, sep = "\t")
subsets = lapply(split(subsets, subsets$Group), function(x) as.character(x[['gene_id']]))

if(is.null(args$bed) & is.null(args$gtf)){
    stop("No BED or GTF file provided")
} else {
    if(!is.null(args$bed)){
        bed = rtracklayer::import.bed(args$bed)
        names(bed) = bed$name
        sapply(names(subsets), function(x) {
                       tmp = bed[names(bed) %in% subsets[[x]]]
                       rtracklayer::export.bed(tmp, file.path(args$out, paste0(basename(file_path_sans_ext(args$bed)), '.', x, '.bed')))
                        })
    } 
    if(!is.null(args$gtf)) {
        gff = rtracklayer::import.gff3(args$gtf)
        names(gff) = gff$ensembl_gene_id
        print(head(gff))
        sapply(names(subsets), function(x) {
                       tmp = gff[names(gff) %in% subsets[[x]]]
                       rtracklayer::export.gff3(tmp, file.path(args$out, paste0(basename(file_path_sans_ext(args$gtf)), '.', x, '.gtf')))
                        })
    }
}

