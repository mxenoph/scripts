#!/usr/bin/env Rscript --slave

# Parsing command line arguments and create output subdirectories# {{{
library = function (...) suppressMessages(base::library(...))
library(argparse)
library(tools)

parser =  ArgumentParser(description="Run differential expression analysis")
parser$add_argument('-d', '--design', metavar = "file", required ='True', type= "character", help = "TSV file with columns File, Library Name, Genotype/Condition")
parser$add_argument('-o', '--out', metavar = "path", type = "character", default = getwd(), help = "Output directory -- all subdirectories will be created here")
parser$add_argument('-c', '--counts', metavar = "path", required = 'True', type = "character", help = "Path to look for htseq-count files")
parser$add_argument('-f', '--factors', action ="store_true", default = FALSE, help = "If set it just prints the scale factors and exits. Default = FALSE")
parser$add_argument('-n', '--fragment', default = 200, help = "Fragment length selected in library prep")
parser$add_argument('-g', '--gtf', default = "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.70.gtf",
                    help= "GTF file. Looks for the rtracklayer generated gtfs for that annotation. If not found it runs get-annotation-rscript.R")
parser$add_argument('-s', '--genome_size', default = FALSE, help= "Chromosome lengths file.")
parser$add_argument('-t', '--targets', default = FALSE, help= "TSV with gene Id column and Bound column generated from annotate-regions.R")
parser$add_argument('--tname', default = FALSE, help= "Name/description for the subsetted/targeted genes")
# e.g "/nfs/research2/bertone/user/mxenoph/common/genome/MM10/MM10.genome"

args = parser$parse_args()

output_path = file.path(args$out, 'deseq')
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

# helper function to plot in preview.pdf when debugging# {{{
preview = function(f) {
    pdf('~/public_html/preview.pdf')
    eval(f)
    dev.off()
}# }}}
# }}}}

x = c('DESeq2', 'dplyr')
lapply(x, library, character.only=T)
source("~/source/Rscripts/deseq-functions.R")
source("~/source/Rscripts/rna-norm-functions.R")
source("~/source/Rscripts/functions.R")
#select = dplyr::select

# dplyr helper function.
add_rownames = function(df, var = 'rowname') {
        stopifnot(is.data.frame(df))
    rowname_df = setNames(data_frame(rownames(df)), var)
        cbind(rowname_df, df)
}

counts_path = args$counts
design = read.delim(args$design, header = TRUE)

# Reading in counts for samples in design
counts = export_counts(counts_path)
lib_indices = match(design$Library, colnames(counts))
counts = counts[, lib_indices[!is.na(lib_indices)]]

# the last level of the Contrast variable will be used over the first level 
# e.g if levels(design$Contrast)
# [1] wt wt ko ko
# Levels: ko wt
# FC = log2(wt/ko) and positive FC corresponds to downregulation in the ko while positive corresponds to upregulation in the ko
# As this is counterintuitive when dealing with wt,ko change the levels order so that positive FC means upreg in ko
# always write files as denominator-nominator.deseq.tsv so you know what comparison has been made
if(all(levels(design$Contrast) == c('ko','wt'))) {
    design$Contrast = factor(design$Contrast, levels = c('wt', 'ko'))
}

cds = DESeqDataSetFromMatrix(counts, design, design = ~Contrast)
cds = estimateSizeFactors(cds)
# colnames are not set in the cds for some puzzling reason so setting them manually
colnames(cds) = colnames(counts)

# Check whether comparison is between treatment or cell type and set the label accordingly# {{{
comparison = design %>% select(Contrast, Condition) %>% unique()
if(comparison %>% select(Condition) %>% unique() %>% length() == 1){
    label = paste(comparison[1,'Condition'],
                  paste(levels(comparison[['Contrast']]), collapse ='-'),
                  sep = '_')
} else {
    label = comparison %>% mutate(Label = paste(Condition, levels(Contrast), sep='_')) %>% .[['Label']] %>% paste(collapse = '-')
}

index = lapply(1:nrow(comparison), function(i){
               design %>% inner_join(comparison[i,], by = c('Contrast', 'Condition')) %>% 
               dplyr::select(Library) %>%
               .[[1]] %>% as.character()
})
names(index) = comparison %>% mutate(Name = paste(Contrast, Condition, sep = "_")) %>% 
                dplyr::select(Name) %>% .[[1]] %>% as.character()
# }}}

write.table(as.data.frame(sizeFactors(cds)),
            file = file.path(output_path, paste0(label, '_size-factors.tsv')),
            row.names = T, col.names = F, quote = F, sep = "\t")

# if set condition will evaluate to true and script will terminate# {{{
if(args$factors) {
    print("Print scale factors and exit.")
    # Move file from args$out/deseq to args$out as this probably is generated for
    # normalising the bigwigs
    file.rename(file.path(output_path, paste0(label, '_size-factors.tsv')),
                file.path(args$out, paste0(label, '_size-factors.tsv')))
    # Removes args$out/deseq as it's now empty
    unlink(output_path, recursive = TRUE)
    quit()
}# }}}

# Counts scaled based on library size
norm_counts = counts(cds, normalized = TRUE)
write.table(norm_counts,
            file = file.path(output_path, paste0(label, '.normalised-counts.tsv')),
            row.names = T, col.names = T, quote = F, sep = "\t")

pdf(file.path(plot_path, paste0(label, '.pdf')), paper='a4')
# Running DE and plotting# {{{
dds = DESeq(cds)
rld = rlog(dds)
source("~/source/Rscripts/deseq-functions.R")
# calls plot_counts from deseq-functions to plot counts distributions,
# pca and correlation heatmap for the contrast selected
plot_counts(counts = norm_counts, design = design, type = 'Normalised by SF')
plot_counts(counts = norm_counts, design = design, type = 'Rlog', rld = rld)

res = results(dds)
plotDispEsts(dds)
plotMA(res, alpha = 0.05)

# http://rpackages.ianhowson.com/bioc/DESeq/man/estimateDispersions.html

# }}}

res = res %>% as.data.frame() %>% add_rownames('Gene') %>% arrange(padj, log2FoldChange)
res %>% select(Gene,
               `Base mean` = baseMean,
               `log2 Fold change` = log2FoldChange,
               `Adjusted p` = padj) %>% head() %>% print()

write.table(res, file.path(output_path, paste0(label, '.deseq.tsv')),
            quote = FALSE, row.names = FALSE, sep = "\t")

significant = res %>% filter(!is.na(padj)) %>% filter(padj < 0.05)
significant = significant %>% mutate(Condition = label)

significant_fc = res %>% filter(padj < 0.05, log2FoldChange >= 2.0 )

write.table(significant, file.path(output_path, paste0(label, '.deseq-de.tsv')),
            quote = FALSE, row.names = FALSE, sep = "\t")


# Check if gene length file exists otherwise create it# {{{
if (! file.exists(paste0(file_path_sans_ext(args$gtf), '.gene_lengths.tsv'))) {
    library(tools)
    # Retrieve annotation from gtf files# {{{
    exons_file = paste0(file_path_sans_ext(args$gtf), '.rtracklayer-exons.gtf')
    
    library(rtracklayer)
    if (file.exists(exons_file)){
        print('Importing file.')
        exons = import.gff3(exons_file)
    } else {
        command = paste("Rscript ~/source/get-annotation-rscript.R -G", args$gtf, '-p <ncores>', collapse = " ")
        stop(paste('Annotation GTFs do not exist. Run `', command, '` first.', collapse = ' '))
    }# }}}

    exons = split(exons, exons$ensembl_id)
    lengths = as.data.frame(sum(width(reduce(exons))))
    colnames(lengths) = 'lengths'
    write.table(add_rownames(lengths, var = 'gene_id'),
                paste0(file_path_sans_ext(args$gtf), '.gene_lengths.tsv'), col.names = T, row.names = F, quote = F, sep = "\t")
} else {
    lengths = read.delim(paste0(file_path_sans_ext(args$gtf), '.gene_lengths.tsv'), head = T, sep = "\t", row.names = 1)
}
# }}}

# For TPMs and FPKMs use the raw counts and not the library size or effective library size normalised counts
# a few posts to consider are https://support.bioconductor.org/p/65683/ and https://www.biostars.org/p/99310/
# but the rational behind non scaled is a scaled library x will give different expression values if scaled against
# library y or libraries y and z
efflengths = compute_effective_length(lengths, args$fragment)

# Ensure that genes in count object are as many and in the same order as in efflengths
# number of genes in efflengths may not match the number of genes in counts even if correct annotation is used
# this is due to filtering out genes on contigs when creating the exon file with get-annotation-rscript.R
if(nrow(counts) != nrow(efflengths)) warning(paste0("Number of genes in counts table does not match number of genes in annotation (counts:",
                                                    nrow(counts), ", annotation:", nrow(efflengths), ", missing:", nrow(counts) - nrow(efflengths), ").\n",
                                                    "Assuming annotation is correct and missing genes are on contigs."))

cnt4expression = counts[rownames(efflengths),]
# Making sure the order of genes in counts table and efflengths is the same
efflengths = efflengths[match(rownames(efflengths), rownames(cnt4expression)), ]

# Nuno's functions countstable2tpm and countstable2rpkm
# give same results as mine https://github.com/nunofonseca/irap/blob/master/aux/R/irap_utils.R
fpkms = round(counts_to_FPKM(cnt4expression, setNames(efflengths[['efflength']], rownames(efflengths))), 2)
fpkms = fpkms %>% as.data.frame() %>% add_rownames('Gene')
write.table(fpkms, file.path(output_path, paste0(label, '.raw-count-fpkms.tsv')),
            quote = FALSE, row.names = FALSE, sep = "\t")

tpms = round(counts_to_TPM(cnt4expression, setNames(efflengths[['efflength']], rownames(efflengths))), 2)
tpms = tpms %>% as.data.frame() %>% add_rownames('Gene')
write.table(tpms, file.path(output_path, paste0(label, '.raw-count-tpms.tsv')),
            quote = FALSE, row.names = FALSE, sep = "\t")
# For plotting FPKMs or TPMs use log2

# Check if ensembl to gene name file exists otherwise print command to create it and exit# {{{
if (! file.exists(paste0(file_path_sans_ext(args$gtf), '.ensembl2gene_name.tsv'))) {
    command = paste("Rscript ~/source/get-annotation-rscript.R -G", args$gtf, '-p <ncores>', collapse = " ")
    warning(paste('Annotation GTFs do not exist. Run `', command, '` first.',
                  'Will now exit script without plotting markers', collapse = ' '))
} else {
    gene_names = read.delim(paste0(file_path_sans_ext(args$gtf), '.ensembl2gene_name.tsv'), head = T, sep = "\t")
}
# }}}

markers = read.table("/nfs/research2/bertone/user/mxenoph/hendrich/markers.txt", header = T, stringsAsFactors = F)
tmp = apply(markers, 2, function(gene) gene_names %>% filter(external_gene_id %in% gene))
# http://stackoverflow.com/questions/28250948/how-to-dplyrinner-join-multi-tbls-or-data-frames-in-r
markers = plyr::join_all(lapply(names(tmp), function(y) tmp[[y]] %>% mutate(state = y)), type = "full", by = "ensembl_gene_id")

first_contrast = setNames(as.character(design[['Contrast']]), design[['Library']])
second_contrast = setNames(gsub(".*_", '', as.character(design[['Library']])), design[['Library']])
plot_fpkm(fpkms = tpms, markers = markers, first_contrast = first_contrast, results = res, second_contrast = second_contrast)

# if this is logical means is set to default which is FALSE
if(!is.logical(args$genome_size)){# {{{
    library(RNASeqPower)
    genome_size = read.delim(args$genome_size, sep = '\t', header = F, row.names = 1)
    depth_of_coverage = apply(counts,2,sum)*100/sum(as.numeric(genome_size[[1]]))
    # power is the fraction of true positives that will be detected
    dataset_power = rnapower(depth = round(median(depth_of_coverage)),
                             # cv is 0.1 for inbred animal lines according to ?rnapower
                             n = 2, cv = 0.1, effect = c(1.25, 1.5, 1.75, 2), alpha = c(0.05, 0.1))
    colnames(dataset_power) = paste0('significance_', colnames(dataset_power))
    dataset_power = add_rownames(as.data.frame(dataset_power), var = 'FC')
    write.table(dataset_power, file.path(output_path, paste0(label, '.DEG-detection-power.tsv')),
                quote = FALSE, row.names = FALSE, sep = "\t")
    
}# }}}

if(!is.logical(args$targets)){# {{{
    targets = read.delim(args$targets, sep='\t')
    subsets =  targets %>% filter(!is.na(Bound)) %>% select(Gene) %>% .[[1]] %>% as.character()
    if (!is.logical(args$tname)) {
        subset_name = args$tname
    } else {
        subset_name = 'In subset'
    }
    plot_density(deseq_res = res, subsets = subsets, subset_name = subset_name, label = label)
} else {
    plot_density(deseq_res = res, label = label)
}# }}}

gene_scores = as.data.frame(res) %>% dplyr::select(Gene, padj)
gene_scores = setNames(gene_scores[[2]], gene_scores[[1]])

# Doing GSEA on de genes and de and db genes
source("~/source/Rscripts/functions.R")
get_set_enrichment(gene_scores = gene_scores, label = label)

dev.off()

inprogress = function(){# {{{
library("biomaRt")
db2query = useMart("ensembl")

mart = useDataset("mmusculus_gene_ensembl", db2query)
genes = getBM(filters="ensembl_gene_id",
              attributes = c("ensembl_gene_id", "entrezgene"),
              values = res %>% dplyr::select(Gene),
              mart = mart)

# For all genes in res (i.e. with !is.na(padj) ) report 1 if de or 0 if not
genes = genes %>% mutate(de = ifelse(ensembl_gene_id %in% (significant %>% dplyr::select(Gene)), 1, 0) )

# keep unique defined entrez ids
genes_set = genes %>% filter(!is.na(entrezgene)) %>% unique() %>% dplyr::select(entrezgene, de)
goseq_set = genes_set %>% .[['de']]
names(goseq_set) = genes_set %>% .[['entrezgene']]

library("org.Mm.eg.db")
library("goseq")

pwf = nullp(goseq_set, "mm10","refGene")
# Perform the enrichment analysis
#GO.wall=goseq(pwf,"dm3","refGene")

# Have a look at the first few results
#print(head(GO.wall))
# Store the enriched GOs
#enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<0.1] # 

mapped_genes = mappedkeys(org.Mm.egPATH)
kegg = as.list(org.Mm.egPATH[mapped_genes])
KEGG = goseq(pwf, gene2cat = kegg)
library(KEGGREST)
enriched_pathways = KEGG[KEGG$over_represented_pvalue < 0.05, ]
enriched_pathways$category = paste("dme", enriched_pathways$category, sep="")
sapply(enriched_pathways$category, function(id) {
       return(keggGet(id)[[1]]$NAME) })
}# }}}









