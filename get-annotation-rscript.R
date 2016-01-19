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

args = parser$parse_args()
ncores = args$threads
# }}}

get_features = function(archive = NULL, organism = NULL, output = NULL){ # {{{
    print('In get_features')
    # Load packages# {{{
    library(GenomicFeatures)
    library(GenomicRanges)
    library(rtracklayer)
    library(parallel)
    library(BiocParallel)
    # }}}

    # formatting dataframes into Granges# {{{
    format_features = function(features, reduce = F, label = NULL){
        tmp = features %>% dplyr::select(-chromosome_name, -strand, -contains("start"), -contains("end"))
        features = with(features, GRanges(seqnames = chromosome_name,
                                          IRanges(features %>% dplyr::select(contains("start")) %>% .[[1]],
                                                  features %>% dplyr::select(contains("end")) %>% .[[1]]),
                                          strand))
        values(features) = tmp
        features = change_seqnames(features)
        features = add_seqlengths_to_gr(features, ucsc_chromosomes)

        features = split(features, features$ensembl_gene_id)
        if(reduce)  features = reduce(features)
        if(!is.null(label)) {
            features = mclapply(features, function(x) {
                                x = sort(x)
                                if(runValue(strand(x)) == '+') {
                                    values(x)[[paste0(label, "_no")]] = paste0(label, "_", 1:length(x))
                                }else{
                                    values(x)[[paste0(label, "_no")]] = paste0(label, "_", length(x):1)
                                }
                                return(x)
                       }, mc.cores = ncores)
        }
        features = unlist(features)
        return(features)
    }
    # }}}

    # Getting chromosome lengths# {{{
    chrom_info = getChromInfoFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                                         dataset = paste0(organism, "_gene_ensembl"), 
                                         id_prefix="ensembl_",
                                         host=archive,
                                         port=80
                                         )
    # As assembly names from Genome Genome reference Consortium (used by ensembl) and NCBI are different
    # I can not use the information in `assemblies` to retrieve the chromosome UCSC chromosome names with 
    # getChromInfoFromUCSC. Assembly names could be matched by parsing https://genome.ucsc.edu/FAQ/FAQreleases.html
    # but that might cause the script to crash with future updates.
    # The inconsistency between chromosome names is solved only for autosomal, sex and mitochondria chromosomes here.
    # Hence genes/transcipts on contigs will be removed further down
    ucsc_chromosomes = chrom_info %>%
                        mutate(chrom = ifelse(chrom %in% c(1:20, 'X', 'Y'), paste0('chr', chrom),
                                              ifelse(chrom == 'MT', 'chrM', as.character(chrom))))
# }}}

    # Transcript biotypes# {{{
    # Based on http://www.ensembl.org/Help/Faq?id=468
    # nontranslating CDS missing from protein_coding as not sure how it's formatted in biomaRt
    # In any case for mouse it is anyway not included in current annotation
    ebiotypes = data.frame(transcript_biotype = c(c('protein_coding', 'IG_gene', 'IG_C_gene', 'IG_D_gene',
                                                    'IG_J_gene', 'IG_LV_gene', 'IG_M_gene', 'IG_V_gene', 
                                                    'IG_Z_gene', 'nonsense_mediated_decay', 'non_stop_decay',
                                                    'polymorphic_pseudogene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene'),
                                                  c('disrupted_domain', 'IG_C_pseudogene', 'IG_J_pseudogene', 'IG_pseudogene',
                                                    'IG_V_pseudogene', 'processed_pseudogene', 'transcribed_processed_pseudogene',
                                                    'unitary_pseudogene', 'unprocessed_pseudogene', 'translated_processed_pseudogene',
                                                    'TR_J_pseudogene', 'pseudogene', 'TR_V_pseudogene'),
                                                  c('3prime_overlapping_ncrna', 'ambiguous_orf', 'antisense', 'antisense_rna',
                                                    'lincRNA', 'ncrna_host', 'processed_transcript', 'sense_intronic',
                                                    'sense_overlapping', 'non_coding', 'retained_intron'),
                                                  c('miRNA', 'miRNA_pseudogene', 'misc_RNA', 'misc_RNA_pseudogene', 'Mt_rRNA',
                                                    'Mt_tRNA', 'rRNA', 'scRNA', 'snlRNA', 'snoRNA', 'snRNA', 'tRNA', 'tRNA_pseudogene'),
                                                  'TEC'),
                            Class = c(rep('protein_coding', 16), rep('pseudogene', 13), rep('long_noncoding', 11),
                                        rep('short_noncoding', 13), 'TEC')
                     )
    # }}}

    ensembl = useMart(host = archive,
                      biomart = 'ENSEMBL_MART_ENSEMBL',
                      dataset = paste0(organism, "_gene_ensembl"))

    ensembl2gene_name = getBM(attributes = c("ensembl_gene_id", "external_gene_id"), mart=ensembl)
    write.table(ensembl2gene_name, paste0(output, ".ensembl2gene_name", ".tsv"), sep="\t", quote=FALSE, col.names=T, row.names=F)
    
    protein_coding_table = getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'transcript_biotype'),
                                 filters ='biotype',
                                 values = ebiotypes[ebiotypes$Class == 'protein_coding', 'transcript_biotype'],
                                 mart = ensembl)
    write.table(protein_coding_table,
               paste0(output, ".protein_coding_ids", ".tsv"), sep="\t", quote=FALSE, col.names=T, row.names=F)

    # One can retrieve transcripts start and ends with biomaRt but will have to retrieve exon length
    # for exons in each transcript to calculate transcript length. transcriptLenghts from GenomicFeatures 
    # correctly handle this. Transcript length is needed to find the canonical transcript
    # transcript databases. Only using txdb to get transcipt length# {{{
    txdb = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL",
                        dataset = paste0(organism, "_gene_ensembl"),
                        transcript_ids = NULL,
                        circ_seqs = DEFAULT_CIRC_SEQS,
                        id_prefix = "ensembl_",
                        host = archive,
                        port = 80,
                        miRBaseBuild = NA)

    transcripts = getBM(attributes = c("ensembl_gene_id",
                                       "chromosome_name",
                                       "transcript_start",
                                       "strand",
                                       "transcript_end",
                                       "transcript_biotype",
                                       "transcript_count",
                                       "ensembl_transcript_id"),
                        mart = ensembl)
    # }}}

    # Canonical transcripts# {{{
    transcript_lengths = transcriptLengths(txdb, with.cds_len = FALSE, with.utr5_len = FALSE, with.utr3_len = FALSE)
    if(nrow(transcripts) != nrow(transcript_lengths)) stop("Different number of transcripts retrieved from biomaRt and GenomicFeatures")
    transcript_lengths = transcript_lengths %>% dplyr::rename(ensembl_transcript_id = tx_name)
    # Definition from http://lists.ensembl.org/pipermail/dev/2012-April/007435.html
    canonical_transcripts = inner_join(transcripts, transcript_lengths, by = "ensembl_transcript_id") %>% dplyr::select(-gene_id)
    canonical_transcripts = canonical_transcripts %>% left_join(ebiotypes, by = 'transcript_biotype') %>%
                            mutate(Class = factor(Class, levels = c('protein_coding', 'long_noncoding', 'short_noncoding', 'pseudogene', 'TEC'))) %>%
                            group_by(ensembl_gene_id, Class) %>%
                            top_n(n=1, wt=tx_len) %>% arrange(Class)

    canonical_transcripts = as.data.frame(do.call(rbind,
                                                  lapply(split(canonical_transcripts, canonical_transcripts$ensembl_gene_id),
                                                         function(x) x[1,]))) %>%
                            dplyr::select(ensembl_gene_id, 6:12) %>% filter(Class %in% c('protein_coding', 'long_noncoding'))

    # }}}

    # All transcripts on chromosomes
    transcripts_gr = format_features(transcripts)

    canonical_transcripts_gr = transcripts_gr[transcripts_gr$ensembl_transcript_id %in% canonical_transcripts$ensembl_transcript_id,]
    print('finding canonical trx')
    # Exons # {{{
    exons = getBM(attributes = c("ensembl_gene_id",
                                 "chromosome_name",
                                 "exon_chrom_start",
                                 "exon_chrom_end",
                                 "strand",
                                 "ensembl_exon_id"),
                  mart=ensembl)

    exons_gr = GRangesList(format_features(exons, reduce = T, label = 'exon'))
    # }}}
    
    print('finding exons')
    # Genes # {{{
    genes = getBM(attributes = c("ensembl_gene_id",
                                 "chromosome_name",
                                 "start_position",
                                 "end_position",
                                 "strand"),
                  mart=ensembl)
    genes_gr = format_features(genes)
    # }}}

    print('finding genes')
    # Introns # {{{
    introns_gr = lapply(names(exons_gr), function(x){
                       features = BiocGenerics::setdiff(genes_gr[x], exons_gr[[x]])
                       features= sort(features)
                       if(length(features) != 0) {
                           if(runValue(strand(features)) == '+') {
                               values(features)[["intron_no"]] = paste0("intron_", 1:length(features))
                           }else{
                               values(features)[["intron_no"]] = paste0("intron_", length(features):1)
                           }
                       }
                       return(features)
                  })

    names(introns_gr) = names(exons_gr)
   # }}}

    print('finding introns')

   # Intergenic # {{{
    genic = genes_gr
    strand(genic) = '*'
    genic = reduce(genic)

    chromosomes = as(seqinfo(genic),'GRanges')
    intergenic = GenomicRanges::setdiff(chromosomes, genic)
    # required for proper convertion to bed
    intergenic$id = 1:length(intergenic)
    # }}}
    
    print('finding intergenic')
    promoter_window = data.frame('upstream' = c(5000, 2000), 'downstream' = c(2000, 500))
    for(x in 1:nrow(promoter_window)){
        # Using trim to deal with chrM (circular) which will give negative start()
        # trim sets the range to start and end at feasible coordinates
        promoter = trim(promoters(transcripts_gr,
                                  downstream = promoter_window[x, 'downstream'],
                                  upstream = promoter_window[x, 'upstream']))
        gene_start = trim(promoters(genes_gr,
                                  downstream = promoter_window[x, 'downstream'],
                                  upstream = promoter_window[x, 'upstream']))
        canonical_tx_promoter = trim(promoters(canonical_transcripts_gr,
                                               downstream = promoter_window[x, 'downstream'],
                                               upstream = promoter_window[x, 'upstream']))
        export.gff3(promoter, paste0(output,
                                     print(paste0(".rtracklayer-", promoter_window[x, 'upstream'],
                                                  "-tss-", promoter_window[x, 'downstream'])),
                                     ".gtf"))
        export.gff3(gene_start, paste0(output,
                                     print(paste0(".rtracklayer-", promoter_window[x, 'upstream'],
                                                  "-gene_start-", promoter_window[x, 'downstream'])),
                                     ".gtf"))
        export.gff3(canonical_tx_promoter, paste0(output,
                                     print(paste0(".rtracklayer-", promoter_window[x, 'upstream'],
                                                  "-canonical_tss-", promoter_window[x, 'downstream'])),
                                     ".gtf"))
   }

    # Writing to file
    # For saving I can't use list so convert to granges object
    exons_gr = unlist(exons_gr)
    exons_gr$ensembl_id = names(exons_gr)
    names(exons_gr) = NULL
    export.gff3(exons_gr, paste0(output, ".rtracklayer-exons", ".gtf"))
    export.gff3(exons_gr[which(exons_gr$exon_no == "exon_1")], paste0(output, ".rtracklayer-first-exon", ".gtf"))
    
    names(genes_gr) = NULL
    export.gff3(genes_gr, paste0(output, ".rtracklayer-genes", ".gtf"))

    introns_gr = unlist(GRangesList(introns_gr))
    introns_gr$ensembl_id = names(introns_gr)
    names(introns_gr) = NULL
    export.gff3(introns_gr, paste0(output, ".rtracklayer-introns", ".gtf"))
    export.gff3(introns_gr[which(introns_gr$intron_no == "intron_1")], paste0(output, ".rtracklayer-first-intron", ".gtf"))
   
    export.gff3(intergenic, paste0(output, ".rtracklayer-intergenic", ".gtf"))
    export.gff3(canonical_transcripts_gr, paste0(output, ".rtracklayer-canonical-transcripts", ".gtf"))
}

# }}}

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

    # Might not be needed at the end# {{{
    assemblies = read_html(args$etable) %>% html_nodes('td') %>% html_text()
    html_tables = read_html(args$etable) %>% html_nodes('table') %>% html_table()
    html_tables = lapply(html_tables, function(x){
                         colnames(x) = gsub("^\\W+$", "Organism", colnames(x))
                         x = x[!grepl("^\\W+$", x[,'Organism'], perl=T),]
                         return(x)
                    })

    assemblies = html_tables[[1]]
    for(i in 2:length(html_tables)){
        assemblies = left_join(assemblies, html_tables[[i]], by = "Organism")
    }
    assemblies = apply(assemblies, 1, function(per_row){
                       for(j in 2:length(per_row)){
                           if(is.na(per_row[j])) per_row[j] = per_row[j-1]
                       }
                       return(per_row)
                    })
    options(stringsAsFactors = T)
    assemblies = as.data.frame(t(assemblies))
    colnames(assemblies) = gsub(".+v", "v", colnames(assemblies))
    # }}}

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
    get_features(archive = archive, organism = organism, output = target)

}
# }}}

main()
