get_annotation <- function (assembly){# {{{
    x <- c('testit', 'rtracklayer', 'GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    assembly <- tolower(assembly)
    database <- read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/assemblies-annotations.config", header=T, stringsAsFactors=F, sep="\t")

    chr <- database[database$assembly == assembly, 'chrom_sizes']
    if(file.exists(chr)){# {{{
        chr <- read.table(chr, row.names = 1)
        assign('chr', chr, envir = globalenv())
    } else {
        assert(paste0(chr, " file does not exist."), !file.exists(chr))
    }# }}}

    ensembl <- database[database$assembly == assembly, 'annotation']

    if(file.exists(ensembl)){# {{{
#        ensembl <- read.table(ensembl, row.names = 1)
        ensembl <- import(ensembl, asRangedData = FALSE)
        assign('ensembl', ensembl, envir = globalenv())
    } else {
        assert(paste0(ensembl, " file does not exist."), !file.exists(ensembl))
    }# }}}
    
    genes <- database[database$assembly == assembly, 'only_genes']

    if(file.exists(genes)){# {{{
        genes <- import(genes, asRangedData = FALSE)
        assign('genes', genes, envir = globalenv())
    } else {
        assert(genes, " file does not exist", !file.exists(genes))
    }# }}}

    tss_window <- database[database$assembly == assembly, 'tss_window']
    if(file.exists(tss_window)){# {{{
        tss_window <- import(tss_window, asRangedData = FALSE)
        assign('tss_window', tss_window, envir = globalenv())

        tss <- GRanges(seqnames(tss_window),
                       IRanges(start = start(tss_window) + (width(tss_window)/2) -1,
                               width = 1),
                       strand = strand(tss_window))
        assign('tss', tss, envir = globalenv())
    } else {
        assert(tss_window, " file does not exist", !file.exists(tss_window))
    }# }}}

}# }}}

get_features <- function(host, dataset, assembly){# {{{
    x <- c('GenomicRanges', 'GenomicFeatures', 'biomaRt', 'rtracklayer')
    lapply(x, suppressMessages(library), character.only=T)
    library('parallel')

    assembly <- tolower(assembly)# {{{

    database <- read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/assemblies-annotations.config",
                           header=T, stringsAsFactors=F, sep="\t")
    chr <- database[database$assembly == assembly, 'chrom_sizes']# }}}
    genome_path <- dirname(chr)

    ensembl<-useMart(host=host, biomart='ENSEMBL_MART_ENSEMBL', dataset=dataset)

    exons <- getBM(attributes=c("ensembl_gene_id",# {{{
                                   "chromosome_name",
                                   "exon_chrom_start",
                                   "exon_chrom_end",
                                   "strand",
                                   "ensembl_exon_id"),
                   mart=ensembl)
    exons <- split(exons, exons$ensembl_gene_id)
    exons <- mclapply(exons, function(x){
                      gr <- with(x, GRanges(seqnames=chromosome_name,
                                            IRanges(exon_chrom_start, exon_chrom_end),
                                            strand))
                      gr <- reduce(gr)
                      gr$exon_no <- paste0('exon_', 1:length(gr))
                      
                      gr <- change_seqnames(gr)
                      gr <- add_seqlengths_to_gr(gr, chr)
                      return(gr)
                   })
    # }}}

    genes <- getBM(attributes=c("ensembl_gene_id",# {{{
                                   "chromosome_name",
                                   "start_position",
                                   "end_position",
                                   "strand"),
                   mart=ensembl)
    genes <- with(genes, GRanges(seqnames = chromosome_name,
                                 IRanges(start_position, end_position),
                                 strand,
                                 names = ensembl_gene_id))
    names(genes) <- genes$names
    genes <- change_seqnames(genes)
    genes <- add_seqlengths_to_gr(genes, chr)
    # }}}

    introns <- mclapply(names(exons), function(x){# {{{
                        setdiff(genes[x], exons[[x]])
                   })
    names(introns) <- names(exons)# }}}

#    protein_coding_db <- makeTranscriptDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",# {{{
#                                                     dataset=dataset,
#                                                     transcript_ids=NULL,
#                                                     circ_seqs=DEFAULT_CIRC_SEQS,
#                                                     filters=as.list(c(biotype="protein_coding")),
#                                                     id_prefix="ensembl_",
#                                                     host=host,
#                                                     miRBaseBuild=NA)# }}}

    genic <- genes# {{{
    strand(genic) <- '*'
    genic <- reduce(genic)

    chromosomes <- as(seqinfo(genic),'GRanges')
    intergenic <- setdiff(chromosomes, genic)# }}}

    # For saving I can't use list so convert to granges object
    exons <- unlist(GRangesList(exons))
    exons$ensembl_id <- names(exons)
    names(exons) <- NULL
    
    introns <- unlist(GRangesList(introns))
    introns$ensembl_id <- names(introns)
    names(introns) <- NULL

    export.gff3(exons, file.path(genome_path, paste0(assembly, '-exons')))
    export.gff3(introns, file.path(genome_path, paste0(assembly, '-introns')))
    export.gff3(intergenic, file.path(genome_path, paste0(assembly, '-intergenic')))

    return(list(exons=exons, introns=introns, genes=genes, genic=genic, intergenic=intergenic))
}# }}}

#feat would be 'exon'# {{{
#match.annot <- function(grange, grfeat, feat){
annotate_range <- function(gr, genomic_feature, label){
    x <- c('GenomicRanges', 'foreach')
    lapply(x, suppressMessages(library), character.only=T)
    library('parallel')
    
    ov <- findOverlaps(gr, genomic_feature)
    # Calculate the width for each overlap
    ov_width <- ranges(ov, ranges(gr), ranges(genomic_feature))
    by_query_range <- splitAsList(ov_width, factor(queryHits(ov)))

    q_attributes <- mclapply(names(by_query_range), function(x, hits, widths, q, s){
                             # For every hit the query got and for which width is in by_query_range
                             # get the index of the subject hits
                             subject_hits <- subjectHits(hits[queryHits(hits) == x])
                             subject_names <- names(s[subject_hits])
                             overlap_width <- width(widths[[x]])
                             query_ov_perc <- overlap_width / width(q[as.numeric(x)])
                             paste(paste0(paste(subject_names, overlap_width, sep=":"),
                                          "(", query_ov_perc, ")"),
                                   collapse="; ")
                   }, ov, by_query_range, gr, genomic_feature)
    names(q_attributes) <- names(by_query_range)
    
    gr <- unlist(GRangesList(mclapply(names(q_attributes), function(x) {
                                      mcols(gr)[as.numeric(x), 'attributes'] <- q_attributes[[x]]
                                      gr[as.numeric(x)]
                                          })
    ))

    return(gr)
}
# }}}


