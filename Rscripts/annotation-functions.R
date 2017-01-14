library = function (...) suppressMessages(base::library(...))

get_annotation <- function (assembly){# {{{
    x <- c('rtracklayer', 'GenomicRanges')
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
        names(genes) <- names(genes) <- genes$gene_id
        assign('genes', genes, envir = globalenv())
    } else {
        assert(genes, " file does not exist", !file.exists(genes))
    }# }}}

    tss_window <- database[database$assembly == assembly, 'tss_window']
    if(file.exists(tss_window)){# {{{
        tss_window <- import(tss_window, asRangedData = FALSE)
        names(tss_window) <- tss_window$gene_id
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

get_features = function(archive = "ensembl", organism = "Mus_musculus",
                        database = "/nfs/research2/bertone/user/mxenoph/common/genome/assemblies-annotations.conf",
                        genome_path = "/nfs/research2/bertone/user/mxenoph/common/genome"){ # {{{
    x = c('GenomicRanges', 'GenomicFeatures', 'biomaRt', 'rtracklayer', 'parallel')
    lapply(x, suppressMessages(library), character.only=T)

    get_eversion = function(html){# {{{
        library(XML)
        html = htmlTreeParse(html, useInternal = TRUE)
        version = unlist(xpathApply(html, '//*/title', xmlValue))
        # Build pattern and extract numeric substring for ensembl version
        pattern = gregexpr("[0-9]+", version)
        version = unlist(regmatches(version, pattern))
        assembly = unlist(xpathApply(html, '//h2', xmlValue))
        assembly = assembly[grep("Genome assembly", assembly)]
        assembly = unlist(strsplit(gsub("Genome assembly: ", "", assembly), " "))[1]
        # Removing .p from genome assembly
        assembly = gsub("\\..*", "", assembly)
        return(list(assembly = assembly, version = version))
    }# }}}

    filename = organism

    # Converting e.g. Mus_musculus to mmusculus
    organism = unlist(strsplit(organism, '_'))
    organism[1] = substr(organism[1], 1, 1)
    organism = tolower(paste(organism, collapse=''))

    # e70 is jan2013 and e67 is may2012# {{{
    if (!grepl("ensembl", archive)) {
        # Retrieve current assembly and version of ensembl
        description = unlist(get_eversion(paste0("http://", archive, ".archive.ensembl.org/", filename, "/Info/Index")))
        print(paste0("Retrieving archived version of ensembl: ", paste(description, collapse = '.')))
        archive = paste0(archive, ".archive.ensembl.org")
        ensembl = useMart(host = archive, biomart = 'ENSEMBL_MART_ENSEMBL', dataset = paste0(organism, "_gene_ensembl"))
    } else {
        # Retrieve current assembly and version of ensembl
        description = unlist(get_eversion(paste0("http://www.ensembl.org/", filename, "/Info/Index")))
        print(paste0("Retrieving current version of ensembl:", paste(description, collapse = '.')))
        ensembl = useMart(biomart = 'ensembl', dataset = paste0(organism, "_gene_ensembl"))
    }# }}}

    filename = paste(filename, paste(description, collapse = '.'), sep = '.')
    database = read.table(database, header = T, stringsAsFactors = F, sep = "\t")
    chromosomes = database[database$GRC == description['assembly'], 'chrom_sizes']
    genome_path = file.path(genome_path,
                            toupper(database[database$GRC == description['assembly'], 'assembly']))

    transcripts = getBM(attributes = c("ensembl_gene_id",# {{{
                                       "chromosome_name",
                                       "transcript_start",
                                       "strand",
                                       "transcript_end",
                                       "transcript_biotype",
                                       "transcript_count",
                                       "ensembl_transcript_id"),
                        mart = ensembl)
    transcripts = split(transcripts, transcripts$ensembl_gene_id)
    transcripts = mclapply(transcripts, function(x){
                           gr = with(x, GRanges(seqnames = chromosome_name,
                                                IRanges(transcript_start, transcript_end),
                                                strand))
                           values(gr) = x %>% select(-chromosome_name, -transcript_start, -transcript_end, -strand)
                           gr = change_seqnames(gr)
                           gr = add_seqlengths_to_gr(gr, chromosomes)
                           return(gr)
                        })
    transcripts = unlist(GRangesList(transcripts))# }}}

    exons = getBM(attributes = c("ensembl_gene_id",# {{{
                                 "chromosome_name",
                                 "exon_chrom_start",
                                 "exon_chrom_end",
                                 "strand",
                                 "ensembl_exon_id"),
                  mart=ensembl)
    exons = split(exons, exons$ensembl_gene_id)
    exons = mclapply(exons, function(x){
                      gr = with(x, GRanges(seqnames=chromosome_name,
                                            IRanges(exon_chrom_start, exon_chrom_end),
                                            strand))
                      gr = reduce(gr)
                      gr$exon_no = paste0('exon_', 1:length(gr))
                      gr = change_seqnames(gr)
                      gr = add_seqlengths_to_gr(gr, chr)
                      return(gr)
                   })
    # }}}

    genes = getBM(attributes = c("ensembl_gene_id",# {{{
                                 "chromosome_name",
                                 "start_position",
                                 "end_position",
                                 "strand"),
                  mart=ensembl)
    genes = with(genes, GRanges(seqnames = chromosome_name,
                                 IRanges(start_position, end_position),
                                 strand,
                                 names = ensembl_gene_id))
    names(genes) = genes$names
    genes = change_seqnames(genes)
    genes = add_seqlengths_to_gr(genes, chr)
    # }}}

    introns = mclapply(names(exons), function(x){ # {{{
                       BiocGenerics::setdiff(genes[x], exons[[x]])
                   })
    names(introns) = names(exons)
    # }}}

    genic = genes # {{{
    strand(genic) = '*'
    genic = reduce(genic)

    chromosomes = as(seqinfo(genic),'GRanges')
#    intergenic = setdiff(chromosomes, genic)
    intergenic = gaps(genic)
    # }}}

    # For saving I can't use list so convert to granges object
    exons = unlist(GRangesList(exons))
    exons$ensembl_id = names(exons)
    names(exons) = NULL
    
    introns = unlist(GRangesList(introns))
    introns$ensembl_id = names(introns)
    names(introns) = NULL

    export.gff3(genes, file.path(genome_path, paste0(assembly, '-genes')))
    export.gff3(exons, file.path(genome_path, paste0(assembly, '-exons')))
    export.gff3(introns, file.path(genome_path, paste0(assembly, '-introns')))
    export.gff3(intergenic, file.path(genome_path, paste0(assembly, '-intergenic')))

    return(list(exons = exons, introns = introns, genes = genes, genic = genic, intergenic = intergenic))
}# }}}

annotate_range = function(gr, genomic_feature, label, return_att = F){ # {{{
    x = c('GenomicRanges', 'foreach')
    lapply(x, suppressMessages(library), character.only=T)
    library('parallel')
    
    ov = findOverlaps(gr, genomic_feature)
    # Calculate the width for each overlap
    ov_width = ranges(ov, ranges(gr), ranges(genomic_feature))
    by_query_range = splitAsList(ov_width, factor(queryHits(ov)))

    q_attributes = mclapply(names(by_query_range), function(x, hits, widths, q, s){
                             # For every hit the query got and for which width is in by_query_range
                             # get the index of the subject hits
                             subject_hits = subjectHits(hits[queryHits(hits) == x])
                             subject_names = names(s[subject_hits])
                             overlap_width = width(widths[[x]])
                             query_ov_perc = overlap_width / width(q[as.numeric(x)])
                             paste(paste0(paste(subject_names, overlap_width, sep=":"),
                                          "(", query_ov_perc, ")"),
                                   collapse="; ")
                   }, ov, by_query_range, gr, genomic_feature)
    names(q_attributes) = names(by_query_range)
    
    mcols(gr)[as.numeric(names(q_attributes)), 'attributes'] = unlist(q_attributes)

    # If true it will return a list of identifiers for the features that overlap
    # the input grange
    if (return_att == T) {
        attribute = unlist(lapply(gr$attributes[!is.na(gr$attributes)],
                                   function(x) unlist(lapply(unlist(strsplit(x, '; ')),
                                                             function(y) unlist(strsplit(y, ':'))[1]))))
        return(list(gr = gr, attribute = unique(attribute)))
    } else {
        return(gr)
    }
}
# }}}

annotate_regions = function(grange, annotations,
                            priority = c('2000-gene_start-500', 'first-exon', 'active_enhancers',
                                         'poised_enhancers', 'exons', 'first-intron', 'introns', 'genes')) { #{{{
    x = c('GenomicRanges', 'foreach')
    sapply(x, library, character.only=T)
    library('parallel')

    stopifnot(is(grange, "GRanges"))
    for (i in priority[1]) {

        ov = findOverlaps(grange, annotations[[i]])
        overlapping_ranges = ranges(ov, ranges(grange), ranges(annotations[[i]]))
        by_query = splitAsList(overlapping_ranges, factor(queryHits(ov)))
    }
    # Calculate the width for each overlap
    ov_width = ranges(ov, ranges(gr), ranges(genomic_feature))
    by_query_range = splitAsList(ov_width, factor(queryHits(ov)))

    q_attributes = mclapply(names(by_query_range), function(x, hits, widths, q, s){
                             # For every hit the query got and for which width is in by_query_range
                             # get the index of the subject hits
                             subject_hits = subjectHits(hits[queryHits(hits) == x])
                             subject_names = names(s[subject_hits])
                             overlap_width = width(widths[[x]])
                             query_ov_perc = overlap_width / width(q[as.numeric(x)])
                             paste(paste0(paste(subject_names, overlap_width, sep=":"),
                                          "(", query_ov_perc, ")"),
                                   collapse="; ")
                   }, ov, by_query_range, gr, genomic_feature)
    names(q_attributes) = names(by_query_range)
    
    mcols(gr)[as.numeric(names(q_attributes)), 'attributes'] = unlist(q_attributes)

    # If true it will return a list of identifiers for the features that overlap
    # the input grange
    if (return_att == T) {
        attribute = unlist(lapply(gr$attributes[!is.na(gr$attributes)],
                                   function(x) unlist(lapply(unlist(strsplit(x, '; ')),
                                                             function(y) unlist(strsplit(y, ':'))[1]))))
        return(list(gr = gr, attribute = unique(attribute)))
    } else {
        return(gr)
    }
}
