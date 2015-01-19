
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

get_transcripts <- function(assembly){# {{{
    x <- c('GenomicRanges', 'GenomicFeatures', 'biomaRt', 'rtracklayer')
    lapply(x, suppressMessages(library), character.only=T)

    host <- 'may2012.archive.ensembl.org'
    dataset <- 'mmusculus_gene_ensembl'

    database <- useMart(host=host, biomart='ENSEMBL_MART_ENSEMBL', dataset=dataset)
    genes <- getBM(attributes=c("chromosome_name", "start_position", "end_position", "strand", "source",
                                "gene_biotype", "ensembl_gene_id", "entrezgene", "mgi_symbol"),
                   mart=ensembl)
}
# }}}

get_features <- function(host, ds, chr.sizes){# {{{
    x <- c('GenomicRanges', 'GenomicFeatures', 'biomaRt', 'rtracklayer')
    lapply(x, suppressMessages(library), character.only=T)

    #ensembl<-useMart(host=host, biomart='ENSEMBL_MART_ENSEMBL', dataset=ds)
    #annotation <- getBM(attributes=c("ensembl_gene_id", "chromosome_name", "strand", "exon_chrom_start", "exon_chrom_end", "ensembl_exon_id"),  mart=ensembl)
    #keep only chromosome
    #annotation <- annotation[grep("^[\\dXYMT]{0,2}$", df$chromosome_name, perl=TRUE),]

    txdb <- makeTranscriptDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset=ds, transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filters="", id_prefix="ensembl_", host=host, miRBaseBuild=NA)

    exons <- exonsBy(txdb, by=c("gene"))
    exons <- unlist(exons)
    exons <- alter.seqnames4granges(exons)

    introns <- intronsByTranscript(txdb, use.names=TRUE)

    # Unlist and remove duplicate introns.
    introns  <- unlist(introns)
    intronsNoDups <- introns[!duplicated(introns)]

    #The select() function can map txid to geneid (and others). See ?select.
    txnames <- names(intronsNoDups)
    map <- select(txdb, keys=txnames, keytype='TXNAME', columns='GENEID')
    map<- unique(map)
    names(intronsNoDups) <- map[match(names(intronsNoDups), map$TXNAME), 'GENEID']

    #This not working as expected
    #idx <- map$GENEID[!is.na(map$GENEID)]
    #intronsByGene <- split(intronsNoDups[!is.na(map$GENEID)], idx)
    #names(intronsByGene) <- unique(idx)
    intronsByGene <- split(intronsNoDups, names(intronsNoDups))
    intronsByGene <- unlist(intronsByGene, use.names=FALSE)
    punion(intronsByGene, intronsByGene)
    intronsByGene <- alter.seqnames4granges(intronsByGene)

    #CAUTION:The whole exon thing this way is a bit faulty-- correct before using

    #Getting the transcripts for all genes
    trx <- transcriptsBy(txdb, by='gene')
    #Combining transcripts ranges for each gene to get genes
    genes <- unlist(range(trx))
    genes <- alter.seqnames4granges(genes)
    #this can be called for mm9 so better not add the mm10 chrom sizes
    genes <- add.seqlengths(genes, chr.sizes)
    genes.nostrand <- genes
    strand(genes.nostrand) <- '*'

    genic <- reduce(genes.nostrand)
    chrom.gr<-as(seqinfo(genic),'GRanges')
    intergenic<-setdiff(chrom.gr,genic)

    genic.str <- reduce(genes)
    chrom.gr.str <- c(chrom.gr, chrom.gr)
    strand(chrom.gr.str[1:length(chrom.gr)])<-'+'
    strand(chrom.gr.str[(length(chrom.gr)+1):length(chrom.gr.str)])<-'-'
    intergenic.str <- setdiff(chrom.gr.str, genic.str)
    return(list(exons=exons, genes=genes, genic=genic, intergenic=intergenic, intergenic.str=intergenic.str, introns=intronsByGene))

}# }}}
