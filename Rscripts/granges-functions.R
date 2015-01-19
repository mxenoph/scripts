source("~/source/Rscripts/functions.R")

#Up and working fine! Could allow to pass the chrom.length as a dataframe too so I don't have to read it up every time. will have to check the class basically -> character it's a file otherwise it should be df
#chrom.length should be the file e.g. "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10_noHead.genome"
# chr should be either a file.name or a dataframe
add_seqlengths_to_gr <- function(gr, chr){# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    # if not data.frame it should be a file pointing 
    if(class(chr) != "data.frame"){
        if( !file.exists(chr) ) { 
            stop(chr, " doesn't exist.")
        } else {
            chr <- read.table(chr, sep="\t", row.names=1)
            seqlengths(gr) <- chr[names(seqlengths(gr)),]
        }
    } else {
        name <- which(apply(chr, 2, function(y) any(grepl("chr\\d+", y, perl=T))))
        # if chr ar not a column
        if (length(name) == 0) {
            if( !any(grepl("chr\\d+", rownames(chr), perl=T)) ) stop("Dataframe does not contain chromosome names.")
            seqlengths(gr) <- chr[names(seqlengths(gr)),]
        } else {
            rownames(chr) <- chr[[name]]
            # removing the column with chr names as added the rownames in the above line
            chr[[name]] <- NULL
            seqlengths(gr) <- chr[names(seqlengths(gr)),]
        }
    }

    return(gr)
}# }}}

alter.seqnames4granges<-function(grange){# {{{
    # List of packages to load
    x <- c('GenomicRanges', 'gtools')
    lapply(x, suppressMessages(library), character.only=T)

    chr2num <- function(grange){
      seqlev <- seqlevels(grange)[grepl("^chr[0-9XYM]", seqlevels(grange))]
      chr <- mixedsort(seqlev)
      new <- mixedsort(gsub('chr', '',seqlev))
      names(new) <- paste("chr", new, sep = "")
      new["chrM"] <- "MT"

      seqlevels(grange, force=TRUE) <- chr
      grange <- keepSeqlevels(grange, chr)
      grange <- renameSeqlevels(grange, new)
      return(grange)
    }

    num2chr <- function(grange){
      seqlev <- seqlevels(grange)[grepl("^[0-9XY]|MT", seqlevels(grange))]
      seqlev[match("MT", seqlev)] <- "M"
      chr <- mixedsort(seqlev)
      new <- mixedsort(gsub("^", 'chr', seqlev))
      names(new) <- mixedsort(seqlev)

      seqlevels(grange, force=TRUE) <- chr
      grange <- keepSeqlevels(grange, chr)
      grange <- renameSeqlevels(grange, new)
      return(grange)
    }

    if(grepl('chr', seqlevels(grange)[1])){
      return(chr2num(grange))
    }
    else{
      return(num2chr(grange))
    }

}# }}}

#host for ensembl70 is 'jan2013.archive.ensembl.org', ds= 'mmusculus_gene_ensembl', chr.sizes=e.g. "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10.chrom.sizes"
get.annotation<-function(host, ds, chr.sizes){# {{{
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

#feat would be 'exon'# {{{
match.annot <- function(grange, grfeat, feat){
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)
    
    overlaps <- findOverlaps(grange, grfeat)
    q <- queryHits(overlaps)
    names(q) <- names(grfeat[subjectHits(overlaps)])
    #for when i thought i might add strand info for the exon/intron
    #names(q) <- as.integer(subjectHits(overlaps))

    #This is an IntegerList that i don't know how to coerce so I get the length as string
    #w <- width(pintersect(grange[queryHits(overlaps)], grfeat[subjectHits(overlaps)]))

    #split the list by feature matching index
    tmp.split <- split(q, q)
    tmp.split <- lapply(tmp.split, function(i){paste(unique(names(i)), collapse=';')})

    #when I thought it would be easier to get strand this way
    #tmp.split <- lapply(tmp.split, function(i){paste(unique(names(grfeat[as.integer(names(i))])), collapse=';')})
    tmp.split <- unlist(tmp.split)
    tmp<-rep('NA', length(grange))

    tmp[as.integer(names(tmp.split))] <- unname(tmp.split)
    metadata <- paste('spanning', feat, sep='.')
    values(grange) <- cbind(values(grange), structure(data.frame(tmp), names=metadata))
    return(grange)
}
# }}}

##################
macs_to_granges <- function(peaks){# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    peaks <- read.table(peaks, header=T, sep="\t")
    gr <- with(peaks, GRanges(
                 seqnames= chr,
                 range= IRanges(start=start, end=end),
                 strand= "*",
                 wd= length,
                 score= X.10.log10.pvalue.,
                 FE= fold_enrichment,
                 summit= start + summit,
                 tags= tags))
                 
    if('FDR...' %in% colnames(peaks)) values(gr)[['fdr']] <- with(peaks, FDR...) # fdr column is not available if the experiment had no control
    return(gr)
}
# }}}


macs2_to_granges <- function(peaks){# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    peaks <- read.table(peaks, header=T, sep="\t")
    gr <- with(peaks, GRanges(
                 seqnames= chr,
                 range= IRanges(start=start, end=end),
                 strand= "*",
                 wd= length,
                 tags= pileup,
                 score= X.log10.pvalue.,
                 FE= fold_enrichment))
                 
    if('X.log10.qvalue.' %in% colnames(peaks)) values(gr)[['fdr']] <- with(peaks, 'X.log10.qvalue.') # fdr column is not available if the experiment had no control
    return(gr)
}
# }}}

name_gr <- function(gr) paste(seqnames(gr), start(gr), end(gr), sep=':')
get_score <- function(gr, value) values(gr)[[value]]

intersect_multi <- function(grl){# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    if(length(grl) < 2) return(grl[[1]])
    gr <- grl[[1]]
    for (j in seq(2, length(grl), 1)){
        gr <- intersect(gr, grl[[j]])
    }
    names(gr) <- name_gr(gr)
    return(gr)
}# }}}

get_hits_scores <- function(query, subj, value) {# {{{
    # List of packages to load
    x <- c('GenomicRanges')
    lapply(x, suppressMessages(library), character.only=T)

    # returns a vector of length(common) containing index for hits in x
    ov <- findOverlaps(query, subj, select='first')
    names(ov) <- names(query)
    # subset based on those that have a hit
    ov[!is.na(ov)] <- values(subj[as.numeric(ov[!is.na(ov)])])[[value]]
    return(ov)
}# }}}


