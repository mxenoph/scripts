source("~/local/functions.R")
library(GenomicRanges)
library(GenomicFeatures)
library(biovizBase)
library(biomaRt)
library(ggplot2)
library(ggbio)
library(rtracklayer)


#Up and working fine! Could allow to pass the chrom.length as a dataframe too so I don't have to read it up every time. will have to check the class basically -> character it's a file otherwise it should be df
#chrom.length should be the file e.g. "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10_noHead.genome"
add.seqlengths <- function(grange, chrom.length){

   if( !file.exists(chrom.length) ) { stop(chrom.length, " doesn't exist.")}
   chrom.length <- read.table(chrom.length, sep="\t", row.names=1)
   seqlengths(grange) <- chrom.length[names(seqlengths(grange)),]
   return(grange)
}

alter.seqnames4granges<-function(grange){
   require('gtools')

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

}

#host for ensembl70 is 'jan2013.archive.ensembl.org', ds= 'mmusculus_gene_ensembl', chr.sizes=e.g. "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10.chrom.sizes"
get.annotation<-function(host, ds, chr.sizes){
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

}

#feat would be 'exon'
match.annot <- function(grange, grfeat, feat){
   require(GenomicRanges)
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

##################
macs2GRanges <- function(peaks){
    require(GenomicRanges)
    peaks <- read.table(peaks, header=T, sep="\t")
    gr <- with(peaks, GRanges(
                 seqnames= chr,
                 range= IRanges(start=start, end=end),
                 strand= "*",
                 wd= length,
                 score= X.10.log10.pvalue.,
                 FE= fold_enrichment,
                 fdr= FDR...,
                 summit= start + summit,
                 tags= tags))
   return(gr)
}


