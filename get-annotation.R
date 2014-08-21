#host for ensembl70 is 'jan2013.archive.ensembl.org', ds= 'mmusculus_gene_ensembl'
get.annotation<-function(host, ds, chromL){
   require('GeonomicFeatures')

   txdb <- makeTranscriptDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset=ds, transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filters="", id_prefix="ensembl_", host=host, miRBaseBuild=NA)

   exons <- exonsBy(txdb, by=c("gene"))
   #Elements are already sorted so number the exons so you can get exon one later
   exons <- endoapply(exons, function(x){
                      d <- data.frame('exon_num'= seq(1:length(x)));
                      values(x)<- cbind(values(x),d);
                      return(x)})
   first.exon <- endoapply(exons, function(x){
                           x[x$exon_num==1,]})
   exons <- unlist(exons)
   first.exon <- unlist(first.exons)
   exons <- alter.seqnames4granges(exons)
   first.exon <- alter.seqnames4granges(first.exon)

   introns <- intronsByTranscript(txdb, use.names=TRUE)

   # Unlist and remove duplicate introns.
   introns  <- unlist(introns)
   intronsNoDups <- introns[!duplicated(introns)]

   #The select() function can map txid to geneid (and others). See ?select.
   txnames <- names(intronsNoDups)
   map <- select(txdb, keys=txnames, keytype='TXNAME', cols='GENEID')
   map<- unique(map)
   names(intronsNoDups) <- map[match(names(intronsNoDups), map$TXNAME), 'GENEID']

   #This not working as expected
   #idx <- map$GENEID[!is.na(map$GENEID)]
   #intronsByGene <- split(intronsNoDups[!is.na(map$GENEID)], idx)
   #names(intronsByGene) <- unique(idx)
   intronsByGene <- split(intronsNoDups, names(intronsNoDups))
   intronsByGene <- unlist(intronsByGene, use.names=FALSE)
   union(intronsByGene)
   intronsByGene <- alter.seqnames4granges(intronsByGene)

   #Getting the transcripts for all genes
   trx <- transcriptsBy(txdb, by='gene')
   #Combining transcripts ranges for each gene to get genes
   genes <- unlist(range(trx))
   genes <- alter.seqnames4granges(genes)
   genes <- add.seqlengths(genes, chromL)
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

   promoters <- promoters(my.annotation$genes, upstream=2000, downstream=500)
   return(list(exons=exons, genes=genes, promoters=promoters, genic=genic, intergenic=intergenic, intergenic.str=intergenic.str, introns=intronsByGene))
}
