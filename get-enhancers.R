#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
suppressMessages(library(argparse))
suppressMessages(library(tools))
source("~/source/Rscripts/annotation-functions.R")

parser <-  ArgumentParser(description="Define enhancer elements")
parser$add_argument('-c', '--config', metavar= "file", required='True', type= "character", help= "config file")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args <- parser$parse_args()
output_path <- file.path(args$out, 'enhancers')
plot_path <- file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

#}}}

# Load packages# {{{
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(biovizBase))
suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer))# }}}
source("~/local/functions.R")

config <- read.table(args$config, sep=",", header=TRUE, stringsAsFactors=TRUE)


#doesn't seperate based on condition
marks <- lapply(levels(config$Mark), function(i){
   expAU <- as.vector(config[config$Mark== i, 2])
   bed <- as.vector(config[config$Mark== i, 4])
   names(bed) <- expAU
   myGRanges <- lapply(bed, function(b){import.bed(con=b)})
})
names(marks) <- levels(config$Mark)

#For all marks take the union from different experiments
red_marks <- lapply(marks, function(i){Reduce(union, i)})

enhancer_marks <- c('p300', 'H3K27ac', 'H3K4me1')
promoter_marks <- 'H3K4me3'

red_marks[[paste(enhancer_marks, collapse='_')]] <- Reduce(intersect, red_marks[enhancer_marks])
#could potentially do a union of any promoter marks if more than one but no need now
red_marks[[paste(promoter_marks, collapse='u')]] <- Reduce(union, red_marks[promoter_marks])

red_marks[[paste(paste(enhancer_marks, collapse='_'), promoter_marks, sep="-")]] <- red_marks[[paste(enhancer_marks, collapse='_')]][! red_marks[[paste(enhancer_marks, collapse='_')]] %over% red_marks[[paste(promoter_marks, collapse='u')]]]

general.crap<-function(){
setwd("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers")

chromL<-read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10_noHead.genome", sep="\t", row.names=1)
mydirs <- list.dirs(path=getwd(), recursive=FALSE)
names(mydirs)<-gsub(".+/", '', mydirs)

#Removing the compendium files as I don't know where that bed comes from
#mydirs<-mydirs[grep("compendium", mydirs, invert=TRUE)]

ListFilesMarks <- sapply(mydirs, function(x){list.files(path=paste(x, "macs/mm10", sep="/"), full.names=TRUE)})

#mod <- lapply(ListFilesMarks, function(x){ list<-lapply(x[!grepl("unlifted", x)], function(i){ list<-read.table(i, colClasses=c(rep(NA, 3), "NULL", NA), col.names=c('chr', 'start', 'end', 'MACS_peak', 'score'))}); names(list)<- gsub("\\..+", "", gsub(".+/", "", x[!grepl("unlifted", x)])); return(list) })
#names(mod)<-names(ListFilesMarks)

mod <- lapply(ListFilesMarks, function(x){ list<-lapply(x[!grepl("unlifted", x)], function(i){ list<-read.table(i, sep="\t")}); names(list)<- gsub("\\..+", "", gsub(".+/", "", x[!grepl("unlifted", x)])); return(list) })
names(mod)<-names(ListFilesMarks)
mod<-lapply(mod, function(x){lapply(x, function(i){ n<-c('chr', 'start', 'end', 'MACS_peak', 'score'); if(ncol(i)==3) {n<-c('chr', 'start', 'end')}; names(i)<-n; if(ncol(i)==3){i$score<-rep(NA, nrow(i))}; return(i)})})


modGR <- lapply(mod, function(x){lapply(x, function(i){with(i, GRanges(chr, IRanges(start, end), rep("*", length(chr)), score=score))})})


#Med12<-modGR$compendium$Med12xYoung

#Creyghton DS 39900 peaks, xiao DS 11190 and intersection 15000
K27ac <- union(modGR$creyghton$H3K27ac, modGR$xiao$H3K27ac)

#Creyghton DS has 58200 peaks, Xiao DS has 1900 and intersection gives 2040...should I still join the two DS?
K4me1 <- union(modGR$creyghton$H3K4me1, modGR$xiao$H3K4me1)

#aEnh <- intersect(modGR$creyghton$p300, K27ac)
aEnh <- modGR$creyghton$p300[modGR$creyghton$p300 %over% K27ac]
aEnh<- aEnh[aEnh %over% Med12]

#Should I be using the DS in 2i or serum? Or separately? Use the union for the identification of enhancers to be as inclusive as possible--OK
mrK4me3 <- union(modGR$marks$H3K4me3c2i, modGR$marks$H3K4me3cserum)
K4me3 <- union(modGR$creyghton$H3K4me3, modGR$xiao$H3K4me3)
K4me3 <- union(K4me3, mrK4me3)

pEnh <- intersect(aEnh, K4me1)
seqlengths(pEnh)<-chromL[names(seqlengths(pEnh)),]
enh<-pEnh[!pEnh %over% K4me3]
#Naming the enhancers identified
names(enh)<-paste("E", 1:length(enh), sep="")

prLike<-pEnh[pEnh %over% K4me3]
names(prLike)<-paste("PrLE", 1:length(prLike), sep="")
###pEnh=enh + prLike
}

#Up and working fine! Could allow to pass the chrom.length as a dataframe too so I don't have to read it up every time. will have to check the class basically -> character it's a file otherwise it should be df
#chrom.length should be the file e.g. "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10_noHead.genome"
add.seqlengths <- function(grange, chrom.length){

   if( !file.exists(chrom.length) ) stop(chrom.length, " doesn't exist.")
   chrom.length <- read.table(chrom.length, sep="\t", row.names=1)
   seqlengths(grange) <- chrom.length[names(seqlengths(grange)),]
   return(grange)
}

alter.seqnames4granges<-function(grange){
   require('gtools')

   #the 2 functions could be combined in an if else statement to assignt he values to chr and new and then just do seqlevels...
   #ToDO: add a line for the chrM and MT conversion
   chr2num <- function(grange){
      chr <- mixedsort(seqlevels(grange))
      new <- mixedsort(gsub('chr', '',seqlevels(grange)))
      names(new) <- paste("chr", new, sep = "")

      seqlevels(grange) <- chr
      grange <- keepSeqlevels(grange, chr)
      grange <- renameSeqlevels(grange, new)
      return(grange)
   }

   num2chr <- function(grange){
      chr <- mixedsort(seqlevels(grange))
      new <- mixedsort(gsub("^", 'chr', seqlevels(grange)))
      names(new) <- mixedsort(seqlevels(grange))

      seqlevels(grange) <- chr
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


gcrap<-function(){
chr.sub <- paste("chr", 1:19, sep = "")
chr.sub<-c(chr.sub, "chrX", "chrY")
new.names<-c(as.character(1:19), 'X', 'Y')
names(new.names) <- paste("chr", new.names, sep = "")

#Could be replaced with a small function to make the code more readable
seqlevels(enh)<-chr.sub
enh<-keepSeqlevels(enh, chr.sub)
enh<-renameSeqlevels(enh, new.names)

seqlevels(pEnh)<-chr.sub
pEnh<-keepSeqlevels(pEnh, chr.sub)
pEnh<-renameSeqlevels(pEnh, new.names)

seqlevels(prLike)<-chr.sub
prLike<-keepSeqlevels(prLike, chr.sub)
prLike<-renameSeqlevels(prLike, new.names)


#Print enhancers granges as bed file so you can use bedtools intersect and work around the loosing the metadata bit
export.bed(enh, "enhancers.bed")
export.bed(enh, "PrLikeEnhancers.bed")



pdf("EnhWidth.pdf")
boxplot(w~s, data=data.frame(s=seqnames(sortSeqlevels(pEnh)), w=width(sortSeqlevels(pEnh))/1000), ylab="Width in Kb", main="pEnhancers")
boxplot(w~s, data=data.frame(s=seqnames(sortSeqlevels(enh)), w=width(sortSeqlevels(enh))/1000), ylab="Width in Kb", main="Enhancers")
boxplot(w~s, data=data.frame(s=seqnames(sortSeqlevels(prLike)), w=width(sortSeqlevels(prLike))/1000), ylab="Width in Kb", main="prLike Enhancers")
dev.off()

fwdFiles <- grep("_fwd", list.files("/nfs/research2/bertone/user/mxenoph/hendrich/bed_files", full.names=TRUE), value=TRUE)
revFiles <- grep("_rev", list.files("/nfs/research2/bertone/user/mxenoph/hendrich/bed_files", full.names=TRUE), value=TRUE)

fwd <- lapply(fwdFiles, function(x){read.table(x, col.names=c('chr', 'start', 'end', 'BEDscore'))})
names(fwd) <-  gsub("_Epi_fwd.+", "", gsub(".+/", "Epi_", fwdFiles))
fwdGR <- lapply(fwd, function(x){with(x, GRanges(chr, IRanges(start, end), rep("+", length(chr)), score=BEDscore))})
names(fwdGR) <-  gsub("_Epi_fwd.+", "", gsub(".+/", "Epi_", fwdFiles))

rev <- lapply(revFiles, function(x){read.table(x, col.names=c('chr', 'start', 'end', 'BEDscore'))})
names(rev) <-  gsub("_Epi_rev.+", "", gsub(".+/", "Epi_", revFiles))
revGR <- lapply(rev, function(x){with(x, GRanges(chr, IRanges(start, end), rep("-", length(chr)), score=BEDscore))})
names(revGR) <-  gsub("_Epi_rev.+", "", gsub(".+/", "Epi_", revFiles))

#both.strands<-lapply(names(rev), function(i){combined<-rbind(fwd$i, rev$i); return(combined))

ensembl70<-useMart(host='jan2013.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
df <- getBM(attributes=c("ensembl_gene_id", "chromosome_name", "strand", "exon_chrom_start", "exon_chrom_end", "ensembl_exon_id"),  mart=ensembl70)
df <- df[grep("^[\\dXYMT]{0,2}$", df$chromosome_name, perl=TRUE),]

df$chromosome_name <- paste("chr", df[,2], sep="")
df[df$chromosome_name=="chrMT", 2]<-rep("chrM", length(df[df$chromosome_name=="chrMT", 2]))
exonsDF<-df
exons <- with(df, GRanges(chromosome_name, IRanges(exon_chrom_start, exon_chrom_end), strand, GeneID=ensembl_gene_id, ExonID=ensembl_exon_id))
seqlevels(exons, force=TRUE)<-chr.sub
exons<-keepSeqlevels(exons, chr.sub)
exons<-renameSeqlevels(exons, new.names)



#Can load TSS info which contains gene_id, strand, gene start and end
#genes<-getBM(attributes=c("ensembl_gene_id", "chromosome_name", "strand", "start_poistion", "end_position"),  mart=ensembl70)


#intronsT<-read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/Mus_musculus.GRCm38.70.intron_conv.gtf", col.names=c("chr", "biotype", "feature", "start", "end", "score", "strand", "frame", "att"), sep="\t")

#This is very slow, so I need to comment it out after one run and just load the object.
#for(i in 1:nrow(intronsT)){l<-unlist(strsplit(as.character(intronsT[i,9]), ';')); intronsT[i,'IntronNo']<-gsub(".+\\s", "", l[2]); intronsT[i,'Gene']<-gsub(".+\\s", "", l[1])}
#introns<-with(intronsT, GRanges(chr, IRanges(start, end), strand,  intronNo=IntronNo, Gene=Gene))

load("/nfs/research2/bertone/user/mxenoph/hendrich/Strand_specific/GeneCoords_TSS_Promoter.Rdata")
tss<-df
tss[tss$strand == '1', 'tss']<- tss[tss$strand == '1','start_position']
tss[tss$strand == '-1', 'tss']<- tss[tss$strand == '-1','end_position']
tss[,'Gene']<-rownames(tss)
tssGR<-with(tss, GRanges(chromosome_name, IRanges(tss, width=1), strand, Gene=Gene))
seqlevels(tssGR, force=TRUE)<-chr.sub

#genesDF<-df[, 1:4]
#genes<-with(genesDF, GRanges(chromosome_name, IRanges(start_position, end_position), strand, GeneID=rownames(genesDF)))


#Need to email GenomicRages author cause promoters() doesn't take into account the strand of the transcript it always returns range (start(x)-width, start(x)+width+1)
txdb<-makeTranscriptDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filters="", id_prefix="ensembl_", host="jan2013.archive.ensembl.org", miRBaseBuild=NA)
#promoters(txdb)
allTx <- transcriptsBy(txdb, by='gene')

#exonsGF<-exonsBy(txdb, by=c("tx", "gene"), use.names=TRUE)
exonsGF<-exonsBy(txdb, by=c("gene"))
exonsGF<-unlist(exonsGF)

introns<-intronsByTranscript(txdb, use.names=TRUE)
introns<-unlist(introns)

genes <- unlist(range(allTx))

seqlevels(genes, force=TRUE)<-chr.sub
genes<-keepSeqlevels(genes, chr.sub)
genes<-renameSeqlevels(genes, new.names)


export.bed(genes, "genes.bed")

#Then run bedtools intersect to get the names of the genes the enhancers overlap with if they do but keeping all the enhancers even if they don't overlap

genes.noStrand<-genes
strand(genes.noStrand) <- '*'
intergenic<-gaps(genes.noStrand)

#Keeps end() for - strand and start() for +.OK!
tss<-resize(genes, 1)
#Might want to take a window instead of one base.


#Not gonna be useful for the TSS
#Futr<-fiveUTRsByTranscript(txdb, use.names=TRUE)

lox_eR<-subsetByOverlaps(enh, fwdGR$Epi_2lox)
seqlevels(fwdGR$Epi_2lox)<-chr.sub
fwdGR$Epi_2lox<-keepSeqlevels(fwdGR$Epi_2lox, chr.sub)
fwdGR$Epi_2lox<-renameSeqlevels(fwdGR$Epi_2lox, new.names)


fwdGR$Epi_2lox[fwdGR$Epi_2lox %over% tss]
check<-intersect(lox_eR, exons)


#need to endoapply these things but it's giving me errors..lots of them
get.IntronExonBound<-function(GRList){
  st<-start(GRList);
  st<-st[2:length(st)];
#  as.vector(seqnames(GRList));
#  strand(GRList);
  #stGR<-GRanges(as.vector(seqnames(i)), IRanges(start= st, width=1), strand=as.vector(strand(i)), end<-end(GRList));
  end<-end(GRList);
  end<-end[1:length(end)-1];
  both<-c(st, end);
  return(both);
}



#intersectScF<-list.files(path=getwd(), pattern=".*intsc.*", full.names=TRUE)
intersectScF<-list.files(path=getwd(), pattern=".*wao.*", full.names=TRUE)
#bedInt<-lapply(intersectScF, function(i){read.table(i, col.names=c('chr', 'start', 'end', 'score', 'strand'))})
bedInt<-lapply(intersectScF, function(i){read.table(i, col.names=c('chr', 'start', 'end', 'name', 'score', 'strand', 'match.chr', 'match.start', 'match.end', 'match.name', 'match.score', 'match.strand', 'length'))})
names(bedInt)<- gsub("\\..*", '', gsub(".*/", '', intersectScF))
bedGR<-lapply(bedInt, function(x){with(x, GRanges(chr, IRanges(start, end), strand=strand, score=score))})
}
