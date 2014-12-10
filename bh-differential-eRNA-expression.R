source("~/local/DEAFunctions.R")
source("~/local/enhancer.functions.R")
source("~/local/granges.functions.R")
require('gtools')
require(GenomicFeatures)
#### IMPORTANT: To get the enh object
load("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/enhancers.mm9.p300etal.Rdata")
#For getting the gene coordinates
load("mm10txdb.Rdata")


eMatx<-mkMatx4enh(getwd())
#Get the row names in order
require('gtools')
eMatx<-eMatx[mixedsort(rownames(eMatx)), ]
colnames(eMatx) <- gsub('_Epi', '',  gsub('.counts', '', colnames(eMatx)))
load("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/enhancers.p300etal.Rdata")
w <- width(enh)
names(w) <- names(enh)
# ordering the width object based on the order of the matrix
w <- w[rownames(eMatx)]

#this is a list..you should be using the $matrix for the heatmap
e.FPKMs <- compute.FPKMS(eMatx, w)


##### DE ###
plus.mat<-eMatx[, c('2lox.plus', '3Flox.plus', '3KO.plus', 'spl2.plus')]
conds <- c(rep("WT_Epi", 2), rep("KO_Epi", 2))

pcds<-newCountDataSet( plus.mat, conds )
pcds <- estimateSizeFactors( pcds )
pcds <- estimateDispersions(pcds)
pres<-nbinomTest(pcds, "WT_Epi", "KO_Epi")
presSig <- pres[!is.na(pres$padj) & pres$padj < .05, ]
presSig<-presSig[ order(presSig$pval), ]

min.mat<-eMatx[, c('2lox.minus', '3Flox.minus', '3KO.minus', 'spl2.minus')]
conds <- c(rep("WT_Epi", 2), rep("KO_Epi", 2))

mcds<-newCountDataSet( min.mat, conds )
mcds <- estimateSizeFactors( mcds )
mcds <- estimateDispersions(mcds)
mres<-nbinomTest(mcds, "WT_Epi", "KO_Epi")
mresSig <- mres[!is.na(mres$padj) & mres$padj < .05, ]
mresSig<-mresSig[ order(mresSig$pval), ]

#ionly used if you get call.de to work in enhancers.functions.R
des.min <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/enhancers_design4DE.minus.txt", header=TRUE, sep='\t', as.is=TRUE, quote="")
des.plus <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/enhancers_design4DE.plus.txt", header=TRUE, sep='\t', as.is=TRUE, quote="")
#####

minus.list <- mresSig$id
plus.list <- presSig$id

#Get the de genes so you can compare the enhancers spanning exons
de.genes <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi-de.txt", quote='', sep="\t", header=TRUE, row.names=1)
de.genes <- de.genes[order(-de.genes$log2FoldChange, -de.genes$baseMean),]

map.up.down <- function(vector){
   unlist( sapply(seq_along(vector), function(i) {
       #in case the enhancer is spannig multiple of the same feature type
       if (grepl(';', vector[i])) {
          ids <- unlist(strsplit(vector[i], ';'))
          ids <- match(ids, rownames(de.genes))
          values <- sapply(seq_along(ids), function(id){
                 if(is.na(ids[id])){v <-'notDE'}
                 else if(de.genes[id, 'log2FoldChange'] > 0) { v <- 'up'}
                 else { v <- 'down'}
          })
         #CC00CC", "#E0E0E0",  "#FF66B2", "#3333FF"
          vector[i] <- paste(values, collapse=';')
          if(length(unique(values)) > 1){
             names(vector[i]) <- "#CC00CC"
          }
          else{
             if(unique(values)== "up") {names(vector[i]) <- "#FF66B2"}
             else if(unique(values)== "down") {names(vector[i]) <- "#3333FF"}
          }
       }
       else if(!vector[i] == "NA"){
          if(!(vector[i] %in% rownames(de.genes))){
            vector[i] <- 'notDE'
            names(vector[i]) <- "#009900"
          }
          else if(de.genes[vector[i], 'log2FoldChange'] >0) {
             vector[i]<-'up'
             names(vector[i]) <- "#FF66B2"
          }
          else{
             vector[i]<-'down'
             names(vector[i]) <- "#3333FF"
          }
       }
source("~/local/granges.functions.R")
require(gtools)
require(GenomicFeatures)
#### IMPORTANT: To get the enh object
#load("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/mm10/enhancers.p300etal.Rdata")
#For getting the gene coordinates
#load("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10.e70.ann.Rdata")

#this doesn't work if you source the script- whic is unusual
load("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/mm9.e67.ann.Rdata")
#load("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/enhancers.mm9.p300etal.Rdata")

eMatx<-mkMatx4enh(paste0(getwd(), '/counts'))
#Get the row names in order
eMatx<-eMatx[mixedsort(rownames(eMatx)), ]
colnames(eMatx) <- gsub('_Epi', '',  gsub('.counts', '', colnames(eMatx)))
w <- width(enh)
names(w) <- names(enh)
# ordering the width object based on the order of the matrix
w <- w[rownames(eMatx)]

#this is a list..you should be using the $matrix for the heatmap
e.FPKMs <- compute.FPKMS(eMatx, w)


##### DE ###
ource("~/local/DEAFunctions.R")
source("~/local/enhancer.functions.R")
source("~/local/granges.functions.R")lus.mat<-eMatx[, c('2lox.plus', '3Flox.plus', '3KO.plus', 'spl2.plus')]
conds <- c(rep("WT_Epi", 2), rep("KO_Epi", 2))

pcds<-newCountDataSet( plus.mat, conds )
pcds <- estimateSizeFactors( pcds )
pcds <- estimateDispersions(pcds)
pres<-nbinomTest(pcds, "WT_Epi", "KO_Epi")
presSig <- pres[!is.na(pres$padj) & pres$padj < .05, ]
presSig<-presSig[ order(presSig$pval), ]

min.mat<-eMatx[, c('2lox.minus', '3Flox.minus', '3KO.minus', 'spl2.minus')]
conds <- c(rep("WT_Epi", 2), rep("KO_Epi", 2))

mcds<-newCountDataSet( min.mat, conds )
mcds <- estimateSizeFactors( mcds )
mcds <- estimateDispersions(mcds)
mres<-nbinomTest(mcds, "WT_Epi", "KO_Epi")
mresSig <- mres[!is.na(mres$padj) & mres$padj < .05, ]
mresSig<-mresSig[ order(mresSig$pval), ]

#ionly used if you get call.de to work in enhancers.functions.R
#des.min <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/enhancers_design4DE.minus.txt", header=TRUE, sep='\t', as.is=TRUE, quote="")
#des.plus <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/enhancers_design4DE.plus.txt", header=TRUE, sep='\t', as.is=TRUE, quote="")
#####

minus.list <- mresSig$id
plus.list <- presSig$id

#Get the de genes so you can compare the enhancers spanning exons
#de.genes <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi-de.txt", quote='', sep="\t", header=TRUE, row.names=1)
#de.genes <- de.genes[order(-de.genes$log2FoldChange, -de.genes$baseMean),]
rsc.mat.minus <- data.frame(# "intergenic" = mk.rsc.matx(map.col.strand, as.vector(as.data.frame(elementMetadata(enh[minus.list])['spanning.intergenic'])[[1]])),
#                      "exons" = mk.rsc.matx(map.exon.intron, as.vector(as.data.frame(elementMetadata(enh[minus.list])['spanning.exons'])[[1]])),
#                      "introns" = mk.rsc.matx(map.exon.intron, as.vector(as.data.frame(elementMetadata(enh[minus.list])['spanning.introns'])[[1]])),
#                      "DE.exons" = de.vectors4enh$spanning.exons$minus,
#                      "DE.introns" = de.vectors4enh$spanning.introns$minus,
                      "genes" = mk.rsc.matx(map.exon.intron, as.vector(as.data.frame(elementMetadata(enh[minus.list])['spanning.genes'])[[1]])),
                      "DE.genes" = de.vectors4enh$spanning.genes$minus,
                      "mi2b.peaks.WT" = mk.rsc.matx(map.bound, as.vector(as.data.frame(elementMetadata(enh[minus.list])['spanning.mi2b.WT.1'])[[1]])),
                      "mi2b.peaks.KO" = mk.rsc.matx(map.bound, as.vector(as.data.frame(elementMetadata(enh[minus.list])['spanning.mi2b.KO.1'])[[1]])),
                      "mbd3.peaks.WT" = mk.rsc.matx(map.bound, as.vector(as.data.frame(elementMetadata(enh[minus.list])['spanning.mbd3.WT.1'])[[1]]))
)
rownames(rsc.mat.minus) <- names(enh[minus.list])

#This works but ain't easy to read with all the rowside colours
draw.enhanc.heatmap(e.FPKMs$matrix[, grep('plus', colnames(e.FPKMs$matrix))], plus.list, 'enhancer.de.Epi.plus', as.matrix(rsc.mat.plus))
draw.enhanc.heatmap(e.FPKMs$matrix[, grep('minus', colnames(e.FPKMs$matrix))], minus.list, 'enhancer.de.Epi.minus', as.matrix(rsc.mat.minus))

rsc.mat.plus.list <- split(rsc.mat.plus, colnames(rsc.mat.plus))
rsc.mat.plus.minus <- split(rsc.mat.minus, colnames(rsc.mat.minus))


#txdb not working in the Rdata saved object but genes and stuff objects are still there
de.genes.gr <- my.annotation$genes[rownames(de.genes)]
values(de.genes.gr) <- data.frame('log2FC'= de.genes$log2FoldChange, 'padj' = de.genes$padj, 'Name' = de.genes$Name)

de.genes.gr <- match.annot(de.genes.gr, mi2b.epi.WT.gr, 'mi2b.WT')
de.genes.gr <- match.annot(de.genes.gr, mi2b.epi.KO.gr, 'mi2b.KO')
de.genes.gr <- match.annot(de.genes.gr, mbd3.epi.WT.gr, 'mbd3.WT')
