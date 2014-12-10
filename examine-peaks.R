require(biomaRt)
require(gridExtra)
require(GenomicRanges)
require(ggplot2)
require(plyr)
require(gridExtra)
require(gsubfn) #for selecting substrings in perl-like way
require(VennDiagram)
require(Vennerable)

source("~/local/functions.R")
load("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/macs/mm9.e67.ann.Rdata")
setwd("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/macs")

e67.mart <- useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
attr <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"), mart=e67.mart, uniqueRows=T)
attr[attr$mgi_symbol=='', 'mgi_symbol'] <- attr[attr$mgi_symbol=='', 'ensembl_gene_id']
attr[duplicated(attr$ensembl_gene_id), 'mgi_symbol'] <- attr[duplicated(attr$ensembl_gene_id),'ensembl_gene_id']

#Get data
db.files <- list.files(pattern= "*diffbind.csv")
names(db.files) <- unlist(lapply(db.files, function(x){ s <- strapplyc(x, ".*_(.*)([[:upper:]].*)ChIP.*", simplify=TRUE); paste(s[2], s[1], sep='_')}))

db.gr <- lapply(db.files, function(f){
       x <- read.table(f, header=TRUE, sep=',')
       x <- with(x, GRanges(Chr, IRanges(Start, End), strand='*', conc.wt= Conc_WT, fold=Fold, pval=p.value, fdr=FDR))
       names(x) <- paste('db', seq(1:length(x)), sep='')
       x
})

peak.files <- list.files(pattern= paste(".*", unique(strapplyc(names(db.gr), "(.*)_", simplify=TRUE)), ".*peaks.bed", sep=''))
names(peak.files) <- gsub('-', '_', gsub("^[0-9]", '', strapplyc(peak.files, "(.*)_.*peaks*")))
tmp <- list.files(pattern="*Epi.*peaks.bed")
names(tmp) <- gsub('-', '_', gsub("^[0-9]", '', strapplyc(tmp, "(.*)..*peaks*")))
peak.files <- c(peak.files, tmp)

mbd3.peak.files <- list.files(pattern="M2.*peaks.bed")
names(mbd3.peak.files) <- gsub('-', '_', gsub("^[0-9]", '', strapplyc(mbd3.peak.files, "(.*)..*peaks")))

mk.peaks.gr <- function(f){
   x <- read.table(f, col.names= c('chr', 'start', 'end', 'macs', 'score.10logp'))
   x <- with(x, GRanges(chr, IRanges(start, end), score=score.10logp, macs=macs))
   names(x) <- x$macs
   x
}

peak.gr <- lapply(peak.files, mk.peaks.gr)
mbd3.peak.gr <- lapply(mbd3.peak.files, mk.peaks.gr)

#mm10
#de.files <- list.files("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/", pattern="WT.*vsKO[^-]*.txt", full.names=TRUE)

de.files <- list.files("/nfs/research2/bertone/user/mxenoph/hendrich/deseq/", pattern="WT.*vsKO[^-]*.txt", full.names=TRUE)
names(de.files) <- gsub('.*WT', 'WTvsKO', gsub('vs.*', '', de.files))
de.df <- lapply(de.files, function(f){
                x <- read.table(f, header=TRUE, sep="\t")
                rownames(x) <- x$id
                x
})

#grS should be a gr for genes
annotate.nearBy <- function(grQ, grS, fun, de=NULL){
   if(fun == 'nearest'){
      d <- as.data.frame(distanceToNearest(grQ, grS))
   }
   else if(fun == 'precedes'){
      v <- precede(grQ, grS)
   }
   else if(fun== 'follows'){
      v <- follow(grQ, grS)
   }

   if(exists("d")) {
      elementMetadata(grQ)['distanceTonear'] <- d$distance
      v <- d$subjectHits
   }
   if(is.null(de)){
      elementMetadata(grQ)[fun] <- names(grS)[v];
      return(grQ);
   }
   elementMetadata(grQ)[fun] <- names(grS)[v];

   FC <- de[elementMetadata(grQ)[[fun]], 'log2FoldChange']
   FC[is.na(FC)] <- 0
   elementMetadata(grQ)[paste(fun, '.logFC', sep='')] <- FC

   padj <- de[elementMetadata(grQ)[[fun]], 'padj']
   elementMetadata(grQ)[paste(fun, '.padj', sep='')] <- padj

   group <- padj; group[group < 0.05] <- '<0.05'; group[group > 0.05] <- '>0.05'; group[is.na(group)] <- 'notExp'
   group <- factor(group, level= c('<0.05','>0.05','notExp'))
   elementMetadata(grQ)[paste(fun, '.fdr', sep='')] <- group

   names  <- de[elementMetadata(grQ)[[fun]], 'Name']
   elementMetadata(grQ)[paste(fun, '.name', sep='')] <- names

   return(grQ)
}

annotate.peak.near <- function(grQ, grS, cutoff){
   p <- as.data.frame(distanceToNearest(grQ, grS))
   elementMetadata(grQ)[p[p$distance <= cutoff, 'queryHits'], paste('near.peak.', cutoff, 'bp.score', sep='')] <- elementMetadata(grS)[p[p$distance <= cutoff, 'subjectHits'], 'score' ]
   elementMetadata(grQ)[p[p$distance > cutoff, 'queryHits'], paste('near.peak.', cutoff, 'bp.score', sep='')] <- 0
   return(grQ)
}

EpiSC <- annotate.nearBy(db.gr$Mi2b_Epi, mm9.e67$genes, 'nearest', de.df$WTvsKO_Epi)
EpiSC <- annotate.peak.near(EpiSC, mbd3.peak.gr$E12.M2.Epi, 200)

serumESC <- annotate.nearBy(db.gr$Mi2b_serum, mm9.e67$genes, 'nearest', de.df$WTvsKO_Lif)
serumESC <- annotate.peak.near(serumESC, mbd3.peak.gr$E12_M2, 200)

ESC2i <- annotate.nearBy(db.gr$Mi2b_2i, mm9.e67$genes, 'nearest', de.df$WTvsKO_2i)
ESC2i <- annotate.peak.near(ESC2i, mbd3.peak.gr$E12_2i_M2, 200)

peak.gr <- sapply(names(peak.gr), function(n){ peak.gr[[n]] <- annotate.nearBy(peak.gr[[n]], mm9.e67$genes, 'nearest')})
mbd3.peak.gr <- sapply(names(mbd3.peak.gr), function(n){ mbd3.peak.gr[[n]] <- annotate.nearBy(mbd3.peak.gr[[n]], mm9.e67$genes, 'nearest')})
nurd <- intersect()

## ggplot2 functions
gg_color_hue <- function(n) {
   hues = seq(15, 375, length=n+1)
   hcl(h=hues, l=65, c=100)[1:n]
}

#modify to count the args and call the right draw venn function
my.venn <- function(areas, intersections, cust.lab, title){
   require(VennDiagram)
   require(RColorBrewer)
   draw.triple.venn.custom <- dget("/homes/mxenoph/local/draw.triple.venn.custom.R")
   col<-brewer.pal(length(areas), 'Set1')

   venn <- draw.triple.venn.custom(
                            area1= areas[1],
                            area2= areas[2],
                            area3= areas[3],
                            n12= intersections[1],
                            n23= intersections[2],
                            n13= intersections[3],
                            n123= intersections[4],
                            fill= col,
                            alpha= 0.5,
                            category= names(areas),
                            euler.d = TRUE,
                            scaled= TRUE,
                            lty = "blank",
                            label.col =rep('grey', (length(areas)*2)+1),
                            cex = 1.5,
                            cat.cex = 1.5,
                            #fontfamily= 2,
                            #cat.fontfamily= rep('arial', length(areas)),
                            #cat.pos = c(),
                            cat.dist = 0.15,
                            #cat.just = list(c(-1, -1), c(1, 1)),
                            ext.pos = 30,
                            ext.dist = -0.05,
                            ext.length = 0.85,
                            ext.line.lwd = 2,
                            ext.line.lty = "dashed",
                            inter.lab = cust.lab,
                            main= title,
                            main.col= 'black'
                            );
   grid.draw(venn);

}

g_legend<-function(a.gplot){
   tmp <- ggplot_gtable(ggplot_build(a.gplot))
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[leg]]
   legend
}

#plot.data could be a grange or a data.frame
peak.multiplot <- function(plot.data){
   c <- gg_color_hue(1)
   plot.data <- as.data.frame(plot.data)
   cor.plot <- ggplot(plot.data, aes(x=fold, y=nearest.logFC))
   cor.plot <- cor.plot + geom_point(aes(colour=nearest.fdr, size=near.peak.200bp.score)) + scale_colour_discrete(drop = FALSE)
   print(cor.plot)
   grid.newpage()
   if (nrow(subset(plot.data, nearest.fdr == '<0.05')) > 0  ){
      cor.plot <- cor.plot + geom_text(data= subset(plot.data, nearest.fdr == '<0.05'), aes(label= nearest.name, position= "identity"), colour=c, size= 4, show_guide=F)
   }
   cor.legend <- g_legend(cor.plot)
   d.boxplot <- ggplot(plot.data, aes(x=nearest.fdr, y=distanceTonear / 1000))
   d.boxplot <- d.boxplot + geom_boxplot(aes(group=nearest.fdr, colour=nearest.fdr, fill=nearest.fdr), alpha=0.3) + theme(legend.position="none") #+coord_flip()

   blankPanel <- grid.rect(gp=gpar(col="white"))
   heights <- c(2, 5)/ sum(2,5)
   widths <- c(2.5/10, 7.5/10)
   p <- grid.arrange(blankPanel, d.boxplot, cor.legend, cor.plot + theme(legend.position="none"),  ncol=2, nrow=2, widths=widths, heights=heights, newpage=FALSE)
   print(p)
}

pdf("DBvsDEplusMbd3.pdf")
peak.multiplot(EpiSC)
grid.newpage()
peak.multiplot(serumESC)
grid.newpage()
peak.multiplot(ESC2i)
dev.off()

###################### draw venns

comp.inter <- function(grL, title){
   mi2b.wtko.o <- findOverlaps(grL[['E12.Chd4']], grL[['G9.Chd4']])
   mbd3.mi2bKO.o <- findOverlaps(grL[['E12.M2']], grL[['G9.Chd4']])
   mbd3.mi2bWT.o <- findOverlaps(grL[['E12.M2']], grL[['E12.Chd4']])

   tmp <- queryHits(mbd3.mi2bWT.o); tmp.1 <- grL[['E12.M2']]; tmp.2 <- queryHits(mbd3.mi2bWT.o);
   if(length(unique(queryHits(mbd3.mi2bWT.o))) > length(unique(subjectHits(mbd3.mi2bWT.o)))) {tmp <- subjectHits(mbd3.mi2bWT.o); tmp.1 <- grL[['E12.Chd4']]; tmp.2 <- subjectHits(mbd3.mi2bWT.o) }
   nurd <- tmp.1[1:length(tmp.1) %in% tmp.2]
   nurd.mi2bKO.o <- findOverlaps(nurd, grL[['G9.Chd4']])

#   nurd <- grL[['E12.M2']][1:length(grL[['E12.M2']]) %in% queryHits(mbd3.mi2bWT.o)]
#   nurd.mi2bKO.o <- findOverlaps(nurd, grL[['G9.Chd4']])
#
   areas <- c(length(grL[['E12.Chd4']]), length(grL[['G9.Chd4']]), length(grL[['E12.M2']]))
   names(areas) <- c('Mi2b.WT','Mi2b.KO','Mbd3')
   intersections <- c(min(length(unique(queryHits(mi2b.wtko.o))), length(unique(subjectHits(mi2b.wtko.o)))),
                      min(length(unique(queryHits(mbd3.mi2bKO.o))), length(unique(subjectHits(mbd3.mi2bKO.o)))),
                      min(length(unique(queryHits(mbd3.mi2bWT.o))), length(unique(subjectHits(mbd3.mi2bWT.o)))),
                      min(length(unique(queryHits(nurd.mi2bKO.o))), length(unique(subjectHits(nurd.mi2bKO.o))))
                      )
   non.uniq <- c(title,
                 paste0('\nmi2b.wt:', length(queryHits(mi2b.wtko.o)) - length(unique(queryHits(mi2b.wtko.o))),
                        '\nmi2b.ko:', length(subjectHits(mi2b.wtko.o)) - length(unique(subjectHits(mi2b.wtko.o)))),
                 '',
                 paste0('\nmbd3:', length(queryHits(mbd3.mi2bWT.o)) - length(unique(queryHits(mbd3.mi2bWT.o))),
                        '\nmi2b.wt:', length(subjectHits(mbd3.mi2bWT.o)) - length(unique(subjectHits(mbd3.mi2bWT.o)))),
                 '',
                 paste0('\nmbd3:', length(queryHits(mbd3.mi2bKO.o)) - length(unique(queryHits(mbd3.mi2bKO.o))),
                        '\nmi2b.ko:', length(subjectHits(mbd3.mi2bKO.o)) - length(unique(subjectHits(mbd3.mi2bKO.o)))),
                 ''
                 )

   my.venn(areas, intersections, non.uniq, title)
}

pdf('testme.pdf')
grL <- GRangesList('E12.Chd4'= peak.gr$E12.Chd4.Epi, 'G9.Chd4'=peak.gr$G9.Chd4.Epi, 'E12.M2'= mbd3.peak.gr$E12.M2.Epi)
comp.inter(grL, ' Epi')
grid.newpage()

grL <- GRangesList('E12.Chd4'= peak.gr$E12_2i_Mi2b, 'G9.Chd4'=peak.gr$G9_2i_Mi2b, 'E12.M2'= mbd3.peak.gr$E12_2i_M2)
comp.inter(grL, ' 2i')
grid.newpage()
grL <- GRangesList('E12.Chd4'= peak.gr$Mi2b_1, 'G9.Chd4'=peak.gr$G9_Mi2b, 'E12.M2'= mbd3.peak.gr$E12_M2)
comp.inter(grL, ' serum')
grid.newpage()

dev.off()
####################

plur.gene.name <- c('nanog', 'pou5f1', 'zfp42', 'fgf4', 'dppa4', 'dppa3', 'lin28a', 'nr0b1', 'esrrb', 'mycn', 'dppa2', 'sall4', 'prdm14', 'dazl', 'tcl1', 'klf2', 'tbx3', 'klf5', 'klf4', 'lin28b', 'crebbp', 'id3', 'utf1', 'id2', 'tert', 'myc')
plur.id <- sapply(plur.gene.name, function(p){as.character(de.df$WTvsKO_Epi$id[grep(paste("^", p, "$", sep=''), tolower(de.df$WTvsKO_Epi$Name))])})
plur.genes <- mm9.e67$genes[plur.id]
values(plur.genes)<- data.frame(gene.name=plur.gene.name)
#check if any of the peaks are near plur.associated genes--better to do the other way round cause that way you don't lose info because another gene is nearest to the peak

nurd.epi <- intersect(peak.gr$E12.M2.Epi, peak.gr$E12.Chd4.Epi)
plr.genes.df <- as.data.frame(plur.genes)
plr.genes.df[plr.genes.df$strand== '+', c('start')] <- plr.genes.df[plr.genes.df$strand== '+', c('start')] - 5000
plr.genes.df[plr.genes.df$strand== '+', c('end')] <- plr.genes.df[plr.genes.df$strand== '+', c('end')] + 500

plr.genes.df[plr.genes.df$strand== '-', c('start')] <- plr.genes.df[plr.genes.df$strand== '-', c('start')] - 500
plr.genes.df[plr.genes.df$strand== '-', c('end')] <- plr.genes.df[plr.genes.df$strand== '-', c('end')] + 5000

plur.genes.ext <- with(plr.genes.df, GRanges(seqnames, IRanges(start, end), strand=strand, gene.name=gene.name))
o <- findOverlaps(nurd.epi, plur.genes.ext)
d <- data.frame(nurd.peak=rep(0, length(plur.genes.ext)))
d[unique(subjectHits(o)), 'nurd.peak'] <- rep(1, length(unique(subjectHits(o))))
values(plur.genes.ext) <- cbind(values(plur.genes.ext), d)

o <- findOverlaps(peak.gr$G9.Chd4.Epi, plur.genes.ext)
d <- data.frame(mi2b.KO=rep(0, length(plur.genes.ext)))
d[unique(subjectHits(o)), 'mi2b.KO'] <- rep(1, length(unique(subjectHits(o))))
values(plur.genes.ext) <- cbind(values(plur.genes.ext), d)

serum.mi2b <- union(peak.gr$Mi2b_1, peak.gr$Mi2b_2)
nurd.serum <- intersect(mbd3.peak.gr$E12_M2, serum.mi2b)
o <- findOverlaps(nurd.serum, plur.genes.ext)
d <- data.frame(nurd.serum.peak=rep(0, length(plur.genes.ext)))
d[unique(subjectHits(o)), 'nurd.serum.peak'] <- rep(1, length(unique(subjectHits(o))))
values(plur.genes.ext) <- cbind(values(plur.genes.ext), d)

complexPeaks <- function(peaks.grL, ){
   require(GenomicRanges)

}

#length(nurd)
#length(peak.gr$E12.M2.Epi)


meta <-  data.frame(nrow=length(plur.genes))
for(x in names(peak.gr)){
   d <-as.data.frame(elementMetadata(annotate.nearBy(plur.genes, peak.gr[[x]], 'nearest')));
   colnames(d) <- paste(x, c('distanceTonear', 'nearest'), sep='.' );
   meta <- cbind(meta, d);
}
for(x in names(db.gr)){
   d <-as.data.frame(elementMetadata(annotate.nearBy(plur.genes, db.gr[[x]], 'nearest')));
   colnames(d) <- paste(x, c('distanceTonear', 'nearest'), sep='.' );
   meta <- cbind(meta, d);
}
values(plur.genes) <- cbind(values(plur.genes), meta, plur.gene.name)


#check for goterm enrichment in the nearest genes
library("org.Mm.eg.db")

s.fac.l <- lapply(list(EpiSC, serumESC), function(gr){
       s.fac <- factor(as.integer(as.data.frame(elementMetadata(gr))[['nearest.fdr']] == '<0.05' & as.data.frame(elementMetadata(gr))[['distanceTonear']] <= 500, level=c(0,1)));
       names(s.fac) <- as.data.frame(elementMetadata(gr))[['nearest']];
       s.fac
})
names(s.fac.l) <- c('EpiSC', 'serumESC')

GO.res.l <- lapply(seq(1:length(s.fac.l)), function(i){
                   info <- names(s.fac.l)[i]
                   GOtermEnr(s.fac.l[[i]], info, "BP", "org.Mm.eg.db")
})
names(GO.res.l) <- names(s.fac.l)


###### Mi2b only targets in ko
mi2b.only.KO.o <- findOverlaps(peak.gr$G9.Chd4.Epi, peak.gr$E12.M2.Epi)
mi2b.only.KO <- peak.gr$G9.Chd4.Epi[!(1:length(peak.gr$G9.Chd4.Epi) %in% queryHits(mi2b.only.KO.o))]

mi2b.only.KO.o <- findOverlaps(mi2b.only.KO, peak.gr$E12.Chd4.Epi)
mi2b.only.KO <- mi2b.only.KO[!(1:length(mi2b.only.KO) %in% queryHits(mi2b.only.KO.o))]

peak.gr$G9.Chd4.Epi <- annotate.nearBy(peak.gr$G9.Chd4.Epi, mm9.e67$genes, 'nearest')
mi2b.f <- factor(as.integer(names(peak.gr$G9.Chd4.Epi) %in% names(mi2b.only.KO)))
names(mi2b.f) <- peak.gr$G9.Chd4.Epi$nearest
mi2b.GO <- GOtermEnr(mi2b.f, 'mi2b.only', "BP", "org.Mm.eg.db")

peak.gr$g9.chd4.epi <- annotate.nearby(peak.gr$g9.chd4.epi, mm9.e67$genes, 'nearest')
mi2b.f <- factor(as.integer(names(peak.gr$g9.chd4.epi) %in% names(mi2b.only.ko)))
names(mi2b.f) <- peak.gr$g9.chd4.epi$nearest
mi2b.go <- gotermenr(mi2b.f, 'mi2b.only', "bp", "org.mm.eg.db")

peak.gr$E12.M2.Epi <- annotate.nearBy(peak.gr$E12.M2.Epi, mm9.e67$genes, 'nearest')
mbd3.f <- factor(as.integer(names(peak.gr$E12.M2.Epi) %in% names(mbd3.only)))
names(mbd3.f) <- peak.gr$E12.M2.Epi$nearest
mbd3.GO <- GOtermEnr(mbd3.f, 'mbd3.only', "BP", "org.Mm.eg.db")

exclusive <- list("mbd3.only.WT"=mbd3.GO, "mi2b.only.KO"= mi2b.GO)
pdf("DBvsDE.GO.pdf")
sapply(names(GO.res.l), function(i){
       t <- tableGrob(GO.res.l[[i]]$table, gpar.coretext =gpar(fontsize=8), gpar.coltext=gpar(fontsize=8), gpar.rowtext=gpar(fontsize=8))
       h <- grobHeight(t)
       w <- grobWidth(t)
       title <- textGrob(i, y=unit(0.5,"npc") + 0.55*h, vjust=0, gp=gpar(fontsize=10))
       gt <- gTree(children=gList(t, title))
       grid.draw(gt)
       grid.newpage()
})
sapply(names(exclusive), function(i){
       t <- tableGrob(exclusive[[i]]$table, gpar.coretext =gpar(fontsize=8), gpar.coltext=gpar(fontsize=8), gpar.rowtext=gpar(fontsize=8))
       h <- grobHeight(t)
       w <- grobWidth(t)
       title <- textGrob(i, y=unit(0.5,"npc") + 0.55*h, vjust=0, gp=gpar(fontsize=10))
       gt <- gTree(children=gList(t, title))
       grid.draw(gt)
       grid.newpage()
})
dev.off()

#selecting db loci based on length (1/4) of overlaps with mbd3 peaks
#d <- findOverlaps(epi.db.gr, mbd3.epi.WT.gr))
#wd <- width(pintersect(epi.db.gr[d$queryHits], mbd3.epi.WT.gr[d$subjectHits]))
#ok <- wd > width(epi.db.gr[queryHits(h)]) / 4
#d <- h[ok]



