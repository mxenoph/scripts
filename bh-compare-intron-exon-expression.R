source("/homes/mxenoph/local/functions.R")
library(ggplot2)

#For matching ensembl ids to gene names later on
ID2genename<-"/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/ensembl_ids_2mgi_symbol"
gene_names<-read.delim(ID2genename, as.is=TRUE, quote="")

#For generating plots on subsets of genes based on their expression--downreg, upreg, de
exonDE<-read.delim("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_EpivsKO_Epi-de.txt", as.is=TRUE, quote="")
exonDown<-exonDE[which(exonDE$baseMeanA > exonDE$baseMeanB), ]
exonUp<-exonDE[which(exonDE$baseMeanA < exonDE$baseMeanB), ]

dir<-"/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts"
#dir<-"/nfs/nobackup2/research/bertone/mxenoph/bh"
exons<-list.files(path=dir, pattern= glob2rx(paste("*Epi*sr.exon",".counts", sep="")), full.names=TRUE)
introns<-list.files(path=dir, pattern= glob2rx(paste("*Epi*sr.intron",".counts", sep="")), full.names=TRUE)

#exons<-list.files(path=dir, pattern= glob2rx(paste("2*Epi*sr.exon",".counts", sep="")), full.names=TRUE)
#introns<-list.files(path=dir, pattern= glob2rx(paste("2*Epi*sr.intron",".counts", sep="")), full.names=TRUE)

nexons<-gsub("\\..*s(s|r)", "", exons)
nexons<-gsub(".*/", "", nexons)

nintrons<-gsub("\\..*s(s|r)", "", introns)
nintrons<-gsub(".*/", "", nintrons)

counts_exons<-lapply(exons, read.table, row.names=1)
matrix_exons<-NULL
for(i in 1:length(counts_exons)){
# stopifnot(all(rownames(counts_exons[[i]])== rownames(counts_exons[[1]])))
 matrix_exons<-cbind(matrix_exons, counts_exons[[i]]$V2)

}
rownames(matrix_exons)<-rownames(counts_exons[[1]])
colnames(matrix_exons)<-nexons


counts_introns<-lapply(introns, read.table, row.names=1)
matrix_introns<-NULL
for(i in 1:length(counts_introns)){
# stopifnot(all(rownames(counts_introns[[i]])== rownames(counts_introns[[1]])))
 matrix_introns<-cbind(matrix_introns, counts_introns[[i]]$V2)

}
rownames(matrix_introns)<-rownames(counts_introns[[1]])
colnames(matrix_introns)<-nintrons

exon.gtf<-read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/Mus_musculus.GRCm38.70.gtf", as.is=TRUE, sep="\t")
intron.gtf<-read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/Mus_musculus.GRCm38.70.intron.gtf", as.is=TRUE, sep="\t")

ilengths<-read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/mouse_ensembl_intron_genelengths.txt")
imatlen<-matrix(ilengths[,2])
rownames(imatlen)<-ilengths[,1]


#get.lengths(intron.gtf, "mouse_ensembl_intron", "intron")
myreturn<-compute.FPKMS(matrix_introns, imatlen)
intron_fpkm<-myreturn$matrix

elengths<-read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/mouse_ensembl_exon_genelengths.txt")
ematlen<-matrix(elengths[,2])
rownames(ematlen)<-elengths[,1]

#get.lengths(exon.gtf, "mouse_ensembl_exon", "exon")
myreturn<-compute.FPKMS(matrix_exons, ematlen)
exon_fpkm<-myreturn$matrix

m<-(ncol(intron_fpkm))-1
for(i in 1:m){
 etmp<-NULL
 etmp<-matrix(exon_fpkm[,i])
 rownames(etmp)<-rownames(exon_fpkm)
 colnames(etmp)<-colnames(exon_fpkm)[i]

 ntmp<-NULL
 ntmp<-matrix(intron_fpkm[,i])
 rownames(ntmp)<-rownames(intron_fpkm)
 colnames(ntmp)<-colnames(intron_fpkm)[i]
 
 tmp<-NULL
 tmp<-merge(etmp, ntmp, by="row.names")
 rownames(tmp)<-tmp[,1]
 tmp<-tmp[, !colnames(tmp) %in% "Row.names"]
 tmp<-log(tmp+1)

 name<-gsub(" ", "", paste(gsub("\\..*", "", colnames(tmp)[1]), ".IntronVsExon.pdf"))
 xlabel<-gsub("\\.counts", "", colnames(tmp)[1])
 ylabel<-gsub("\\.counts", "", colnames(tmp)[2])

 pdf(name)
  p<-qplot(tmp[,1], tmp[,2], type="point", xlab=xlabel, ylab=ylabel)
  p<-p + aes(xmax=round(max(tmp)), ymax=(round(max(tmp))))
  p<-p + geom_point(alpha=1/10)
  p<-p + geom_smooth(method= "loess", size=1.5)
  print (p)
 dev.off() 

 #Adding 3 column --Names-- with the gene name for that id
 tmp<-cbind(tmp, Names=gene_names[match(rownames(tmp), gene_names$ensembl_gene_id), "mgi_symbol"])
 #Writing the normalized FPKMs to a table
 write.table(tmp, file=gsub(".pdf", ".txt", name), sep="\t", quote=FALSE, row.names=TRUE)
 #Also write the features for which intron FPKM is higher than exon
 write.table(tmp[tmp[,2] > tmp[,1] & tmp[,1] >1 & tmp[,2] >3,], file=gsub("IntronVsExon.pdf", "hIntronVsExon.txt", name), sep="\t", quote=FALSE, row.names=TRUE)


 tmp0<-tmp[rownames(tmp) %in% exonDE$id, !colnames(tmp) %in% "Names"]
 pdf(gsub(" ", "", paste(gsub("\\..*", "", colnames(tmp0)[1]), ".IntronVsExonRedDE.pdf")))
  p<-qplot(tmp0[,1], tmp0[,2], type="point", xlab=xlabel, ylab=ylabel)
  p<-p + aes(xmax=round(max(tmp0)), ymax=(round(max(tmp0))))
  p<-p + geom_point(alpha=1/10)
  p<-p + geom_smooth(method= "loess", size=1.5)
  print (p)
 dev.off()

 #Adding 3 column --Names-- with the gene name for that id
 tmp0<-cbind(tmp0, Names=gene_names[match(rownames(tmp0), gene_names$ensembl_gene_id), "mgi_symbol"])
 #Writing the normalized FPKMs to a table
 write.table(tmp0, file=gsub("IntronVsExon.pdf", "IntronVsExonRedDE.txt", name), sep="\t", quote=FALSE, row.names=TRUE)
 #Also write the features for which intron FPKM is higher than exon
 write.table(tmp0[tmp0[,2] > tmp0[,1] & tmp0[,1] >1 & tmp0[,2] >3,], file=gsub("IntronVsExon.pdf", "hIntronVsExonRedDE.txt", name), sep="\t", quote=FALSE, row.names=TRUE)


 tmp1<-tmp[rownames(tmp) %in% exonUp$id, !colnames(tmp) %in% "Names"]
 pdf(gsub(" ", "", paste(gsub("\\..*", "", colnames(tmp1)[1]), ".IntronVsExonRedUp.pdf")))
  p<-qplot(tmp1[,1], tmp1[,2], type="point", xlab=xlabel, ylab=ylabel)
  p<-p + aes(xmax=round(max(tmp1)), ymax=(round(max(tmp1))))
  p<-p + geom_point(alpha=1/10)
  p<-p + geom_smooth(method= "loess", size=1.5)
  print (p)
 dev.off()


 tmp2<-tmp[rownames(tmp) %in% exonDown$id, !colnames(tmp) %in% "Names"]
 pdf(gsub(" ", "", paste(gsub("\\..*", "", colnames(tmp2)[1]), ".IntronVsExonRedDown.pdf")))
  p<-qplot(tmp2[,1], tmp2[,2], type="point", xlab=xlabel, ylab=ylabel)
  p<-p + aes(xmax=round(max(tmp2)), ymax=(round(max(tmp2))))
  p<-p + geom_point(alpha=1/10)
  p<-p + geom_smooth(method= "loess", size=1.5)
  print (p)
 dev.off()

}

