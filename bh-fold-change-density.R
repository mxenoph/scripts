
all.data <- read.table("/nfs/research2/bertone/user/remco/hendrich/chipseq/peaks/more/all.data.txt",header=T,as.is=TRUE,sep="\t")
twoi.data <- read.table("/nfs/research2/bertone/user/remco/hendrich/chipseq/peaks/2i/all.data_2i.txt",header=T,as.is=TRUE,sep="\t")

pdf("FC_upon_KO-Serum.pdf")
plot(density(all.data$LogFC),col="blue",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="All genes")
lines(density(all.data$LogFC[!is.na(all.data$pvalWTKO) & all.data$pvalWTKO < 0.05]),col="red",lwd=2)
legend("topright",legend=c("All","pval<0.05"),fill=c("blue","red"))
dev.off()

pdf("FC_upon_KO-Serum-NuRDbound.pdf")
plot(density(all.data$LogFC[all.data$nurd=="bound"]),col="blue",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="NuRD bound genes")
lines(density(all.data$LogFC[all.data$nurd=="bound" & !is.na(all.data$pvalWTKO) & all.data$pvalWTKO < 0.05]),col="red",lwd=2)
legend("topright",legend=c("All","pval<0.05"),fill=c("blue","red"))
dev.off()

pdf("FC_upon_KO-2i.pdf")
plot(density(twoi.data$LogFC),col="blue",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="All genes")
lines(density(twoi.data$LogFC[!is.na(twoi.data$pvalWTKO) & twoi.data$pvalWTKO < 0.05]),col="red",lwd=2)
legend("topright",legend=c("All","pval<0.05"),fill=c("blue","red"))
dev.off()

pdf("FC_upon_KO-2i-NuRDbound.pdf")
plot(density(twoi.data$LogFC[twoi.data$nurd=="bound"]),col="blue",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="NuRD bound genes")
lines(density(twoi.data$LogFC[twoi.data$nurd=="bound" & !is.na(twoi.data$pvalWTKO) & twoi.data$pvalWTKO < 0.05]),col="red",lwd=2)
legend("topright",legend=c("All","pval<0.05"),fill=c("blue","red"))
dev.off()

### add 2i mbd3/nurd binding

all.data$nurd.2i <- twoi.data$nurd[match(all.data$Gene,twoi.data$Gene)]
all.data$all.mbd3.2i <- twoi.data$all.mbd3[match(all.data$Gene,twoi.data$Gene)]

pdf("FC_upon_KO-Serum-NuRDbound_in_2i.pdf")
plot(density(all.data$LogFC[all.data$nurd.2i=="bound"]),col="blue",lwd=2,xlim=c(-4.3,4.3),xlab="Log fold change",main="NuRD bound genes in 2i")
lines(density(all.data$LogFC[all.data$nurd.2i=="bound" & !is.na(all.data$pvalWTKO) & all.data$pvalWTKO < 0.05]),col="red",lwd=2)
legend("topright",legend=c("All","pval<0.05"),fill=c("blue","red"))
dev.off()

