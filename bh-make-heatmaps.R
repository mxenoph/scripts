source("~/source/Rscripts/functions.R")
source("~/source/Rscripts/deseq-functions.R")

setwd("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/TestMe/")

dir <- getwd() 

# Matrix with raw read counts for all samples
counts <- make_matrix(dir)

conds <- gsub(".*lox", "WT", colnames(counts))
conds <- gsub("(3KO)|(spl2)", "KO", conds)

cds <- newCountDataSet( counts, conds )
cds <- estimateSizeFactors( cds )
scaleFactors <- sizeFactors( cds )

length <- read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/mouse_ensembl_exon_genelengths.txt", sep="\t", row.names= 1)

# Compute FPKMs
fpkms <- compute.FPKMS.v2(counts, length, scaleFactors)

cols <- NULL
cols[grep("WT_2i", conds)] <- "#000080" 
cols[grep("WT_Epi", conds)] <- "#FF4500"
cols[grep("WT_Lif", conds)] <- "#228B22" 
cols[grep("WT_N2B27", conds)] <- "#4169E1" 
cols[grep("WT_NoLif", conds)] <- "#32CD32" 
cols[grep("KO_2i", conds)] <- "#483D8B" 
cols[grep("KO_Epi", conds)] <- "#FF8C00" 
cols[grep("KO_Lif", conds)] <- "#556B2F" 
cols[grep("KO_N2B27", conds)] <- "#6A5ACD" 
cols[grep("KO_NoLif", conds)] <- "#DAA520" 

leg <- c("#000080", "#4169E1", "#483D8B", "#6A5ACD", "#228B22", "#32CD32", "#556B2F", "#DAA520", "#FF4500", "#FF8C00")
names(leg) <- c("WT_2i", "WT_N2B27", "KO_2i", "KO_N2B27", "WT_Lif", "WT_NoLif", "KO_Lif", "KO_NoLif", "WT_Epi", "KO_Epi")

# Plot a pca
plot.pca(fpkms, cols, 'NuRD_WT&KO', leg)

deFiles <- list.files("/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/", pattern="*-de.txt", full.name=TRUE) 
mynames <- gsub("-de.*", "", deFiles)
mynames <- gsub(".*/", "", mynames) 
deGenes <- lapply(deFiles, read.table, sep="\t", header=T, row.names=1) 
names(deGenes) <- mynames 


for(i in 1:length(deGenes)){
 cond1 <- gsub("vs.*", "", names(deGenes)[i])
 cond2 <- gsub(".*vs", "", names(deGenes)[i])

 if (gsub("_.*", "", cond1) != "KO") { 
  col <- c(paste("2lox", gsub(".*_", "_", cond1), sep=''), paste("3Flox", gsub(".*_", "_", cond1), sep= '')) 
 } else {
  col <- c(paste("spl2", gsub(".*_", "_", cond1), sep=''), paste("3KO", gsub(".*_", "_", cond1), sep= ''))
 }
 if (gsub("_.*", "", cond2) != "KO") { 
  col2 <- c(paste("2lox", gsub(".*_", "_", cond2), sep=''), paste("3Flox", gsub(".*_", "_", cond2), sep= '')) 
 } else{
  col2<-c(paste("spl2", gsub(".*_", "_", cond2), sep=''), paste("3KO", gsub(".*_", "_", cond2), sep= ''))
 }

 matrix<-fpkms[,c(col,col2)]

 deGenes[[i]] <- deGenes[[i]][order(deGenes[[i]]$pval),]
 # sets NA for missing MGIs
 deGenes[[i]][deGenes[[i]]==""] <- NA
 deGenes[[i]] <- subset(deGenes[[i]], select= 8)   
 top50 <- rownames(deGenes[[i]])[1:50]

 draw.heatmap(matrix,top50, names(deGenes)[i], deGenes[[i]])
}

save(fpkms, file= "FPKMs.Rdata")

comparisons <- grep("WT.*WT.*", names(deGenes), value=T)
#Venn diagrams
library("VennDiagram")

for(i in comparisons){
 iKO <- gsub("WT", "KO", i)
 
 set1 <- deGenes[[i]]
 set2 <- deGenes[[iKO]]
 s1<-rownames(set1)
 s2<-rownames(set2)
 #alternative cross<-intersect(s1,s2)
 cross<-s2[s2 %in% s1]

 cond <- gsub("WT_", "", i)
 pdf(paste("VennWT-KO", cond, ".pdf", sep= ''), width= 12, height= 8)
 par(mar=c(2.1, 8.1, 2.1, 8.1))
 par(oma=c(1.5,1,1.5,1.5))

 venn.plot <- draw.pairwise.venn(
 area1 = length(s1),
 area2 = length(s2),
 cross.area = length(cross),
 scaled = TRUE,
 category = c(i, iKO),
 fill = c("blue", "red"),
 aplha = 0.5,
 lty = "blank",
 cex = 1.5,
 cat.cex = 1.5,
 cat.pos = c(285, 105),
 cat.dist = 0.09,
 cat.just = list(c(-1, -1), c(1, 1)),
 ext.pos = 30,
 ext.dist = -0.05,
 ext.length = 0.85,
 ext.line.lwd = 2,
 ext.line.lty = "dashed"
 );
 grid.draw(venn.plot);

 dev.off()
}

save(deGenes, file= "deGenes.Rdata")
