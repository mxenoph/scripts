#!/usr/bin/env Rscript

source("~/source/Rscripts/functions.R")# {{{
source("~/source/Rscripts/granges-functions.R")
x <- c('Repitools', 'GenomicRanges', 'rtracklayer', 'plyr', 'ggplot2', 'gridExtra', 'reshape')
lapply(x, suppressMessages(library), character.only=T)# }}}

load("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/mm9.e67.ann.Rdata")
gene_assoc_IDs <- read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/MM9.maps/mm9.e67.associated.geneIDs.txt", header=TRUE, as.is=TRUE, sep="\t")

args <- commandArgs(trailingOnly=TRUE)
peaks_1 <- read.table(args[1], header= TRUE, sep="\t")
peaks_2 <- read.table(args[2], header=TRUE, sep="\t")

common_targets <- as.vector(unlist(read.table(args[3])))
de <- args[4]
# Names for args[1] and args[2] separated by ;. e.g. 'SUZ12;MBD3'
pnames <- unlist(strsplit(args[5], ';'))
descr <- args[6]
out_dir <- args[7]

# Check if 'plots' dir exists in the wdr and create if not
plot_dir <- file.path(out_dir, 'plots', '')
dir.create(plot_dir)

peaks_1 <- "/nfs/research2/bertone/user/mxenoph/hendrich/chip/marks_2012/macs/2i-Suz12_peaks.xls"
peaks_2 <- "/nfs/research2/bertone/user/mxenoph/hendrich/chip/bh/chipseq/macs/7E12-2i-M2_peaks.xls"
de <- "/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal/WT_2ivsKO_2i.txt"
common_targets <- as.vector(unlist(read.table("/nfs/research2/bertone/user/mxenoph/hendrich/ChEA/M2-Chd4.2i.SUZ12_CommonTargets.txt")))
pnames <- unlist(strsplit('SUZ12;MBD3', ';'))
##

# TODO: Check if files are xls before using macs2GRanges# {{{
peaks <- lapply(c(peaks_1, peaks_2), function(x){
#peaks <- lapply(args[1:2], function(x){
                # For macs xls files
                if( grepl("\\.xls$", x) ){
                    macs2GRanges(x)
                }else if( grepl("\\.bed$", x) ){
                    import(x, format='bed')
                }
})
names(peaks) <- pnames
genes <- my.annotation$genes
extended_genes <- resize(promoters(genes, upstream=2000), width=width(genes)+2500)# }}}

# select for the common genes
keep <- gene_assoc_IDs[toupper(gene_assoc_IDs[["Associated.Gene.Name"]]) %in% common_targets, 'Ensembl.Gene.ID']
extended_genes <- extended_genes[keep,]

peaks_map_genes <- lapply(seq_along(peaks), function(ind){
                          x <- peaks[[ind]]
                          over <- findOverlaps(x, extended_genes)
                          x <- x[queryHits(over)]
                          df <- data.frame('Name'= names(extended_genes[subjectHits(over)]))
                          values(x) <- cbind(values(x), df)
                          x <- split(x, x$Name)

                          if(length(x) < length(common_targets)){
                              warning(paste0("Not all input targets have a peak in ", names(peaks)[ind]))
                              x
                          } #shouldn't it be greater than?
                          else if ( length(x) < length(extended_genes)) {
                              # Because one associated gene name can be mapped to more than
                              # one unique ensembl id
                              x
                          }
})
names(peaks_map_genes) <- pnames
# Returned genes are more than the common targets used as input. WHY ???
# Answer not all genes reported by ChEA to be targets of both MBD3 and 
# SUZ12 have a SUZ12 peak in the marks dataset
# 

# It's a bit weird that by intersect() you don't get the min number of 
# targets of the 2 ChIP. Look into it
dist <- lapply(intersect(names(peaks_map_genes[[1]]), names(peaks_map_genes[[2]])), function(name){
       summit_1 <<- sort(unlist(peaks_map_genes[[1]][[name]]$summit))
       summit_2 <<- sort(unlist(peaks_map_genes[[2]][[name]]$summit))

       min_dist <- NULL
       max_dist <- NULL
       if(length(summit_1) < length(summit_2)){
          small <- summit_1
          big <- summit_2
       } else {
           small <- summit_2
           big <- summit_1
       }

       for(i in 1:length(small)){
           tmp <- min(abs(big - small[i]))
           if (i ==1 || tmp < min_dist ) min_dist <- tmp

           tmp <- max(abs(big - small[i]))
           if (i ==1 || tmp > max_dist ) max_dist <- tmp
       }
       return(data.frame('Shortest.distance'= min_dist, 'Longest.distance' = max_dist, row.names=name))
})
tmp <- rbind.fill(dist)
rownames(tmp) <- unlist(lapply(dist, row.names))
dist <- tmp

# Subsets based on the 75th percentile
perc75 <- dist[dist$Shortest.distance < quantile(dist$Shortest.distance)[4], 'Shortest.distance']
de <- deseq2vect(de)
detar <- dist[rownames(dist) %in% rownames(de$up) | rownames(dist) %in% rownames(de$down), 'Shortest.distance', drop=F]
detar[rownames(detar) %in% rownames(de$up), 'DE'] <- 'up'
detar[rownames(detar) %in% rownames(de$down), 'DE'] <- 'down'

#detar75 <- detar[detar$Shortest.distance < quantile(detar$Shortest.distance)[4],]
# Subset based on the 75 percentile for the whole set not just the DE
detar75 <- detar[detar$Shortest.distance < quantile(dist$Shortest.distance)[4],]

pdf(paste0(plot_dir, paste(pnames, collapse='-'), ".",  descr, '.pdf'), paper='a4')
#p <- ggplot(data= melt(dist), aes(x=factor(variable), y=value/1000)) 
#p <- p + geom_violin(scale="count", colour="black", fill="black", adjust=0.5)

p <- ggplot(data=melt(dist[, 'Shortest.distance']), aes(x=value/1000))
p <- p + geom_histogram()
p <- p + theme_bw()
p <- p + theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background= element_blank(),
          panel.border= element_blank(),
          axis.line= element_line())


x <- ggplot(data=melt(detar75), aes(x=value/1000, fill=DE))
x <- x + geom_density(alpha=.3)
#x <- x + geom_histogram()
x <- x + theme_bw()
x <- x + theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background= element_blank(),
          panel.border= element_blank(),
          axis.line= element_line(),
          legend.position= "top")

grid.arrange(ggplotGrob(p), ggplotGrob(p %+% melt(perc75)), ggplotGrob(x),
            nrow=1, ncol=3,
            main= paste0("Distance (min) between ", paste(pnames, collapse="-"), " in ", descr))

grid.newpage()

de$df[is.na(de$df$padj), 'padj'] <- 1
de$df[, 'de'] <- de$df$padj < 0.05
de$df[, 'common'] <- rownames(de$df) %in% gene_assoc_IDs[toupper(gene_assoc_IDs[['Associated.Gene.Name']]) %in% common_targets, 'Ensembl.Gene.ID']

dev.off()

pp <- list('a'=p, 'b'=p %+% melt(perc75), 'c'=x)
save(pp, file=paste0("PlotObj-", descr, ".Rda"))

# do.call("grid.arrange", c(pp, ncol=2))

#selecting loci based on length (1/4) of overlaps
#wd <- width(pintersect(peaks[[1]][over$queryHits], genes[o$subjectHits]))
#ok <- wd > width(peaks[[1]][queryHits(h)]) / 4
#d <- over[ok]


