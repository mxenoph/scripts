#!/usr/bin/env Rscript

# List of packages to load# {{{
x <- c('GenomicRanges', 'rtracklayer')
lapply(x, suppressMessages(library), character.only=T)
source("~/local/rangeoverlapper.R")
source("~/source/Rscripts/plotting-functions.R")
# }}}

config <- read.table("peaks.conf.txt", sep=",", header=TRUE, stringsAsFactors=TRUE)# {{{
chrom.length<-"/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/mm9.chrom.sizes"
#add the add.seqlengths function to a general script for sourcing it quickly

peaks <- lapply(levels(config$Treatment), function(i){
                sampleID <- paste(as.vector(config[config$Treatment== i, 'Condition']), as.vector(config[config$Treatment== i, 'Factor']), as.vector(config[config$Treatment== i, 'Replicate']), sep='_')
                bed <- as.vector(config[config$Treatment== i, 5])
                names(bed) <- sampleID
                myGRanges <- lapply(bed, function(b){import.bed(con=b)})
})
names(peaks) <- paste('cond', levels(config$Treatment), sep='_')

pdf("peaks.length.pdf")
lapply(names(peaks), function(treat){
       lapply(names(peaks[[treat]]), function(i){
              name<-paste(paste(treat, i, sep='_'), 'pdf', sep='.')
              boxplot(w~s, data=data.frame(s=seqnames(sortSeqlevels(peaks[[treat]][[i]])), w=width(sortSeqlevels(peaks[[treat]][[i]]))/1000), ylab="Width in Kb", main=paste("Peaks in", paste(treat, i, sep='_'), sep=" "))
              abline(h=c(0.5, 1),  col=c("blue","red"), lwd=2, lty="solid")
              })
})
dev.off()# }}}

#pdf("peaksEpiVsESVenn.pdf")# {{{
#Compare the peaks found in Epi with those in ES
lapply(seq_along(peaks)[1:length(peaks)-1], function(treat){
       v <- seq(treat+1, length(peaks))
       lapply(names(peaks[[treat]]), function(cond){

              toMatch <- paste(unlist(strsplit(cond, '_'))[1:2], collapse='_')
                            lapply(v, function(i){
                                 matching <- grep(toMatch, names(peaks[[i]]))


                                 if(length(matching) > 1){
                                    lapply(matching, function(x){
                                           label1 <-  paste(names(peaks[treat]), names(peaks[[treat]][cond]), sep="-")
                                           label2 <- paste(names(peaks[i]), names(peaks[[i]][x]), sep="-")
                                           id <- paste(label1, label2, sep="-")
                                           overlaps <- findOverlaps(peaks[[treat]][[cond]], peaks[[i]][[x]], type="any", select="first", minoverlap=100L)

                                           #myol <- olRanges(peaks[[treat]][[cond]], peaks[[i]][[x]], output="gr")
                                           #myol[elementMetadata(myol)[, "OLpercQ"] > 50]
                                           #myol[elementMetadata(myol)[, "OLpercS"] > 50]

                                           pdf(paste(id, '.pdf', sep=''))
                                           myVenn(length(peaks[[treat]][[cond]]), length(peaks[[i]][[x]]), sum(!is.na(overlaps)), c(label1, label2))
                                           dev.off()

                                    })
                                 }
                                 else{
                                   overlaps <- findOverlaps(peaks[[treat]][[cond]], peaks[[i]][[matching]], type="any", select="first", minoverlap=100L)
                                   label1 <-  paste(names(peaks[treat]), names(peaks[[treat]][cond]), sep="-")
                                   label2 <- paste(names(peaks[i]), names(peaks[[i]][matching]), sep="-")
                                   id <- paste(label1, label2, sep="-")


q
                                   print(id)
                                   pdf(paste(id, '.pdf', sep=''))
                                   #myol <- olRanges(peaks[[treat]][[cond]], peaks[[i]][[x]], output="gr")
                                   #myol[elementMetadata(myol)[, "OLpercQ"] > 50]
                                   #myol[elementMetadata(myol)[, "OLpercS"] > 50]

                                   plot_venn(length(peaks[[treat]][[cond]]), length(peaks[[i]][[matching]]), sum(!is.na(overlaps)), c(label1, label2))
                                   dev.off()
                                 }
                           })
         })
})
#dev.off()
# }}}
