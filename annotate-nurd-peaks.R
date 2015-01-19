#!/usr/bin/env Rscript

# List of packages to load# {{{
x <- c('Repitools', 'ChIPpeakAnno', 'biomaRt', 'rtracklayer', 'org.Mm.eg.db')
lapply(x, suppressMessages(library), character.only=T)
source("~/source/Rscripts/granges-functions.R")# }}}

# Read in arguments -- replace with argsparse# {{{
if(TRUE){
args <- commandArgs(trailingOnly=TRUE)
#This has to be the prefixes before _peaks.xls separated by ; e.g.: 'E12.M2.Epi;E12.Chd4.Epi' always give mbd3 first (bacause of merging after)
conf <- read.table(args[1], header=T, as.is=T, quote="", colClasses=rep('character',4))
assembly <- args[2]
#This has to be the path, can't be $PWD when you run it with bsub
out <- args[3]
setwd(out)
}else{
   setwd("/nfs/research2/bertone/user/mxenoph/hendrich/chip/")
   conf <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/chip/NuRD.Epi.conf", header=T, as.is=T, quote="", colClasses=rep('character',4))
   assembly <- 'mm9'
}

if(grepl('^mm', assembly, ignore.case=TRUE)){
   ds <- 'mmusculus_gene_ensembl'

   if(grepl('mm10', assembly, ignore.case=TRUE)){
      load("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10.e70.ann.Rdata")
      chr.sizes <- "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10.chrom.sizes"
      host <- 'jan2013.archive.ensembl.org'
   } else if(grepl('mm9', assembly, ignore.case=TRUE)){
      load("/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/mm9.e67.ann.Rdata")
      chr.sizes <- "/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/mm9.chrom.sizes"
      host <- 'may2012.archive.ensembl.org'
   }
} else{
   stop("This script currently runs just for mouse data\n")
}# }}}

# Read in peaks from xls --needs to be updated# {{{
peaks <- lapply(conf$file, function(f){
                if(grepl('xls', f)){
                    p <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
                    p <- macs_to_granges(p)
                    p <- p[order(-values(p)$score)]
                    #p <- add.seqlengths(p, chr.sizes)
                } else{
                    p <- read.table(f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
                    p <- with(p, GRanges(V1, IRanges(start=V2, end=V3)))
                }
                return(p)
})
names(peaks) <- conf$factor# }}}

complex <- findOverlappingPeaks(as(peaks[[2]], "RangedData"), as(peaks[[1]], "RangedData"),# {{{
                                maxgap=200, select='first',
                                NameOfPeaks1=names(peaks)[1], NameOfPeaks2=names(peaks)[2],
                                annotate=1)

#cause I don't know what it's doing
#complex.anno <- annotatePeakInBatch(complex$MergedPeaks, AnnotationData=as(my.annotation$genes, "RangedData"))
export.bed(complex$MergedPeaks, con=paste0(paste(paste(names(peaks), collapse='-'), unique(conf$condition), sep="."), "_cmpl_peaks.bed"))
stop('Reached end')# }}}

# complex annotation# {{{
#this won't do cause if middle of the peak not in promoter range even if peak overlapping promoter I will not get it when I do the screening for peaks overlapping promoters
#complex.anno <- annotatePeakInBatch(complex$MergedPeaks, AnnotationData=as(my.annotation$genes, "RangedData"), PeakLocForDistance='middle', FeatureForDistance='TSS')

#complex.anno.2 <- annotatePeakInBatch(complex$MergedPeaks, AnnotationData=as(my.annotation$genes, "RangedData"), output='overlapping')

complex.anno <- addGeneIDs(as(complex.anno, "RangedData"), "org.Mm.eg.db", c("symbol", "ensembl") )
rownames(complex.anno) <- gsub(" ", ':',rownames(complex.anno))

#Convert from rangeddata to granges so you can add metadata
complex.anno <- as(complex.anno, "GRanges")
#0 means it doesn't bind the gene in the rowname
values(complex.anno)[,'score'] <- rep(0, nrow(values(complex.anno)))

#instead of end-start +500  could also use shortestDistance < 500 for peaks downstream

values(complex.anno)[!(complex.anno$insideFeature == 'upstream' & abs(complex.anno$distancetoFeature) > 2000)
                     & !(complex.anno$insideFeature == 'downstream' & complex.anno$distancetoFeature > (complex.anno$end_position - complex.anno$start_position) + 500), 'score' ] <- 1


ensembl <- useMart(host=host, biomart='ENSEMBL_MART_ENSEMBL', dataset=ds)
mapping <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "mgi_symbol"),  mart=ensembl)
#handle duplicated ensembl ids
mapping <- mapping[!duplicated(mapping$ensembl_gene_id),]

complex.anno.df <- as.data.frame(complex.anno)
#selecting for features/genes that have a peak at 2000 from the TSS and 500 from the TTS, or the peak resides in them
target.genes <- complex.anno.df[!(complex.anno.df$insideFeature == 'upstream' & abs(complex.anno.df$distancetoFeature) > 2000)
                                & !(complex.anno.df$insideFeature == 'downstream' & complex.anno.df$distancetoFeature > (complex.anno.df$end_position - complex.anno.df$start_position) + 500)
                                & !is.na(complex.anno.df$symbol), 'symbol']

#getting both genes if peak is overlapping 2 or more; ChIPpeakAnno reports these as separated by ';'
target.genes <- unlist(strsplit(as.vector(target.genes), ';'))
target.genes <- unique(target.genes)

#Save genes for ChEA analysis
write.table(target.genes,
            paste(paste(names(peaks), collapse='-'), unique(conf$condition), "target.genes", sep="."),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")


##### for annotation plots
complex.anno.df[,'Type'] <- factor(rep('NA', nrow(complex.anno.df)), levels=c('Body', 'Promoter', '3 prime', 'Enhancer', 'NA'))
complex.anno.df[,'Feature'] <- factor(rep('NA', nrow(complex.anno.df)), levels=c('Genic','Intergenic','NA'))


features <- list('Body'= with(complex.anno.df, (insideFeature=='overlapStart' & abs(distancetoFeature) > 500) | insideFeature== 'inside' | insideFeature== 'overlapEnd' | insideFeature== 'includeFeature'),
                 'Promoter'= with(complex.anno.df, (insideFeature=='upstream' & abs(distancetoFeature) < 2000) | (insideFeature== 'overlapStart' & shortestDistance <= 500)),
                 '3 prime'= with(complex.anno.df, (insideFeature== 'downstream' & distancetoFeature < (end_position - start_position + 500) )))

#Because when I iterate the list promoter peaks that also span the Body
for (type in names(features)){
   ind <- grep(TRUE, features[[type]])
   if( ! identical(ind, integer(0)) ) complex.anno.df[ind, 'Type'] <- factor(rep(type, length(ind)))
}

complex.anno.df[complex.anno.df$Type == 'NA', 'Feature'] <- 'Intergenic'
complex.anno.df[complex.anno.df$Type != 'NA', 'Feature'] <- 'Genic'

#Adding the genic, intergenic info to the grange so I can check for enhancer
values(complex.anno) <- cbind(values(complex.anno), complex.anno.df[, c('Type', 'Feature')])

load("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/enhancers.mm9.p300etal.Rdata")
enh <- reduce(enh)

values(complex.anno)[,'enhancer'] <- complex.anno %over% enh
complex.anno.df <- as.data.frame(complex.anno)
complex.anno.df[complex.anno.df$Feature == 'Intergenic' & complex.anno.df$enhancer== TRUE, 'Type'] <- 'Enhancer'

#check first if directory is there so you don't overwrite
dir.create(file.path(getwd(), 'plots'), showWarnings = FALSE)

pdf(paste0(paste(getwd(), 'plots/', sep='/'), paste(names(peaks),   collapse='-'), unique(conf$condition), '.GenFeat.pdf'))
p <- ggplot(complex.anno.df, aes(x=Feature, fill=Type))
p <- p + geom_bar(data=complex.anno.df, aes(y = (..count..)/sum(..count..)), width=.5) + scale_fill_grey(start = .4, end = .9, na.value = "grey50") + guides(fill= guide_legend(reverse=TRUE))
p <- p + ylab(paste0('Fraction of ', paste(names(peaks), collapse='-'), ' peaks in ', unique(conf$condition)))
p <- p + theme_bw() + theme(panel.grid.major.x = element_blank(),
                            panel.grid.major.y = element_blank(),
                            panel.background= element_blank(),
                            panel.border= element_blank(),
                            axis.line= element_line(),
                            axis.line.x= element_blank(),
                            axis.ticks= element_blank())
print(p)
dev.off()

#tss <- resize(my.annotation$genes, width=1, fix='start')
#tts <- resize(my.annotation$genes, width=1, fix='end')
#exons <- my.annotation$exons
#exonsPerGene <- split(my.annotation$exons, names(exons))
#firstExons <- lapply(names(exonsPerGene), function(gene){sort(exonsPerGene[[genei]][1])})
#firstexon<- endoapply(exonsPerGene, function(geneExons){sort(geneExons)[1]})
#save(firstexon, "/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/1stexon/Rdata")
# }}}
