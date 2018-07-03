# bsub -Ip /homes/ralser/source/R-3.2.0/bin/R
library (ChIPpeakAnno)
library(data.table)
library(Vennerable)
library(GenomicRanges)
library(gridExtra)
library(gplots)
require(VennDiagram)
library("RColorBrewer")

path = "/nfs/research2/bertone/user/ralser/CHIPseq/Hendrich/macs2_mm10/"
setwd(path)

## Define promoter and gene body positions
annotation = read.delim(file = "/nfs/research2/bertone/user/ralser/annotation/Ens75_genes+transscripts.txt",header = T)
genes = annotation[,c(1,10,9,18,19)]
genes = unique(genes)
annotation = annotation[,c(1,10,4,5,6,7,8,9)]
annotation = unique(annotation)
colnames(annotation) = c("gene_id","gene_name","chr","gene_start","gene_end","transcript_start","transcript_end","strand")
annotation$chr = paste0("chr",annotation$chr)
togrep = paste("chr",c(seq(1,19,1),"X","Y"),sep = "")
annotation = annotation[c(grep(paste(togrep,collapse = "|"),annotation$chr)),]

annotation$gene_id = as.character(annotation$gene_id)
annotation$promoter_start = 0
annotation$promoter_start = ifelse(annotation$strand==1,annotation$transcript_start-2001,annotation$transcript_end-501) ## bed file format 0-based
annotation$promoter_end = 0
annotation$promoter_end = ifelse(annotation$strand==1,annotation$transcript_start+500,annotation$transcript_end+2000)

## find first (if strand +) or last (if strand -) promoter positions
data.table(subset(annotation,strand==1))->dt.forward
dt.forward[,.SD[which.min(promoter_start)],by = gene_id][]->promoter_first
data.frame(promoter_first)->promoter_first
data.table(subset(annotation,strand!=1))->dt.reverse
dt.reverse[,.SD[which.max(promoter_start)],by = gene_id][]->promoter_last
data.frame(promoter_last)->promoter_last
annotation_promoter = rbind(promoter_first,promoter_last)

## find gene body coordinates (bed format 0-based)
annotation_promoter$gene_body_start = 0
annotation_promoter$gene_body_start = ifelse(annotation_promoter$strand==1,annotation_promoter$gene_start+501,annotation_promoter$gene_start)
annotation_promoter$gene_body_end = 0
annotation_promoter$gene_body_end = ifelse(annotation_promoter$strand==1,annotation_promoter$gene_end,annotation_promoter$gene_end-502)

## gene body and promoter
promcoord = data.frame(annotation_promoter[,c("chr","promoter_start","promoter_end")],paste(annotation_promoter$gene_id,"prom",sep = ":"))
promcoord = unique(promcoord)
colnames(promcoord) = c("chr","start","end","id")
bodycoord = data.frame(annotation_promoter[,c("chr","gene_body_start","gene_body_end")],paste(annotation_promoter$gene_id,"gene_body",sep = ":"))
bodycoord = unique(bodycoord)
colnames(bodycoord) = c("chr","start","end","id")
## remove very short genes 500bp (gene body positions) 
bodycoord = subset(subset(bodycoord,start<end))

coord = rbind(bodycoord,promcoord)
coord = coord[ do.call(order, coord), ]
coord_Ranged = BED2RangedData(coord,header = F)
coord_GRanges = toGRanges(coord_Ranged, format = "RangedData")



annotate_peaks <- function(rangedData,name) {
        annotatedSample = annotatePeakInBatch(rangedData, AnnotationData = coord_Ranged,output = "overlapping")
        annotatedSample = as.data.frame(annotatedSample)
        annotatedSample = annotatedSample[,c("space","start","end","width","names", "peak","strand","feature","start_position","end_position","insideFeature","distancetoFeature","shortestDistance","fromOverlappingOrNearest")]
        annotatedSample = subset(annotatedSample,!is.na(insideFeature))
        annotatedSample = cbind(annotatedSample, read.table(text = as.character(annotatedSample$feature), sep = ":"))
        colnames(annotatedSample)[c(15,16)] = c("gene_id","status")

        ## find interval
        annotatedSampleinterval_start = 0
        annotatedSample$interval_start = ifelse(annotatedSample$start>annotatedSample$start_position,annotatedSample$start,annotatedSample$start_position)
        annotatedSample$interval_end = 0
        annotatedSample$interval_end = ifelse(annotatedSample$end<annotatedSample$end_position,annotatedSample$end,annotatedSample$end_position)

        ## if a peak overlaps with the prom and gene body of a gene, select the prom
        annotatedSample$peak_gene_id = paste(annotatedSample$peak,annotatedSample$gene_id,sep = " ")
        index = annotatedSample[duplicated(annotatedSample[,c("peak_gene_id")]),"peak_gene_id"]
        annotatedSampleDuplicated = annotatedSample[annotatedSample$peak_gene_id %in% index,]
        annotatedSamplenotDuplicated = annotatedSample[!annotatedSample$peak_gene_id %in% index,]
        annotatedSample = rbind(annotatedSamplenotDuplicated,subset(annotatedSampleDuplicated,status=="prom"))

        ## a peak could be associated to two genes
        annotatedSample = annotatedSample[,-c(19)]
        assign(paste0(name,"_annotated"),annotatedSample,envir = .GlobalEnv)
        assign(paste0(name,"_gene_body"),subset(annotatedSample,status=="gene_body"),envir = .GlobalEnv)
        assign(paste0(name,"_prom"),subset(annotatedSample,status=="prom"),envir = .GlobalEnv)
}

## Define NuRD bound genes
## 2i condition
Mbd3 = read.table("7E12_M2_2i_peaks.narrowPeak")
Mbd3 = Mbd3[,c(1:4)]
Mbd3_Ranged = BED2RangedData(Mbd3,header = F)
Mbd3_GRanges = toGRanges(Mbd3_Ranged, format = "RangedData")

Chd4 = read.table("7E12_Chd4_2i_peaks.narrowPeak")
Chd4 = Chd4[,c(1:4)]
Chd4_Ranged = BED2RangedData(Chd4,header = F)
Chd4_GRanges = toGRanges(Chd4_Ranged, format = "RangedData")

overlap = findOverlapsOfPeaks(Mbd3_GRanges,Chd4_GRanges, maxgap = 0)

Chd4_peaks = overlap$peaklist[["Chd4_GRanges"]]
Chd4_peaks = as.data.frame(unlist(Chd4_peaks))
Chd4_peaks$peaks = as.character(Chd4_peaks$peakNames)
Chd4_peaks = Chd4_peaks[,-c(4,5,6)]
Chd4_peaks$seqnames = paste0("chr",Chd4_peaks$seqnames)
Chd4_peaks$peaks = gsub("Chd4_GRanges__","",Chd4_peaks$peaks)
Chd4_peaks$peaks = gsub("X","",Chd4_peaks$peaks)
Chd4_peaks$peaks = gsub("\\.","-",Chd4_peaks$peaks)


Mbd3 = overlap$peaklist[["Mbd3_GRanges"]]
Mbd3_peaks = as.data.frame(unlist(Mbd3))
Mbd3_peaks$peaks = as.character(Mbd3_peaks$peakNames)
Mbd3_peaks = Mbd3_peaks[,-c(4,5,6)]
Mbd3_peaks$seqnames = paste0("chr",Mbd3_peaks$seqnames)
Mbd3_peaks$peaks = gsub("Mbd3_GRanges__","",Mbd3_peaks$peaks)
Mbd3_peaks$peaks = gsub("X","",Mbd3_peaks$peaks)
Mbd3_peaks$peaks = gsub("\\.","-",Mbd3_peaks$peaks)

mergedPeaks = overlap$peaklist[[3]]
mergedPeaks = as.data.frame(unlist(mergedPeaks))
mergedPeaks = mergedPeaks[,-c(4:6)]
mergedPeaks$seqnames = paste0("chr",mergedPeaks$seqnames)


NuRD_positions = mergedPeaks
NuRD_positions = data.frame(mergedPeaks,ID = paste0("NuRD_",c(1:nrow(mergedPeaks))))
NuRD_Ranged = BED2RangedData(NuRD_positions,header = F)
NuRD_GRanges = toGRanges(NuRD_Ranged, format = "RangedData")

annotate_peaks(NuRD_Ranged,"NuRD")
NuRD_2i_annotated = get("NuRD_annotated")
NuRD_2i_annotated =subset(NuRD_2i_annotated,status=="prom")
NuRD_2i_genes=data.frame(gene_id=NuRD_2i_annotated[,"gene_id"])
NuRD_2i_genes=unique(NuRD_2i_genes)

## SL condition
Mbd3 = read.table("7E12_Mbd3_SL_peaks.narrowPeak")
Mbd3 = Mbd3[,c(1:4)]
Mbd3_Ranged = BED2RangedData(Mbd3,header = F)
Mbd3_GRanges = toGRanges(Mbd3_Ranged, format = "RangedData")

Chd4 = read.table("7E12_Chd4_SL_peaks.narrowPeak")
Chd4 = Chd4[,c(1:4)]
Chd4_Ranged = BED2RangedData(Chd4,header = F)
Chd4_GRanges = toGRanges(Chd4_Ranged, format = "RangedData")


overlap = findOverlapsOfPeaks(Mbd3_GRanges,Chd4_GRanges, maxgap = 0)#,connectedPeaks="min")

Chd4_peaks = overlap$peaklist[["Chd4_GRanges"]]
Chd4_peaks = as.data.frame(unlist(Chd4_peaks))
Chd4_peaks$peaks = as.character(Chd4_peaks$peakNames)
Chd4_peaks = Chd4_peaks[,-c(4,5,6)]
Chd4_peaks$seqnames = paste0("chr",Chd4_peaks$seqnames)
Chd4_peaks$peaks = gsub("Chd4_GRanges__","",Chd4_peaks$peaks)
Chd4_peaks$peaks = gsub("X","",Chd4_peaks$peaks)
Chd4_peaks$peaks = gsub("\\.","-",Chd4_peaks$peaks)


Mbd3 = overlap$peaklist[["Mbd3_GRanges"]]
Mbd3_peaks = as.data.frame(unlist(Mbd3))
Mbd3_peaks$peaks = as.character(Mbd3$peakNames)
Mbd3_peaks = Mbd3_peaks[,-c(4,5,6)]
Mbd3_peaks$seqnames = paste0("chr",Mbd3_peaks$seqnames)
Mbd3_peaks$peaks = gsub("Mbd3_GRanges__","",Mbd3_peaks$peaks)
Mbd3_peaks$peaks = gsub("X","",Mbd3_peaks$peaks)
Mbd3_peaks$peaks = gsub("\\.","-",Mbd3_peaks$peaks)
mergedPeaks = overlap$peaklist[[3]]
mergedPeaks = as.data.frame(unlist(mergedPeaks))
mergedPeaks = mergedPeaks[,-c(4:6)]
mergedPeaks$seqnames = paste0("chr",mergedPeaks$seqnames)

NuRD_positions = mergedPeaks
NuRD_positions = data.frame(mergedPeaks,ID = paste0("NuRD_",c(1:nrow(mergedPeaks))))
NuRD_Ranged = BED2RangedData(NuRD_positions,header = F)
NuRD_GRanges = toGRanges(NuRD_Ranged, format = "RangedData")

annotate_peaks(NuRD_Ranged,"NuRD")
NuRD_SL_annotated = get("NuRD_annotated")
NuRD_SL_annotated =subset(NuRD_SL_annotated,status=="prom")
NuRD_SL_genes=data.frame(gene_id=NuRD_SL_annotated[,"gene_id"])
NuRD_SL_genes=unique(NuRD_SL_genes)


## calculate gene length
gene_annotation=annotation[,c("gene_id","gene_name","chr","strand","gene_start","gene_end")]
gene_annotation=unique(gene_annotation)
gene_annotation$gene_length=(gene_annotation$gene_end-gene_annotation$gene_start)+1
genes_long=subset(gene_annotation,gene_length>=800)


## add expression
expression = read.table(file="/nfs/research2/bertone/user/ralser/RNAseq/Brian/expression/20150811_corr_Expression_SL_and_2i",header=T,sep="\t")
deseq_2i = read.table(file="/nfs/research2/bertone/user/ralser/RNAseq/Brian/expression/20150811_2i_DESeq2",header=T,sep="\t")
colnames(deseq_2i)[9:14]=paste0(colnames(deseq_2i)[9:14],"_2i")
deseq_SL = read.table(file="/nfs/research2/bertone/user/ralser/RNAseq/Brian/expression/20150811_SL_DESeq2",header=T,sep="\t")
colnames(deseq_SL)[9:14]=paste0(colnames(deseq_SL)[9:14],"_SL")
expression = merge(expression,deseq_2i)
expression = merge(expression,deseq_SL)

expression$Status_2i=""
expression$Status_2i = ifelse(expression$padj_2i<=0.05 & expression$log2FoldChange_2i<0, "UP in KO",expression$Status_2i)
expression$Status_2i = ifelse(expression$padj_2i<=0.05 & expression$log2FoldChange_2i>0, "DOWN in KO",expression$Status_2i)

expression$Status_SL=""
expression$Status_SL = ifelse(expression$padj_SL<=0.05 & expression$log2FoldChange_SL<0, "UP in KO",expression$Status_SL)
expression$Status_SL = ifelse(expression$padj_SL<=0.05 & expression$log2FoldChange_SL>0, "DOWN in KO",expression$Status_SL)

## read gene signal tables

create.table <- function(org) {
  files <- dir(pattern=glob2rx(paste("*",org,"*",".txt",sep="")))
  names<-sub("_gene_signal.txt","",files)
  count.list <- lapply(files, read.table, header=T)
  colnames(count.list[[1]])[3:4] <- paste(names[1],".",colnames(count.list[[1]])[3:4],sep="")
  count.matrix <-  count.list[[1]]
  if(length(count.list) == 1) return(count.matrix)
  for(i in 2:length(count.list)) {
    colnames(count.list[[i]])[3:4] <- paste(names[i],".",colnames(count.list[[i]])[3:4],sep="")
    count.matrix <- cbind(count.matrix, count.list[[i]][match(count.list[[1]]$gene.id,count.list[[i]]$gene.id),3:4])
  }
  rownames(count.matrix) <- count.list[[1]]$gene.id
  count.matrix
}

setwd("/nfs/research2/bertone/user/ralser/CHIPseq/Hendrich/gene_signal_mm10")
vals <- create.table("_gene_signal")
colnames(vals)[1:2]=c("gene_id","gene_name")
vals <- vals[vals$gene_id %in% expression$gene_id,]


##### ALL comparision
setwd("/nfs/research2/bertone/user/ralser/CHIPseq/Hendrich/polII_mm10/plots/")

pausing.index.sample <- function(data,sample) {
  prom=data[,grep(paste(sample,".prom.max",sep=""),colnames(data))]
  body=data[,grep(paste(sample,".body.max",sep=""),colnames(data))]
  totals=prom+body		
  index <- round(prom/(totals+0.001),digits=2)
  index[totals < 20] <- NA
  return(index)
}


######## calculate travelling ratio per sample 
samples=unique(unlist(lapply(strsplit(colnames(vals[,-c(1:2)]),"\\."), function(x) x[1])))

pi.table=data.frame(vals[,1:2])
for(t in 1:length(samples))
{
	print(t)
	pis=pausing.index.sample(vals,as.character(samples[t]))
	pi.table=cbind(pi.table,pis)
	colnames(pi.table)[ncol(pi.table)]=samples[t]
}

colnames(pi.table)=gsub("-","_",colnames(pi.table))

## generate PI density and cumulative plots
pi.plot=function(data,samples,conditions)
{
	comparison=data[,c(1,2)]
	for(s in 1:length(samples))
	{
		comparison=merge(comparison,data[,c(1,2,grep(paste("^",samples[s],"$",sep=""),colnames(data)))])
		
	}
	    comparison.long=merge(comparison,genes_long[,c("gene_id","gene_name")])
	
    comparison.expression=merge(comparison,expression[,c("gene_id","gene_name","padj_SL","padj_2i","Status_SL","Status_2i")])
		
	
        for(s in 1:length(conditions))
        {
                assign(paste0("comparison.nurd.condition",s),merge(data[,c(1,2,grep(paste("^",samples[s],"$",sep=""),colnames(data)))],get(paste0("NuRD_",conditions[s],"_genes"))))
		assign(paste0("comparison.nurd.long.condition",s),merge(get(paste0("comparison.nurd.condition",s)),genes_long[,c("gene_id","gene_name")]))

		## de genes
 	        de=comparison.expression[,c(1,2,grep(paste("^",samples[s],"$",sep=""),colnames(comparison.expression)),grep(paste0("padj_",conditions[s]),colnames(comparison.expression)),grep(paste0("Status_",conditions[s]),colnames(comparison.expression)))]
		de=subset(de,de[,4]<=0.05)
		assign(paste0("comparison.de.condition",s),de)
		assign(paste0("comparison.down.condition",s),subset(de,de[,5]=="DOWN in KO"))
		assign(paste0("comparison.up.condition",s),subset(de,de[,5]=="UP in KO"))
		## de nurd genes

		de.nurd=merge(de,get(paste0("NuRD_",conditions[s],"_genes")))
		assign(paste0("comparison.de.nurd.condition",s),de.nurd)
		assign(paste0("comparison.down.nurd.condition",s),subset(de.nurd,de.nurd[,5]=="DOWN in KO"))
                assign(paste0("comparison.up.nurd.condition",s),subset(de.nurd,de.nurd[,5]=="UP in KO"))

        }
		
		## add de genes
		cols=c("blue","red","green","magenta")
		name=paste(samples, collapse = "_vs_")
                
		# Density plots
		## all
                pdf(paste("PI_density_",name,".pdf",sep=""))
                plot(density(comparison[,3],na.rm = TRUE),col=cols[1],main="Pausing Index",ylim=c(0,3))
		for(s in 4:ncol(comparison))
		{
                	lines(density(comparison[,s],na.rm = TRUE),col=cols[s-2])
		}
                legend("topleft",legend=colnames(comparison)[c(3:ncol(comparison))], fill= cols[1:length(samples)], bty="n")
                dev.off()

		## long genes
                pdf(paste("PI_density_",name,"_Long.pdf",sep=""))
                plot(density(comparison.long[,3],na.rm = TRUE),col=cols[1],main="Pausing Index",ylim=c(0,3))
                for(s in 4:ncol(comparison))
                {
                        lines(density(comparison.long[,s],na.rm = TRUE),col=cols[s-2])
                }
                legend("topleft",legend=colnames(comparison.long)[c(3:ncol(comparison.long))], fill= cols[1:length(samples)], bty="n")
                dev.off()

		## NuRD bound
		pdf(paste("PI_density_",name,"_NuRD_bound.pdf",sep=""))
                plot(density(comparison.nurd.condition1[,3],na.rm = TRUE),col=cols[1],main="Pausing Index NuRD bound",ylim=c(0,3))
		for(s in 2:length(samples))
		{
                	lines(density(get(paste0("comparison.nurd.condition",s))[,3],na.rm = TRUE),col=cols[s])
		}
                legend("topleft",legend=samples, fill=cols[1:length(samples)], bty="n")
                dev.off()

		## NuRD bound long genes
                pdf(paste("PI_density_",name,"_NuRD_bound_Long.pdf",sep=""))
                plot(density(comparison.nurd.long.condition1[,3],na.rm = TRUE),col=cols[1],main="Pausing Index NuRD bound",ylim=c(0,3))
                for(s in 2:length(samples))
                {
                        lines(density(get(paste0("comparison.nurd.long.condition",s))[,3],na.rm = TRUE),col=cols[s])
                }
                legend("topleft",legend=samples, fill=cols[1:length(samples)], bty="n")
                dev.off()
		
	        ## de genes
                pdf(paste("PI_density_",name,"_DE.pdf",sep=""))
                plot(density(comparison.de.condition1[,3],na.rm = TRUE),col=cols[1],main="Pausing Index DE genes",ylim=c(0,3))
                for(s in 2:length(samples))
                {
                        lines(density(get(paste0("comparison.de.condition",s))[,3],na.rm = TRUE),col=cols[s])
                }
                legend("topleft",legend=samples, fill=cols[1:length(samples)], bty="n")
                dev.off()

                ## down in ko
                pdf(paste("PI_density_",name,"_DOWN-in-KO.pdf",sep=""))
                plot(density(comparison.down.condition1[,3],na.rm = TRUE),col=cols[1],main="Pausing Index Down in KO",ylim=c(0,3))
                for(s in 2:length(samples))
                {
                        lines(density(get(paste0("comparison.down.condition",s))[,3],na.rm = TRUE),col=cols[s])
                }
                legend("topleft",legend=samples ,fill=cols[1:length(samples)], bty="n")
                dev.off()

                ## up in ko
                pdf(paste("PI_density_",name,"_UP-in-KO.pdf",sep=""))
                plot(density(comparison.up.condition1[,3],na.rm = TRUE),col=cols[1],main="Pausing Index Up in KO",ylim=c(0,3))
                for(s in 2:length(samples))
                {
                        lines(density(get(paste0("comparison.up.condition",s))[,3],na.rm = TRUE),col=cols[s])
                }
                legend("topleft",legend=samples, fill=cols[1:length(samples)], bty="n")
                dev.off()

                ## de nurd bound
                pdf(paste("PI_density_",name,"_DE_and_NuRD_bound.pdf",sep=""))
                plot(density(comparison.de.nurd.condition1[,3],na.rm = TRUE),col=cols[1],main="Pausing Index DE NuRD bound",ylim=c(0,3))
                for(s in 2:length(samples))
                {
                        lines(density(get(paste0("comparison.de.nurd.condition",s))[,3],na.rm = TRUE),col=cols[s])
                }
                legend("topleft",legend=samples ,fill=cols[1:length(samples)], bty="n")
                dev.off()

                ## down in ko and nurd bound
                pdf(paste("PI_density_",name,"_DOWN-in-KO_and_NuRD_bound.pdf",sep=""))
                plot(density(comparison.down.nurd.condition1[,3],na.rm = TRUE),col=cols[1],main="Pausing Index Down in KO NuRD bound",ylim=c(0,3))
                for(s in 2:length(samples))
                {
                        lines(density(get(paste0("comparison.down.nurd.condition",s))[,3],na.rm = TRUE),col=cols[s])
                }
                legend("topleft",legend=samples ,fill=cols[1:length(samples)], bty="n")
                dev.off()

                ## up in ko and nurd bound
                pdf(paste("PI_density_",name,"_UP-in-KO_and_NuRD_bound.pdf",sep=""))
                plot(density(comparison.up.nurd.condition1[,3],na.rm = TRUE),col=cols[1],main="Pausing Index Up in KO NuRD bound",ylim=c(0,3))
                for(s in 2:length(samples))
                {
                        lines(density(get(paste0("comparison.up.nurd.condition",s))[,3],na.rm = TRUE),col=cols[s])
                }
                legend("topleft",legend=samples ,fill=cols[1:length(samples)], bty="n")
                dev.off()


		# Cumulative density                    
		## all 
                pdf(paste("PI_cumulative_density_",name,".pdf",sep=""))
                plot(ecdf(comparison[!is.na(comparison[,3]),3]),col=cols[1],main="Pausing Index",xlab="Pausing index",ylab="Cumulative distribution")
		for(s in 4:ncol(comparison))
		{
                	lines(ecdf(comparison[!is.na(comparison[,s]),s]),col=cols[s-2])
		}
		legend("topleft",legend=colnames(comparison)[c(3:ncol(comparison))], fill= cols[1:length(samples)], bty="n")
                dev.off()

		## long genes
                pdf(paste("PI_cumulative_density_",name,"_Long.pdf",sep=""))
                plot(ecdf(comparison.long[!is.na(comparison.long[,3]),3]),col=cols[1],main="Pausing Index",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 4:ncol(comparison))
                {
                        lines(ecdf(comparison.long[!is.na(comparison.long[,s]),s]),col=cols[s-2])
                }
                legend("topleft",legend=colnames(comparison.long)[c(3:ncol(comparison.long))], fill= cols[1:length(samples)], bty="n")
                dev.off()

		## NuRD bound
                pdf(paste("PI_cumulative_density_",name,"_NuRD_bound.pdf",sep=""))
                plot(ecdf(comparison.nurd.condition1[!is.na(comparison.nurd.condition1[,3]),3]),col=cols[1],main="Pausing Index NuRD bound",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 2:length(samples))
                {
			tmp=get(paste0("comparison.nurd.condition",s))
			lines(ecdf(tmp[!is.na(tmp[,3]),3]),col=cols[s])
                }
                legend("topleft",legend=samples, fill=cols[1:length(samples)], bty="n")
                dev.off()

		## NuRD bound long
                pdf(paste("PI_cumulative_density_",name,"_NuRD_bound_Long.pdf",sep=""))
                plot(ecdf(comparison.nurd.long.condition1[!is.na(comparison.nurd.long.condition1[,3]),3]),col=cols[1],main="Pausing Index NuRD bound",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 2:length(samples))
                {
                        tmp=get(paste0("comparison.nurd.long.condition",s))
                        lines(ecdf(tmp[!is.na(tmp[,3]),3]),col=cols[s])
                }
                legend("topleft",legend=samples ,fill=cols[1:length(samples)], bty="n")
                dev.off()
		
		## de genes
                pdf(paste("PI_cumulative_density_",name,"_DE.pdf",sep=""))
                plot(ecdf(comparison.de.condition1[!is.na(comparison.de.condition1[,3]),3]),col=cols[1],main="Pausing Index DE genes",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 2:length(samples))
                {
                        tmp=get(paste0("comparison.de.condition",s))
                        lines(ecdf(tmp[!is.na(tmp[,3]),3]),col=cols[s])
                }
                legend("topleft",legend=samples, fill=cols[1:length(samples)], bty="n")
                dev.off()

		## down in ko
                pdf(paste("PI_cumulative_density_",name,"_DOWN-in-KO.pdf",sep=""))
                plot(ecdf(comparison.down.condition1[!is.na(comparison.down.condition1[,3]),3]),col=cols[1],main="Pausing Index Down in KO",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 2:length(samples))
                {
                        tmp=get(paste0("comparison.down.condition",s))
                        lines(ecdf(tmp[!is.na(tmp[,3]),3]),col=cols[s])
                }
                legend("topleft",legend=samples, fill=cols[1:length(samples)], bty="n")
                dev.off()

		## up in ko
                pdf(paste("PI_cumulative_density_",name,"_UP-in-KO.pdf",sep=""))
                plot(ecdf(comparison.up.condition1[!is.na(comparison.up.condition1[,3]),3]),col=cols[1],main="Pausing Index Up in KO",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 2:length(samples))
                {
                        tmp=get(paste0("comparison.up.condition",s))
                        lines(ecdf(tmp[!is.na(tmp[,3]),3]),col=cols[s])
                }
                legend("topleft",legend=samples ,fill=cols[1:length(samples)], bty="n")
                dev.off()

                ## de and nurd bound
                pdf(paste("PI_cumulative_density_",name,"_DE_and_NuRD_bound.pdf",sep=""))
                plot(ecdf(comparison.de.nurd.condition1[!is.na(comparison.de.nurd.condition1[,3]),3]),col=cols[1],main="Pausing Index DE NuRD bound",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 2:length(samples))
                {
                        tmp=get(paste0("comparison.de.nurd.condition",s))
                        lines(ecdf(tmp[!is.na(tmp[,3]),3]),col=cols[s])
                }
                legend("topleft",legend=samples, fill=cols[1:length(samples)], bty="n")
                dev.off()
		
                ## down in ko and nurd bound
                pdf(paste("PI_cumulative_density_",name,"_DOWN-in-KO_and_NuRD_bound.pdf",sep=""))
                plot(ecdf(comparison.down.nurd.condition1[!is.na(comparison.down.nurd.condition1[,3]),3]),col=cols[1],main="Pausing Index Down in KO NuRD bound",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 2:length(samples))
                {
                        tmp=get(paste0("comparison.down.nurd.condition",s))
                        lines(ecdf(tmp[!is.na(tmp[,3]),3]),col=cols[s])
                }
                legend("topleft",legend=samples ,fill=cols[1:length(samples)], bty="n")
                dev.off()

                ## up in ko and nurd bound
                pdf(paste("PI_cumulative_density_",name,"_UP-in-KO_and_NuRD_bound.pdf",sep=""))
                plot(ecdf(comparison.up.nurd.condition1[!is.na(comparison.up.nurd.condition1[,3]),3]),col=cols[1],main="Pausing Index Up in KO NuRD bound",xlab="Pausing index",ylab="Cumulative distribution")
                for(s in 2:length(samples))
                {
                        tmp=get(paste0("comparison.up.nurd.condition",s))
                        lines(ecdf(tmp[!is.na(tmp[,3]),3]),col=cols[s])
                }
                legend("topleft",legend=samples ,fill=cols[1:length(samples)], bty="n")
                dev.off()

}




 # samples
 
#"2lox_S5P"     "2lox_S5P_2"   "2lox_8WG"    
#"2lox2i_8WG"   "2lox2i_S5P"   "7E12_2P"      "7E12_8WG"     "7E12_S5P"
# "7E122i_8WG"   "7E122i_8WG_2" "7E122i_S5P"   "7E122i_S5P_2" "7G9_2P"      
# "7G9_8WG"      "7G9_S5P"      "7G92i_S5P"    "7G92i_S5Pa"   "7G92i_S5Pb"  
# "C9_8WG"       "Spl_2i_8WG"   "Spl_2i_S5P"   "Spl_8WG"      "Spl_S5P"     


pi.plot(pi.table,c("2lox_S5P","Spl_S5P"),c("SL","SL"))
pi.plot(pi.table,c("2lox_S5P_2","Spl_S5P"),c("SL","SL"))
pi.plot(pi.table,c("2lox_S5P","2lox2i_S5P"),c("SL","2i"))
pi.plot(pi.table,c("2lox_8WG","2lox2i_8WG"),c("SL","2i"))
pi.plot(pi.table,c("7E12_S5P","7G9_S5P","7E122i_S5P","7G92i_S5P"),c("SL","SL","2i","2i"))



