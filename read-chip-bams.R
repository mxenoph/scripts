#!/homes/mxenoph/local/bin/Rscript --vanilla

#Turning off warning messages for loading the packages-globally
options(warn=-1)
suppressMessages(require(Repitools))
#Turn warnings back on
options(warn=0)

args <- commandArgs(trailingOnly=TRUE)
#Directory to search for the chip and input bam files
bam.dir <- args[1]
prefix <- args[2]

paths <- list.files("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/bowtie", pattern="*E12_2i.*sort.bam", full.name=T)
paths <- c(paths, list.files("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/bowtie", pattern="*E12.*Input.bam", full.name=T))
grL <- BAM2GRangesList(paths)
names(grL) <- paths

#save(grL, file= "/nfs/nobackup2/research/bertone/mxenoph/bh/ChIP-Epi/mm9/bowtie/bam2grList.Rdata")
save(grL, file= "/nfs/nobackup2/research/bertone/mxenoph/bh/ChIP-Epi/mm9/bowtie/2ibam2grList.Rdata")


#paths <- list.files("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/bowtie", pattern="*E12_2i.*sort.bam", full.name=T)

#grL <- BAM2GRangesList(paths)
#names(grL) <- paths
#save(grL, file= "/nfs/nobackup2/research/bertone/mxenoph/bh/ChIP-Epi/mm9/bowtie/2ibam2grList.Rdata")
