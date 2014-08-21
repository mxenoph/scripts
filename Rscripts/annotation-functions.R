
annotationInfo <- function (assembly){
    conf <- read.table("/nfs/research2/bertone/user/mxenoph/genome_dir/assemblies-annotations.config", header=T, stringsAsFactors=F, sep="\t")
    chr_size <<- conf[conf$assembly == assembly, 'chrom_sizes']
    annotation <<- conf[conf$assembly == assembly, 'annotation']
 }
