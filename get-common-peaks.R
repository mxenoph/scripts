#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")

parser <-  ArgumentParser(description="Find common peaks for TFs")
parser$add_argument('-f', '--first', metavar= "file", required='True', type= "character", help= "Config file with (system, factor, genotype, condition, replicate-index,file) columns for complex/TF. e.g. chip/config/oct4.conf")
parser$add_argument('-s', '--second', metavar= "file", required='True', type= "character", help= "Config file with (system, factor, genotype, condition, replicate-index,file) columns for complex/TF. e.g. chip/config/oct4.conf")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-l', '--label', type= "character", default=format(Sys.time(), "%d%b%Y"), help= "Give a label for the common peaks examined")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args <- parser$parse_args()
output_path <- file.path(args$out, 'targets', args$label)
plot_path <- file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)

# function in annotation-functions.R
get_annotation(tolower(args$assembly))
#}}}

# Packages# {{{
# List of packages to load
x <- c('GenomicRanges', 'rtracklayer', 'pheatmap', 'RColorBrewer', 'foreach', 'parallel', 'bigmemory', 'biganalytics', 'NMF')
lapply(x, suppressMessages(library), character.only=T)
# }}}

intersect_multi <- function(grl){# {{{
    gr <- grl[[1]]
    for (j in seq(2, length(grl), 1)){
        gr <- intersect(gr, grl[[j]])

    }
    names(gr) <- name_gr(gr)
    return(gr)
}# }}}

get_hits_scores <- function(query, subj, value) {# {{{
    # returns a vector of length(common) containing index for hits in x
    ov <- findOverlaps(query, subj, select='first')
    names(ov) <- names(query)
    # subset based on those that have a hit
    ov[!is.na(ov)] <- values(subj[as.numeric(ov[!is.na(ov)])])[[value]]
    return(ov)
}# }}}

parse_peaks <- function(config){# {{{
    source("~/source/Rscripts/granges-functions.R")
    # This assumes it's mouse data. If not this has to be changed
    chr <- c(paste0('chr', 1:19), 'chrX', 'chrY')
    
    y <- mclapply(as.vector(config[,'file']), function(files){
                  #bed <- import.bed(files)
                  bed <- macs_to_granges(files)
                  bed <- subset(bed, seqnames(bed) %in% chr)
                  names(bed) <- name_gr(bed)
                  return(bed)
             })
    names(y) <- with(config, paste(system, factor, genotype, condition, replicates, rownames(config), sep=':'))
    
    descriptors <- unique(gsub(":[0-9]{1,2}$", '', names(y)))
    pooled <- sapply(descriptors, function(x){
                     i <- grep(x, names(y))
                     if(length(i) == 1) {
                         return(y[[i]])
                     }else {
                         gr <- intersect_multi(y[i])
                         for(j in colnames(values(y[[1]])) ) values(gr)[j] <- rep(NA, length(gr))
                         
                         scores <- mclapply(y[i], function(x) get_hits_scores(gr, x, 'FE'))
                         names(scores) <- names(y[i])
                         scores <- do.call(cbind.data.frame, scores)
                         scores[,'mean'] <- apply(scores,1,mean)
                         values(gr)['FE'] <- scores[['mean']]
                         #apply(scores,1,sd)
                         return(gr)
                     }
             })
    return(list(y, pooled))
}# }}}

plotting <- function(){# {{{
    cool_cols <- colorRampPalette(c('aliceblue','darkcyan'))(100)
    heat_cols <- colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(50)
    pdf(file.path(plot_path, paste0(args$label, '.pdf')), paper='a4')

    aheatmap(pdata[,-ncol(pdata)], scale="row", Rowv=NA, Colv=NA,
             col=heat_cols, border_color=NA, annRow=side_anno)

    # Performs kmeans clustering and aggregates the rows in each cluster
    pheatmap(mat, #trace= "none", density.info= "none",
             color=cool_cols, kmeans_k=50,
             clustering_distance_cols= "binary", clustering_method= "ward",
             cluster_cols=TRUE, cluster_rows=TRUE,
             border_color=NA, #annotation_colors= , annotation_legend=F,
             fontsize=10, fontsize_row=7, show_rownames=FALSE, show_colnames=TRUE)

    dev.off()
}# }}}

bigcor <- function(x, nblocks = 10, verbose = TRUE, ...)# {{{
{
    library(ff, quietly = TRUE)
    library(doMC)

    #detecting how many cores on the machine and setting the limit
    if(ncore=="all"){
        ncore = parallel:::detectCores()
        registerDoMC(cores = ncore)
    } else{
        registerDoMC(cores = ncore)
    }

    NCOL <- ncol(x)
     
    ## test if ncol(x) %% nblocks gives remainder 0
    if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
     
    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
    corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
     
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
     
    ## create all unique combinations of blocks
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)
     
    ## iterate through each block combination, calculate correlation matrix
    ## between blocks and store them in the preallocated matrix on both
    ## symmetric sides of the diagonal
    for (i in 1:nrow(COMBS)) {
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]
        if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()
        COR <- cor(x[, G1], x[, G2])#, ...)
        corMAT[G1, G2] <- COR
        corMAT[G2, G1] <- t(COR)
        COR <- NULL
    }
     
    size <- floor(nrow(corMAT)/nblocks)
    residual <- nrow(corMAT) %% nblocks
    
    distMAT <- ff(vmode = "single", dim = dim(corMAT))
    pointer <- 1
    for(i in 1:nblocks){
        current_range <- pointer:(pointer+size-1)
        distMAT[current_range, current_range] <- (1-corMAT[current_range, current_range])/2
        pointer <- pointer + size
    }
    if(residual !=0) {
        current_range <- pointer:(pointer+residual-1)
        distMAT[current_range, current_range] <- (1-corMAT[current_range, current_range])/2
    }

    gc()
    return(corMAT)
}
# }}}

main <- function(){# {{{

    chr <- c(paste0('chr', 1:19), 'chrX', 'chrY')
    # get protein coding genes
    genes <- import(genes_gtf, asRangedData=F)
    genes_all <- import(genes_gtf, asRangedData=F)
    genes <- subset(genes, seqnames(genes) %in% chr & source == 'protein_coding')

    first <- read.table(args$first, header=T, sep="\t")
    second <- read.table(args$second, header=T, sep="\t")

    name_gr <- function(gr) paste(seqnames(gr), start(gr), end(gr), sep=':')
    get_score <- function(gr, value) values(gr)[[value]]
    x <- parse_peaks(first)
    pooled_x <- GRangesList(x[[2]])
    y <- parse_peaks(second)
    pooled_y <- GRangesList(y[[2]])

    x_int <- intersect_multi(pooled_x)
    names(x_int) <- name_gr(x_int)
    y_int <- intersect_multi(pooled_y)
    names(y_int) <- name_gr(y_int)

    factor_peaks <- function(pooled){# {{{
        x  <- reduce(unlist(pooled))
        names(x) <- name_gr(x)

        values(x) <- cbind(values(x),
                           as.data.frame(mclapply(pooled,
                                                  function(p) get_hits_scores(x, p,'FE'))))

        # have to parse to data.frame cause granges 
        #dataframe not exactly the same as r dataframe I think
        meta <- as.data.frame(values(x))
        meta[is.na(meta)] <- 0
        values(x) <- meta
        return(x)
    }# }}}

    x_all <- factor_peaks(pooled_x)
    y_all <- factor_peaks(pooled_y)

    diagnostics <- function(meta){# {{{
        meta <- as.data.frame(meta)
        #meta <- as.data.frame(t(scale(t(as.matrix(meta)))))

        meta[,'mean'] <- rowMeans(meta)
        meta[,'mad'] <- apply(as.matrix(meta), 1, mad)
        meta[,'sd'] <- apply(as.matrix(meta), 1, sd)
        meta <- meta[order(-meta[,'mad']),]

        return(meta)
    }# }}}

    x_meta <- diagnostics(values(x_all))
    y_meta <- diagnostics(values(y_all))

    clusters <- rbind(c(10,0,0), c(10,10,0), c(10,10,10),
                      c(0,0,10), c(0,10,10), c(10,0,10))
    test <- rbind(as.matrix(x_meta[1,1:3]), clusters)

    common <- intersect(x_all, y_all)
    # Remove ranges that are only 10 bp long: 10bp chosen because of TF motifs being 10bp on average
    common <- common[width(common) >= 10]
    names(common) <- name_gr(common)

    x_scores <- mclapply(pooled_x, function(x) get_hits_scores(common, x, 'FE'))
    y_scores <- mclapply(pooled_y, function(x) get_hits_scores(common, x, 'FE'))

    mat <- matrix(, nrow= length(common), ncol=length(pooled_x)+length(pooled_y))
    rownames(mat) <- names(common)
    colnames(mat) <- c(names(pooled_x), names(pooled_y))

    # Fill the matrix with the scores for the hits
    foreach (x=1:ncol(mat)) %do%
    {
        if (colnames(mat)[x] %in% names(x_scores)) {
            mat[,x] <- as.integer(x_scores[[colnames(mat)[x]]])
        } else if (colnames(mat)[x] %in% names(y_scores)) {
            mat[,x] <- as.integer(y_scores[[colnames(mat)[x]]])
        }
        mat[is.na(mat[,x]),x] <- 0
    }

    # get number of clusters for kmeans 
    # based on http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters/15376462#15376462
    #wss <- (nrow(mat)-1) *sum(apply(mat, 2, var))
    #for (i in 2:50) wss[i] <- sum(kmeans(mat, centers=i)$withinss)
    #plot(1:length(wss), wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

    mat <- log2(mat)
    # remove -inf values (basically where there is no peak)
    mat[is.infinite(mat) == TRUE] <- -1
    mat <- as.big.matrix(mat)
    ans <- bigkmeans(mat, 5, nstart=30, iter.max=100)

    # cannot coerce bigmatrix to dataframe directly -- make matrix first
    pdata <- as.data.frame(as.matrix((mat)))
    pdata[['cluster']] <- factor(as.character(ans$cluster))
    pdata <- data.matrix(pdata[with(pdata, order(cluster)),])
    side_anno <- as.character(pdata[,ncol(pdata),drop=T])

    plotting()

}# }}}

pca <- function(matrix){# {{{
    x <- c('mclust', 'GGally')
    mclapply(x, suppressMessages(library), character.only=T)

    pca <- prcomp(matrix)
    mclust_clusters <- mclustBIC(matrix)
    clusters <- summary(mclust_clusters, data=matrix)$classification
    
    plot_data <- as.data.frame(pca$x)
    plot_data['clusters'] <- as.factor(unname(clusters[match(names(clusters),
                                                             rownames(plot_data))]))

    pdf("ggpairs_pca.pdf")
    g <- ggpairs(plot_data, color="clusters", axisLabels="internal")
    # necessary to print like this otherwise pdf() does not work with ggpairs
    print(g)
    dev.off()

    return(list(pca=pca, clusters=clusters))
}

# }}}
# check if .Rda object in output directory already
rda_object <- file.path(output_path, args$label)
if(file.exists(rda_object){
   load(rda_object)
   #No need to pass arguments as those are global and already loaded
   plotting() 
}else{
    main()
}

#reached memory limit
#mat_apclus <- apcluster(negDistMat(r=2), mat)
#save(mat_apclus, file='mat_apclus.Rdata')
#quit()

#oct4_gain <- mat[,'oct4-2i-2lox'] == 0 & mat[,'mbd3-2i'] == 1 & mat[,'mbd3-EpiSCs'] == 0  & mat[,'oct4-EpiLCs'] == 1# {{{
#oct4_loss <- mat[,'oct4-2i-2lox'] == 1 & mat[,'mbd3-2i'] == 0 & mat[,'mbd3-EpiSCs'] == 1  & mat[,'oct4-EpiLCs'] == 0
#
## check if the peaks overlap genes etc 
#gain <- genes[subjectHits(findOverlaps(oct4_all[oct4_gain], genes)), 'gene_id']
#loss <- elementMetadata(genes[unique(subjectHits(findOverlaps(oct4_all[oct4_loss], genes)))])[['gene_id']]
#gain <- elementMetadata(genes[unique(subjectHits(findOverlaps(oct4_all[oct4_gain], genes)))])[['gene_id']]
#
#oct4_over_gene <- as.integer(overlapsAny(oct4_all, genes))
#names(oct4_over_gene) <- names(oct4_all)
#
#dev_time <- unlist(mclapply(c('0hr', '2i', 'SL', '16hr', 'EpiLCs', 'EpiSCs'), function(x){colnames(mat)[grep(x, colnames(mat))]}))
#pdata <- mat[,dev_time]
#pdata <- pdata[do.call(order, as.data.frame(pdata)),]
#pdata <- cbind(pdata, 'overlaps_gene'=as.integer(oct4_over_gene[rownames(pdata)]))
#
#cool_cols <- colorRampPalette(c('aliceblue','darkcyan'))(100)
#
#pdf("plots/oct4matrix.pdf")
#pheatmap(pdata, #trace= "none", density.info= "none",
#         color=cool_cols,
#         clustering_distance_cols= "binary", clustering_method= "ward",
#         cluster_cols=FALSE, cluster_rows=FALSE,
#         border_color=NA,
#         fontsize=10, fontsize_row=7, show_rownames=FALSE, show_colnames=TRUE)
#dev.off()
#
#pdata_rev <- mat[,rev(dev_time)]
#pdata_rev <- pdata_rev[do.call(order, as.data.frame(pdata_rev)),]
#pdata_rev <- cbind(pdata_rev, 'overlaps_gene'=as.integer(oct4_over_gene[rownames(pdata_rev)]))
#
#
#pdf("plots/oct4matrix-rev.pdf")
#pheatmap(pdata_rev, #trace= "none", density.info= "none",
#         color=cool_cols,
#         clustering_distance_cols= "binary", clustering_method= "ward",
#         cluster_cols=FALSE, cluster_rows=FALSE,
#         border_color=NA,
#         fontsize=10, fontsize_row=7, show_rownames=FALSE, show_colnames=TRUE)
#dev.off()
#
#source("~/source/Rscripts/functions.R")
#library(topGO)
#go <- factor(as.integer(elementMetadata(genes)[['gene_id']] %in% gain))
#names(go) <- elementMetadata(genes)[['gene_id']]
#
#go <- GOtermEnr(go, 'oct4-gain', "BP", "org.Mm.eg.db")
#
#library(gridExtra)
#
#pdf(file.path('plots',"oct-gain_GO.pdf"), paper='a4')
#p <- ggplot(go$table,aes(y=Significant/Annotated,x=Term,fill=weight01Fisher)) + geom_bar(stat="identity") + scale_fill_gradient(low="red",high="yellow")   + coord_flip()
#p
#dev.off()
#
#markers <- read.table("/nfs/research2/bertone/user/mxenoph/hendrich/markers.txt", header=T, stringsAsFactors=F)
#markers <- lapply(as.list(markers), function(x){
#                  x <- x[!is.na(x)]
#                  id <- mclapply(x, function(n) values(genes[values(genes)[['gene_name']] == n ])[['gene_id']] )
#                  names(id) <- x
#                  unlist(id)
#         })# }}}
