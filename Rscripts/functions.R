######functions.R#######

##########################################
get.lengths <- function(gtf,name,feature) {# {{{
    #by Remco
    # List of packages to load
    x <- c('GenomicRanges')
    mclapply(x, suppressMessages(library), character.only=T)

    colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr")
    gtf$g_id = gsub(".*gene_id ","",gtf$attr)
    gtf$g_id = gsub(";.*","",gtf$g_id)
    gtf$t_id = gsub(".*transcript_id ","",gtf$attr)
    gtf$t_id = gsub(";.*","",gtf$t_id)
    gtf = gtf[gtf$feature ==feature,]
    gene_ranges <- GRanges(seqnames = gtf$seqname, ranges = IRanges(start = gtf$start, end = gtf$end, names = gtf$g_id),
                        strand = gtf$strand, source = gtf$source)
    gene_list <- split(gene_ranges, as.factor(names(gene_ranges)))

    reduced_list = endoapply(gene_list, reduce)

    lengths <- as.data.frame(sapply(reduced_list,function(x){sum(width(x))}))

    write.table(lengths, file = paste(name,"_genelengths.txt",sep = ""), col.names = FALSE, sep = "\t", quote = FALSE)
    return(list(matrix = lengths))
}# }}}

#######################################i#
compute.FPKMS<-function(counts, length){# {{{
    #Give it a matrix for exons or introns and a matrix with the combined lengths of the respective features.
    total <- NULL
    for(c in 1:ncol(counts)){
        total <- cbind(total, sum(counts[,c]))
    }
    total <- total/10^6

    fpkm <- merge(counts, length, by= "row.names")
    rownames(fpkm) <- fpkm[,1]
    remove <- "Row.names"
    fpkm <- fpkm[, !colnames(fpkm) %in% remove]
    myncol <- ncol(fpkm)
    m <- (ncol(fpkm))-1
    for(i in 1:m){
        fpkm[,i] <- fpkm[,i] / ((fpkm[,myncol] / 10^3) * (total[1,i]))
    }
    return(list(matrix = fpkm, total = total))
}# }}}

draw.heatmap<-function(matrix, list, name, mgi) {# {{{
#by Remco
    # List of packages to load
    x <- c('gplots', 'RColorBrewer')
    mclapply(x, suppressMessages(library), character.only=T)

    list <- list[list %in% rownames(matrix)]
    width <- length(list) / 5

    if (length(list) > 1) {
    list.matrix <- matrix[list,]
    plot.data <- as.matrix(log(list.matrix+1))

    plot.data <- plot.data-rowMeans(plot.data)

    # colours
    heatcol <- colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(255)
    # rownames to be set as MGI symbols
    rnames <- as.vector(mgi[rownames(plot.data),])
    # For genes that an MGI symbol not available use ensembl id
    rnames[which(is.na(rnames))] <- rownames(plot.data)[which(is.na(rnames))]

    #plot
    if (nrow(plot.data) > 1) {
      pdf(paste("Heatmap_", name,".pdf", sep= ""), height= max(width,5), width= 7)
      par(xpd = NA, mar= c(5,5,5,16))
      par(oma= c(1.5,2,1,1))

      heatmap.2(plot.data, trace= "none", density.info= "none",
                symm= FALSE,cexRow= .5, cexCol= 1, labRow=rnames, Colv= FALSE, Rowv= FALSE,
                main= name, keysize=1.5, col= heatcol, dendrogram= "none")

      heatmap.2(plot.data, trace= "none", density.info= "none",
                symm= FALSE, cexRow= .5, cexCol= 1, labRow= rnames,
                Colv= FALSE, Rowv= TRUE, keysize= 1.5, main= name, col= heatcol, dendrogram= "row")

      #legend.names<-colnames(plot.data)
      #legend("bottomleft",legend=legend.names,fill=rainbow(4), bty="n", cex=0.4)

      dev.off()
    }
    }
}
# }}}


#######################################
plot.pca <- function(matrix, cols, filename, leg, pch) {# {{{
    # List of packages to load
    x <- c('scatterplot3d')
    mclapply(x, suppressMessages(library), character.only=T)
    #by Remco

    pdf(paste("PCA_", filename,".pdf",sep= ""), width= 12, height= 8)
    par(xpd = NA, mar= c(5,5,5,16))
    pca= princomp(matrix)
    plot(pca$loadings, main= filename, pch=pch)
    points(pca$loadings, pch= pch, cex=1, lwd=2.5, col= cols)
    #text(pca$loadings,colnames(matrix),pos=3,cex=0.8)
    legend("topright", inset= c(-0.4,0), legend=leg, bty= "n", pch=pch, cex=1, col=cols)
    dev.off()

    pdf(paste("3D_PCA_", filename,".pdf",sep= ""), width= 12, height= 8)
    par(xpd = NA, mar= c(5,5,5,16))
    pca = prcomp(matrix)
    s3d <- scatterplot3d(pca$rotation[,1:3], type= "h",color= cols, pch=pch, main= filename)
    s3d.coords <- s3d$xyz.convert(pca$rotation[,1], pca$rotation[,2], pca$rotation[,3])
    text(s3d.coords$x, s3d.coords$y,     # x and y coordinates
       labels= colnames(matrix),       # text to plot
       pos= 3, cex= .6)

    dev.off()
}# }}}


compute.FPKMS.v2 <- function(counts, length, scaleFactors){# {{{
    #Give it a matrix for exons or introns and a matrix with the combined lengths of the respective features.
    #Correct for the scale factors- t() transposes the matrix
    corrMatrix <- t(t(counts) / scaleFactors)
    FPKMs <- corrMatrix / colSums(corrMatrix) * 1000000
    FPKMs <- round(FPKMs / length$V2 * 1000, digits= 4)
    return(FPKMs)
}# }}}

#The geneUniverse are factors indicating which gees are interesting in the Universe
#e.g.GOtermEnr(Univ, Sel, "org.Mm.eg.db", "BP")
GOtermEnr <- function(geneUniverse, info, ontology, orgDB){# {{{
    # List of packages to load
    x <- c('topGO', 'ggplot2')
    mclapply(x, suppressMessages(library), character.only=T)

    GOdata <- new("topGOdata",
                    description = info, ontology = ontology,
                    allGenes = geneUniverse,
                    nodeSize = 10,
                    annot = annFUN.org, mapping = orgDB, ID = "Ensembl")


    #SOS
    #Ask Remco if fisher is the appropriate test to use and what are the gene-wise scores one uses for the KS test
    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher" )
    #Adjusting for multiple testing
    score(resultFisher) <- p.adjust(score(resultFisher), method="BH")
    #resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks" )

    #Make a histogram of the p values returned by the test statistic
    # hist(score(resultFisher), 50, xlab = "p-values")

    #Summarizing and ordering the results; add ranksOf="whatever" if you want to get them ranked by a specific test
    #allres <- GenTable(GOdata, weight01Fisher = resultFisher, weight01KS = resultKS, orderBy = "weight01Fisher", topNodes = 20)
    allRes <- GenTable(GOdata, weight01Fisher = resultFisher, orderBy = "weight01Fisher", topNodes = 20)
    allRes$weight01Fisher <- as.numeric(allRes$weight01Fisher)
    #allRes$Term <- paste(allRes$Term,allRes$GO.ID)
    allRes$Term <- factor(allRes$Term, levels=allRes$Term[order(allRes$Significant/allRes$Annotated,decreasing=F)])
    #fn.prefix should be the name of the pair you are comparing
    # printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

    pdf(paste0(info,"GO.pdf"), paper='a4')
    ggplot(allRes,aes(y=Significant/Annotated,x=Term,fill=weight01Fisher)) + geom_bar(stat="identity", aes(order=desc(weight01Fisher))) + scale_fill_gradient(low="red",high="yellow")   + coord_flip()
    dev.off()

    return(list(GOdata = GOdata, table = allRes, resFis = resultFisher))

}# }}}

###Writting granges object to gtf
###To DO need to make it take the names of element metadata and print that
granges.to.gtf <- function (gr, gtf, feature) {# {{{
    x <- c('GenomicRanges')
    mclapply(x, suppressMessages(library), character.only=T)

    write.data <- sprintf  ("%s\tgranges.to.gtf\t%s\t%d\t%d\t.", seqnames (gr), feature, start (gr), end (gr))
    strands <- as.character (strand (gr))
    strands[strands=='*'] <- '.'
    write.data <- sprintf ("%s\t%s\t.", write.data, strands)
    if(feature %in% c("exon", "intron", "gene")){
     groups <- sprintf ("gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\";",
                     elementMetadata (gr)$Gene, elementMetadata (gr)$Symbol, elementMetadata (gr)$Gene)
    }
    else{
     groups <- sprintf ("%s_id \"%s\"; spanning.exons \"%s\"; spanning.introns \"%s\"; spanning.intergenic \"%s\"; spanning.intergenic.both.strands \"%s\";",
                     feature, names(gr), elementMetadata (gr)$spanning.exons, elementMetadata (gr)$spanning.introns, elementMetadata (gr)$spanning.intergenic, elementMetadata (gr)$spanning.intergenic.both.strands)
     #tmp <- paste(unlist(lapply(1:length(names(elementMetadata(gr))), function(i){tmp <- paste(names(elementMetadata(gr))[i], "%s;", sep=' ' )})), collapse=' ')
     #tmp <-  paste("%s_id %s;", tmp, sep=' ')
     #groups <- sprintf (tmp, feature, elementMetadata(gr)[1], elementMetadata(gr)[2], elementMetadata(gr)[3] )
    }
    write.data <- sprintf ("%s\t%s", write.data, groups)

    writeLines (write.data,gtf)
}# }}}

drawVenn <- function(sets){# {{{
    x <- c('Vennerable', 'RColorBrewer')
    mclapply(x, suppressMessages(library), character.only=T)
    
    vn <- Venn(sets)
    c_vn <- compute.Venn(vn)
    theme_vn <- VennThemes(c_vn)
    cols <- brewer.pal(9, "Set1")
    
    for (i in 1:length(theme_vn[["Set"]])) {
        theme_vn[["Set"]][[i]]$lwd <- length(theme_vn[["Set"]]) + 1 - i
        theme_vn[["Set"]][[i]]$col <- cols[i]

        theme_vn[["SetText"]][[i]]$lwd <- length(theme_vn[["SetText"]]) + 1 - i
        theme_vn[["SetText"]][[i]]$col <- cols[i]
    }

    for (i in 1:length(theme_vn[["FaceText"]])) {
        theme_vn[["Face"]][[i]]$col <- "white"
        theme_vn[["Face"]][[i]]$fill <- "white"

        # Identifier for face e.g. 00001 is Set1 for venn of 5 sets
        # If the current text does not corresponds to an intersection, give it the same colour as the Set
        if( length(grep(1, unlist(strsplit(names(theme_vn[["FaceText"]])[i], '')))) == 1){
            set <- paste0('Set', grep(1, unlist(strsplit(names(theme_vn[["FaceText"]])[i], ''))))
            theme_vn[["FaceText"]][[i]]$col <- theme_vn[["SetText"]][[set]]$col
        }
        else{
            theme_vn[["FaceText"]][[i]]$col <- colors()[295]
            # TODO: email the developer about how to control FaceText colour and position; especially the 100-1. Give examples
#            print(names(theme_vn[["FaceText"]])[i])
        }
    }

    vn_plot <- plot(vn, type="ChowRuskey", doWeights=TRUE, show=list(SetLabels=TRUE, FaceLabels=FALSE), gp=theme_vn)
}# }}}

# TODO: move to functions.R and source that file
deseq2vect <- function(f, padj=0.05){# {{{
    raw <- df <- read.table(f, header=TRUE, sep="\t", row.names='id')

    df <- df[!is.na(df$padj),]
    #df[df$log2FoldChange== -Inf, 'log2FoldChange'] <- min(df[is.finite(df$log2FoldChange), 'log2FoldChange'])
    df[df$log2FoldChange== -Inf, 'log2FoldChange'] <- -5
    #df[df$log2FoldChange== Inf, 'log2FoldChange'] <- max(df[is.finite(df$log2FoldChange), 'log2FoldChange'])
    df[df$log2FoldChange== Inf, 'log2FoldChange'] <- 5

    up <- down <- de <- df[df$padj < padj, ]
    up <- up[up$log2FoldChange > 0, ]; down <- down[down$log2FoldChange < 0, ]
    v <- setNames(df$padj, rownames(df)); v[is.na(v)] <- 1; v <- v < padj

    return(list('raw'=raw, 'df'= df, 'de'=de, 'up'=up, 'down'=down, 'binary'=v))
}
# }}}

#ds is the dataset e.g. "mmusculus_gene_ensembl" ----NEEDS EDITING
ensembl2entrezIDs<- function(ds){
    x <- c('biomaRt')
    mclapply(x, suppressMessages(library), character.only=T)
    ensembl<-useMart("ensembl", dataset=ds)
}
