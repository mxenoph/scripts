######functions.R#######
library = function (...) suppressMessages(base::library(...))

##########################################
get.lengths = function(gtf,name,feature) {# {{{
    #by Remco
    # List of packages to load
    x = c('GenomicRanges')
    mclapply(x, suppressMessages(library), character.only=T)

    colnames(gtf) = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr")
    gtf$g_id = gsub(".*gene_id ","",gtf$attr)
    gtf$g_id = gsub(";.*","",gtf$g_id)
    gtf$t_id = gsub(".*transcript_id ","",gtf$attr)
    gtf$t_id = gsub(";.*","",gtf$t_id)
    gtf = gtf[gtf$feature ==feature,]
    gene_ranges = GRanges(seqnames = gtf$seqname, ranges = IRanges(start = gtf$start, end = gtf$end, names = gtf$g_id),
                        strand = gtf$strand, source = gtf$source)
    gene_list = split(gene_ranges, as.factor(names(gene_ranges)))

    reduced_list = endoapply(gene_list, reduce)

    lengths = as.data.frame(sapply(reduced_list,function(x){sum(width(x))}))

    write.table(lengths, file = paste(name,"_genelengths.txt",sep = ""), col.names = FALSE, sep = "\t", quote = FALSE)
    return(list(matrix = lengths))
}# }}}

#######################################i#
compute.FPKMS=function(counts, length){# {{{
    #Give it a matrix for exons or introns and a matrix with the combined lengths of the respective features.
    total = NULL
    for(c in 1:ncol(counts)){
        total = cbind(total, sum(counts[,c]))
    }
    total = total/10^6

    fpkm = merge(counts, length, by= "row.names")
    rownames(fpkm) = fpkm[,1]
    remove = "Row.names"
    fpkm = fpkm[, !colnames(fpkm) %in% remove]
    myncol = ncol(fpkm)
    m = (ncol(fpkm))-1
    for(i in 1:m){
        fpkm[,i] = fpkm[,i] / ((fpkm[,myncol] / 10^3) * (total[1,i]))
    }
    return(list(matrix = fpkm, total = total))
}# }}}

draw.heatmap=function(matrix, list, name, mgi) {# {{{
#by Remco
    # List of packages to load
    x = c('gplots', 'RColorBrewer')
    mclapply(x, suppressMessages(library), character.only=T)

    list = list[list %in% rownames(matrix)]
    width = length(list) / 5

    if (length(list) > 1) {
    list.matrix = matrix[list,]
    plot.data = as.matrix(log(list.matrix+1))

    plot.data = plot.data-rowMeans(plot.data)

    # colours
    heatcol = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(255)
    # rownames to be set as MGI symbols
    rnames = as.vector(mgi[rownames(plot.data),])
    # For genes that an MGI symbol not available use ensembl id
    rnames[which(is.na(rnames))] = rownames(plot.data)[which(is.na(rnames))]

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

      #legend.names=colnames(plot.data)
      #legend("bottomleft",legend=legend.names,fill=rainbow(4), bty="n", cex=0.4)

      dev.off()
    }
    }
}
# }}}


#######################################
plot.pca = function(matrix, cols, filename, leg, pch) {# {{{
    # List of packages to load
    x = c('scatterplot3d')
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
    s3d = scatterplot3d(pca$rotation[,1:3], type= "h",color= cols, pch=pch, main= filename)
    s3d.coords = s3d$xyz.convert(pca$rotation[,1], pca$rotation[,2], pca$rotation[,3])
    text(s3d.coords$x, s3d.coords$y,     # x and y coordinates
       labels= colnames(matrix),       # text to plot
       pos= 3, cex= .6)

    dev.off()
}# }}}


compute.FPKMS.v2 = function(counts, length, scaleFactors){# {{{
    #Give it a matrix for exons or introns and a matrix with the combined lengths of the respective features.
    #Correct for the scale factors- t() transposes the matrix
    corrMatrix = t(t(counts) / scaleFactors)
    FPKMs = corrMatrix / colSums(corrMatrix) * 1000000
    FPKMs = round(FPKMs / length$V2 * 1000, digits= 4)
    return(FPKMs)
}# }}}

#The geneUniverse are factors indicating which gees are interesting in the Universe
#e.g.GOtermEnr(Univ, Sel, "org.Mm.eg.db", "BP")
GOtermEnr = function(...) get_set_enrichment(...)

get_set_enrichment = function(gene_scores, ontology = "BP", organism = "org.Mm.eg.db", label = date(), padj = 0.05, bound){# {{{
    # List of packages to load
    x = c('topGO', 'dplyr', 'ggplot2', organism)
    lapply(x, suppressMessages(library), character.only=T)
    gene_scores[is.na(gene_scores)] = 1
    diff_genes = function(gene_scores, padj = 0.05) {
        return(gene_scores < padj)
    }
    
    if(missing(bound)){# {{{
        # perform GSEA on de genes
        go_data = new("topGOdata",
                     description = label, ontology = ontology,
                     allGenes = gene_scores, geneSel = diff_genes,
                     # prune GO hierarchy from trms with less than 10 annotated genes
                     nodeSize = 30,
                     annot = annFUN.org, mapping = organism, ID = "Ensembl")
    } else {
        # For performing GSEA not only on genes passing that are also bound by one criterion
        diff_and_bound_genes = function(gene_scores, subset = bound, padj = 0.05) {
            return(gene_scores < padj & names(gene_scores) %in% subset)
        }
        go_data = new("topGOdata",
                     description = label, ontology = ontology,
                     allGenes = gene_scores, geneSelectionFun = diff_and_bound_genes,
                     nodeSize = 10,
                     annot = annFUN.org, mapping = organism, ID = "Ensembl")
    }# }}}

    if( length(sigGenes(go_data)) < 10 ) {
        logs = "There are no feasible significant genes"
        print(logs)
        return(logs = logs)
    }

    result_fisher = runTest(go_data, algorithm = "weight01", statistic = "fisher" )
    result_ks = runTest(go_data, algorithm = "weight01", statistic = "ks" )

    # Plotting distribution of p values for both tests# {{{
    tmp = as.data.frame(score(result_fisher))
    colnames(tmp) = 'pval'

    p = ggplot(tmp, aes(x = pval)) + geom_histogram(binwidth = 0.02)
    p = p + xlab('raw p-values') + theme_classic()
    
    print(p + ggtitle('Fisher\'s exact test'))

    tmp = as.data.frame(score(result_ks))
    colnames(tmp) = 'pval'
    
    print(p %+% tmp + ggtitle('Kolmogorov-Smirnov test'))
# }}}
#    uncorrected = GenTable(go_data,
#                       Fisher = result_fisher,
#                       Kolmogorov_Smirnov = result_ks,
#                       orderBy = "Kolmogorov_Smirnov", topNodes = 20)
#    head(uncorrected)
#    uncorrected = uncorrected %>% mutate(Ratio = Significant / Annotated)
#    uncorrected = uncorrected %>% mutate_each(funs(as.numeric), 3:9)
#
#    p = ggplot(uncorrected, aes(y = Ratio, x = Term, fill = Kolmogorov_Smirnov))
#    p = p + geom_bar(stat = "identity") + scale_fill_gradient(low = "red", high = "yellow")
#    p = p + ylab('Significant / Annotated') + coord_flip() + theme_classic()
#    print(p + geom_text(aes(y = Ratio,
#                            x = Term,
#                            label = uncorrected[[6]], hjust = 1), size = 5) + ggtitle(paste0(label, ", uncorrected")))
#    
    # Adjusting for multiple testing
    score(result_fisher) = p.adjust(score(result_fisher), method="BH")
    # Adjusting the scores from ks test is not recommended as the multiple testing theory does not
    # directly apply

    # Plotting distribution of p values for both tests, after fdr correction # {{{
    tmp = as.data.frame(score(result_fisher))
    colnames(tmp) = 'pval'
    p = ggplot(tmp, aes(x = pval)) + geom_histogram(binwidth = 0.02)
    p = p + xlab('adjusted p-values') + theme_classic()
    print(p + ggtitle('Fisher\'s exact test'))

# }}}

    #Summarizing and ordering the results; add ranksOf="whatever" if you want to get them ranked by a specific test
    results = GenTable(go_data,
                       Fisher = result_fisher,
                       Kolmogorov_Smirnov = result_ks,
                       orderBy = "Kolmogorov_Smirnov", topNodes = 20)

    results = results %>% mutate(Ratio = Significant / Annotated)
    results = results %>% mutate_each(funs(as.numeric), 3:9)
    significant = results %>% filter(Kolmogorov_Smirnov < 0.05)

    # Makes no sense to get the GO graph if there are no significant terms
    if((significant %>% nrow()) > 0) {
        showSigOfNodes(go_data, score(result_fisher), firstSigNodes = 10, useInfo = 'all')
        go_id = results %>% dplyr::select(GO.ID) %>% .[[1]]
        # Error with chip = organism, see https://support.bioconductor.org/p/68661/
        # still not solved so do not run for the moment
        #gene_table = printGenes(go_data, whichTerms = go_id[1], chip = organism, numChar = 40)
    }

    p = ggplot(results, aes(y = Ratio, x = Term, fill = Kolmogorov_Smirnov)) + geom_bar(stat = "identity", aes(order=desc(Kolmogorov_Smirnov)))
    p = p + scale_fill_gradient(low="red",high="yellow") + coord_flip() + ggtitle(label) + theme_classic()
    print(p)

    return(list(go_data = go_data, table = results))

}# }}}

###Writing granges object to gtf
###To DO need to make it take the names of element metadata and print that
granges.to.gtf = function (gr, gtf, feature) {# {{{
    x = c('GenomicRanges')
    mclapply(x, suppressMessages(library), character.only=T)

    write.data = sprintf  ("%s\tgranges.to.gtf\t%s\t%d\t%d\t.", seqnames (gr), feature, start (gr), end (gr))
    strands = as.character (strand (gr))
    strands[strands=='*'] = '.'
    write.data = sprintf ("%s\t%s\t.", write.data, strands)
    if(feature %in% c("exon", "intron", "gene")){
     groups = sprintf ("gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\";",
                     elementMetadata (gr)$Gene, elementMetadata (gr)$Symbol, elementMetadata (gr)$Gene)
    }
    else{
     groups = sprintf ("%s_id \"%s\"; spanning.exons \"%s\"; spanning.introns \"%s\"; spanning.intergenic \"%s\"; spanning.intergenic.both.strands \"%s\";",
                     feature, names(gr), elementMetadata (gr)$spanning.exons, elementMetadata (gr)$spanning.introns, elementMetadata (gr)$spanning.intergenic, elementMetadata (gr)$spanning.intergenic.both.strands)
     #tmp = paste(unlist(lapply(1:length(names(elementMetadata(gr))), function(i){tmp = paste(names(elementMetadata(gr))[i], "%s;", sep=' ' )})), collapse=' ')
     #tmp =  paste("%s_id %s;", tmp, sep=' ')
     #groups = sprintf (tmp, feature, elementMetadata(gr)[1], elementMetadata(gr)[2], elementMetadata(gr)[3] )
    }
    write.data = sprintf ("%s\t%s", write.data, groups)

    writeLines (write.data,gtf)
}# }}}

drawVenn = function(sets){# {{{
    x = c('Vennerable', 'RColorBrewer')
    mclapply(x, suppressMessages(library), character.only=T)
    
    vn = Venn(sets)
    c_vn = compute.Venn(vn)
    theme_vn = VennThemes(c_vn)
    cols = brewer.pal(9, "Set1")
    
    for (i in 1:length(theme_vn[["Set"]])) {
        theme_vn[["Set"]][[i]]$lwd = length(theme_vn[["Set"]]) + 1 - i
        theme_vn[["Set"]][[i]]$col = cols[i]

        theme_vn[["SetText"]][[i]]$lwd = length(theme_vn[["SetText"]]) + 1 - i
        theme_vn[["SetText"]][[i]]$col = cols[i]
    }

    for (i in 1:length(theme_vn[["FaceText"]])) {
        theme_vn[["Face"]][[i]]$col = "white"
        theme_vn[["Face"]][[i]]$fill = "white"

        # Identifier for face e.g. 00001 is Set1 for venn of 5 sets
        # If the current text does not corresponds to an intersection, give it the same colour as the Set
        if( length(grep(1, unlist(strsplit(names(theme_vn[["FaceText"]])[i], '')))) == 1){
            set = paste0('Set', grep(1, unlist(strsplit(names(theme_vn[["FaceText"]])[i], ''))))
            theme_vn[["FaceText"]][[i]]$col = theme_vn[["SetText"]][[set]]$col
        }
        else{
            theme_vn[["FaceText"]][[i]]$col = colors()[295]
            # TODO: email the developer about how to control FaceText colour and position; especially the 100-1. Give examples
#            print(names(theme_vn[["FaceText"]])[i])
        }
    }

    vn_plot = plot(vn, type="ChowRuskey", doWeights=TRUE, show=list(SetLabels=TRUE, FaceLabels=FALSE), gp=theme_vn)
}# }}}

# Sets should be a named list containing the 
draw_venn_diagram = function(sets){# {{{
    x = c('Vennerable', 'RColorBrewer')
    mclapply(x, suppressMessages(library), character.only=T)
    
    cols = brewer.pal(9, "Set1")

    c1 = Venn(SetNames = c('NuRD', 'MTA1 FLAG', 'MTA2 GFP'), Weight = c(`000`=0, `100` = 4563, `110`=1071, `111` = 5597, `101`=3136, `010`=167, `011`=45, `001`=237))
    c2 = Venn(SetNames = c('NuRD', 'MTA1', 'MTA2'), Weight = c(`000`=0, `100` = 32544, `110`=10524, `111` = 19820, `101`=6606, `010`=5947, `011`=1079, `001`=1541))
    c3 = Venn(SetNames = c('chd4', 'MTA1 FLAG', 'MTA2 GFP'), Weight = c(`000`=0, `100` = 115218, `110`=13457, `111` = 20573, `101`=7096, `010`=2472, `011`=360, `001`=514))
    c4 = Venn(SetNames = c('Mbd3', 'MTA1 FLAG', 'MTA2 GFP'), Weight = c(`000`=0, `100` =37722, `110`=10805, `111` = 19968, `101`=6760, `010`=5504, `011`=937, `001`=1396))
    plot(c1,  doWeights = TRUE, type = "circles")
}
# }}}

# TODO: move to functions.R and source that file
deseq2vect = function(f, padj=0.05){# {{{
    raw = df = read.table(f, header=TRUE, sep="\t", row.names='id')

    df = df[!is.na(df$padj),]
    #df[df$log2FoldChange== -Inf, 'log2FoldChange'] = min(df[is.finite(df$log2FoldChange), 'log2FoldChange'])
    df[df$log2FoldChange== -Inf, 'log2FoldChange'] = -5
    #df[df$log2FoldChange== Inf, 'log2FoldChange'] = max(df[is.finite(df$log2FoldChange), 'log2FoldChange'])
    df[df$log2FoldChange== Inf, 'log2FoldChange'] = 5

    up = down = de = df[df$padj < padj, ]
    up = up[up$log2FoldChange > 0, ]; down = down[down$log2FoldChange < 0, ]
    v = setNames(df$padj, rownames(df)); v[is.na(v)] = 1; v = v < padj

    return(list('raw'=raw, 'df'= df, 'de'=de, 'up'=up, 'down'=down, 'binary'=v))
}
# }}}

#ds is the dataset e.g. "mmusculus_gene_ensembl" ----NEEDS EDITING
ensembl2entrezIDs= function(ds){
    x = c('biomaRt')
    mclapply(x, suppressMessages(library), character.only=T)
    ensembl=useMart("ensembl", dataset=ds)
}
