library("DESeq2")
library("plyr")
library("ggplot2")
library("gplots")
library("gridBase")
library("gridExtra")
library("RColorBrewer")
library("reshape")

# green, hotpink, lime, purple, light pinkish, dark red, mustard, brown, grey, blue
default_colors = c("#76C2A0","#A94777",
                   "#87CE4F", "#B25FCD",
                   "#C6AEA9", "#C05038",
                   "#B99F47", "#494235",
                   "#A2AFB9", "#6F7FB5")
names(default_colors) = c('green', 'hotpink',
                           'lime', 'purple',
                           'lightpinkish', 'darkred',
                           'mustard', 'brown',
                           'grey','blue')

# More elegant box plot.- should be moved to plotting functions and source that file
boxplot = function (..., col) {
    # Plot very narrow, minimal boxes.
    pars = list(boxwex = 0.4, staplewex = 0,
                medlwd = 1, whisklty = 1,
                frame.plot = F, outpch = 20, outcex = 0.3)
    if (! missing(col))
        pars = c(pars, list(boxcol = col, whiskcol = col, outcol = col))
    graphics::boxplot(..., pars = pars)
}

# Color scheme for divergent colors.
divergent_colors = colorRampPalette(c('#603D71', 'white', '#A4B962'))(50)
divergent_warm_colors = colorRampPalette(rev(brewer.pal(9,"RdBu")))(50)

# Less intrusive defaults for heatmap.2.
heatmap.2 = function (...)
        gplots::heatmap.2(..., trace = 'none', density.info = 'none',
                          col = divergent_colors, margins = c(12,12))

### FUNCTIONS ###
# create a master file with all counts from the counts directory# {{{
export_counts = function(counts_path) {
    if (file.exists(file.path(counts_path, 'experiment-counts.tsv'))) {
        #check.names=FALSE ensures that no X is put in front of colnames starting with a number
        counts = read.table(file.path(counts_path, 'experiment-counts.tsv'), row.names = 1, header = T, sep = "\t", check.names = FALSE)
    } else {
        counts_files = list.files(counts_path, pattern = "\\.counts$", full.name = T)
        counts_per_sample <- mclapply(counts_files, read.delim, row.names = 1, header = F)
        names(counts_per_sample) = basename(file_path_sans_ext(counts_files))
        
        # Ensure that all files have the same identifiers before merging
        check = unlist(lapply(2:length(counts_per_sample)-1, function(i){
                              s=i+1
                              unlist(lapply(s:length(counts_per_sample), function(j) {
                                            check = identical(rownames(counts_per_sample[[i]]), rownames(counts_per_sample[[j]]))
                                            names(check) = paste0(names(counts_per_sample[i]), ':', names(counts_per_sample[j]))
                                            return(check) })) }))
        if (any(check == FALSE)) {print(check); stop('Count files do not contain the same order and number feature IDs. Please inspect the files and rerun.')}
        
        counts = do.call(cbind.data.frame, counts_per_sample)
        names(counts) = basename(file_path_sans_ext(counts_files))
        write.table(counts, file.path(counts_path, 'experiment-counts.tsv'), row.names = T, col.names = T, quote = F, sep = "\t")
    }
    # tested against DESeqDataSetFromHTSeqCount the size factors are computed only on the genes, not the 
    # ambiguous and no features
    counts = add_rownames(counts, var = 'Gene')
    # Cleaning up htseq-count output
    counts = counts %>% filter(! Gene %in% c('__no_feature',
                                            '__ambiguous', 
                                            '__too_low_aQual',
                                            '__not_aligned',
                                            '__alignment_not_unique')) %>%
    droplevels()
    rownames(counts) = counts[['Gene']]
    counts = counts[,-1]
    return(counts)
}# }}}

# Used in plot_counts() to resetmar and mfrow# {{{
reset_par = function(){
    # Reset par to defaults in the parent environment
    eval(quote(par(mfrow = c(1,1))), parent.frame())
    eval(quote(par(oma = c(0,0,0,0))), parent.frame())
    eval(quote(par(mar = c(5, 4, 4, 2) + 0.1)), parent.frame())
}# }}}

# plot counts # {{{
plot_counts = function(eset, counts, design,
                        # colors should be a named vector
                       contrast = "Contrast", colors,
                       scaled = T, type = NULL, rld){

    library(gtools)
    # If missing create eset from counts# {{{
    if(missing(eset)){
        library(EDASeq)
        print('Building eset object from counts...')
        contrasts = design[match(colnames(counts), as.character(design$Library)), contrast]
        unique_contrast = unique(contrasts)
        contrasts = factor(contrasts, levels = unique_contrast)
        eset = newSeqExpressionSet(as.matrix(counts),
                                   phenoData = data.frame(contrasts, row.names = colnames(counts)))
    } else {
        tmp = DESeq2::counts(eset)
        tmp = colnames(tmp)
        # contrasts is used in plotRLE for the colors
        contrasts = design[match(tmp, as.character(design$Library)), contrast]
        unique_contrast = unique(contrasts)
        contrasts = factor(contrasts, levels = unique_contrast)
    } # }}}

    # Set colours if missing# {{{
    # If color not passed in the call will use the default colors
    if(missing(colors)){
        colors = default_colors
        colors = colors[1:length(unique_contrast)]
        names(colors) = unique_contrast
    }# }}}

    format_RLE = function(){# {{{
        # Plot very narrow, minimal boxes.
        pars = list(boxwex = 0.4, staplewex = 0,
                    medlwd = 1, whisklty = 1,
                    bty = "n", outpch = 20, outcex = 0.3,
                    # handle colors
                    boxcol = colors[contrasts], whiskcol = colors[contrasts], outcol = colors[contrasts])
        plotRLE(eset, outline = FALSE, ylim =c(-0.5, 0.5), col = colors[contrasts],
                las = 2, cex = 0.7, pars = pars, frame.plot = F)
        title(main = type, ylab = "Relative Log Expression")
        legend('topright', bty = "n", legend = names(colors[unique_contrast]), fill = colors[unique_contrast])
    }# }}}
   
    format_cntbox = function(){# {{{
        # Use global minimum value > 0 as pseudocount.
        eps = min(counts[counts > 0])
        boxplot(counts + eps, las = 2, frame.plot = F, log = 'y', col = colors[contrasts])
        title(main = type, ylab = "Counts")
    }# }}}

    format_rld_pca = function(){# {{{
        library(ggplot2)

        theme = theme_set(theme_bw())
        theme = theme_update(legend.position = "bottom",
                             panel.border = element_blank(),
                             axis.line.x = element_line(),
                             axis.line.y = element_line(),
                             panel.grid.major.x = element_blank())

        keep_groups = as.data.frame(colData(rld)) %>% select(-c(Library, File, sizeFactor)) %>% colnames()
        (data = DESeq2::plotPCA(rld, intgroup = keep_groups, returnData=TRUE))
        variance = round(100 * attr(data, "percentVar"))
        p = ggplot(data, aes(PC1, PC2, color = mixedsort(get(contrast)))) + geom_point(size = 3)
        p = p + scale_color_manual(values = colors)
        p = p + labs(x = paste0("PC1: ", variance[1],"% variance"), 
                     y = paste0("PC2: ", variance[2],"% variance"),
                     title = "PCA on rlog(counts)",
                     color = contrast)
        print(p)
    }# }}}

    # Calculates disctance from rld object by default otherwise calculates and plots correlation# {{{# {{{
    format_dheatmap = function(mat, corr = FALSE) {
        library(pheatmap)
        if(corr){
            distance_matrix = cor(counts, method = 'spearman')
            distances = as.dist(distance_matrix)
            label = "(Spearman correlation)"

            #design is inherited by parent function
            annotation = design %>% select(-c(Library, File))
            rownames(annotation) = colnames(design[['Library']])
        } else{
            distances = dist(t(assay(mat)))
            distance_matrix = as.matrix(distances)
            label = "(rlog - distance)"
            annotation = as.data.frame(colData(mat)) %>% select(-c(Library, File, sizeFactor))
            divergent_colors = rev(divergent_colors)
            
#            library(PoiClaClu)
#            poison_distances = PoissonDistance(t(counts))
#            distance_matrix = as.matrix(poison_distances$dd)

        }# }}}
        
        if(any(!sapply(annotation, is.numeric))){
            annotation = annotation %>% mutate_each_(funs(factor),
                                                    names(which(!sapply(annotation, is.numeric))))
        }

        anno_colors = list(colors)
        names(anno_colors) = contrast
        annotation = annotation[,contrast, drop=F]
        # Reverting colours so that smallest distance matches the 1 in correlation heatmap
        pheatmap(distance_matrix,
                 clustering_distance_rows = distances,
                 clustering_distance_cols = distances,
                 col = divergent_colors,
                 border_col = NA,
                 annotation_col = annotation,
                 annotation_colors = anno_colors,
                 main = paste("Sample-to-sample distance", label))
        reset_par()
    }# }}}

    layout(matrix(c(1, 2), nrow = 1), widths = c(0.5, 0.5))
    format_cntbox()
    format_RLE()

    # Creating all library combinations
    if(!missing(rld)){
        print('Rlog transformed data provided...')
        format_rld_pca()
        reset_par()
        format_dheatmap(rld)
    } else {
        format_dheatmap(counts, corr = TRUE)
        
    }
    reset_par()
    return()
} 
# }}}

# plot_fpkm()# {{{
# results is the full results table from DESeq, markers should have a column called ensembl_gene_id
# contrast should be a named vector where values are contrast e.g. wt/ko and names are library names
# markers should have a column called state
plot_fpkm = function(fpkms, markers, first_contrast = NULL, second_contrast = NULL, results = NULL, label = "TPM"){
    library(reshape2)
    source("~/source/Rscripts/ggplot-functions.R")

    get_basic_p = function(f_subset, f_contrast=NULL, s_contrast = NULL, default_colors, results = NULL){# {{{

        # format f_subset# {{{
        if(! is.null(f_contrast) & ! is.null(s_contrast)) {
            print("Two contrasts provided. The first one will be used to colour groups and the second one to make in facets.")
            f_subset = f_subset %>% mutate(First_contrast = variable) %>%
                mutate(First_contrast = ifelse(First_contrast %in% names(f_contrast),
                                               f_contrast[First_contrast], First_contrast)) %>%
                mutate(Second_contrast = variable) %>%
                mutate(Second_contrast = ifelse(Second_contrast %in% names(s_contrast), 
                                                s_contrast[Second_contrast], Second_contrast))
        } else if(! is.null(f_contrast)) {
            f_subset = f_subset %>% mutate(First_contrast = variable) %>%
                mutate(First_contrast = ifelse(First_contrast %in% names(f_contrast),
                                               f_contrast[First_contrast], First_contrast))
        } # }}}
        
        # If no contrast provided only plotting based on state which is set as grid
        gg = "#A2AFB9"
        if(any(!is.null(f_contrast), !is.null(s_contrast))){

            if(!is.null(results)){
                p = ggplot(f_subset, aes(x = external_gene_id, y = value, color = First_contrast, shape = FDR))
                p = p + scale_shape_manual(values = setNames(c(19, 1), c('< 0.05', '>= 0.05')))
            } else{
                p = ggplot(f_subset, aes(x = external_gene_id, y = value, color = First_contrast))
            }
            # http://stackoverflow.com/questions/34734218/geom-errorbar-no-stat-called-stathline
            p = p + stat_summary(fun.y = mean, aes(ymin =..y.., ymax =..y.., color = First_contrast),
                                 geom = 'errorbar', width = 0.3, size = 1)
            # if first contrast is wt/ko then set the colours to default colors used in plot counts
            if(all(sort(levels(as.factor(f_contrast))) == c("ko", "wt"))) gg = setNames(default_colors[1:2], c("wt", "ko"))
            
            if(!is.null(second_contrast)) {
                p = p + facet_grid(Second_contrast ~ state)
            } else {
                p = p + facet_grid(~state)
            }
        } else {
            if(!is.null(results)){
                p = ggplot(f_subset, aes(x = external_gene_id, y = value, shape = FDR))
                p = p + scale_shape_manual(values = setNames(c(19, 1), c('< 0.05', '>= 0.05')))
            } else{
                p = ggplot(f_subset, aes(x = external_gene_id, y = value))
            }
            p = p + facet_grid(~state)
        }
        p = p + geom_point(size = 2) + scale_color_manual(values = gg)
        
        txt_angle = 20
        if(f_subset %>% dplyr::select(Gene) %>% unique() %>% .[[1]] %>% length() > 12) txt_angle = 90

        p = p + theme_bw() + theme(strip.text.x = element_text(size = 11, colour = "black"),
                                   axis.text.x = element_text(angle = txt_angle, size = 11)) 
        p = p + xlab('') + ylab(label)
        return(p)

    }# }}}

    markers = markers %>% dplyr::rename(Gene = ensembl_gene_id)
    fpkms_subset = fpkms %>% inner_join(markers, by = "Gene")
    if(!is.null(results)){
        fpkms_subset = fpkms_subset %>% inner_join(as.data.frame(results) %>% dplyr::select(Gene, padj), by = "Gene")
        fpkms_subset = fpkms_subset %>% mutate(FDR = ifelse(padj < 0.05, '< 0.05', '>= 0.05')) %>% dplyr::select(-padj)
        fpkms_subset = reshape2::melt(fpkms_subset, id = c("Gene", "external_gene_id", "state", "FDR"))
    } else {
        fpkms_subset = reshape2::melt(fpkms_subset, id = c("Gene", "external_gene_id", "state"))
    }
    # Not doing log transformation with ggplot because I do not get the axis to look like I want it to despite:
    # http://stackoverflow.com/questions/11214012/set-only-lower-bound-of-a-limit-for-ggplot
    # p = p + scale_y_continuous(trans = "log10",  name = label, limits = c(1, NA), breaks = scales::pretty_breaks())
    fpkms_subset = fpkms_subset %>% mutate(value = log10(value))
    label = paste0("log10(", label, ")")

    # If groups are too big then split the dataframe and plot separately
    if(nrow(markers) > 10){
        fpkms_list = split(fpkms_subset, fpkms_subset$state)
        lapply(fpkms_list, function(f_subset){
                   print(get_basic_p(f_subset = f_subset, f_contrast = first_contrast,
                                     s_contrast = second_contrast, default_colors = default_colors, results = results))
                                   })
    } else {
        lapply(fpkms_list, function(f_subset){
                   print(get_basic_p(f_subset = f_subset, f_contrast = first_contrast,
                                  s_contrast = second_contrast, default_colors = default_colors, results = results))})
    }

}# }}}

# Plotting density# {{{
plot_density = function(deseq_res, subsets = NULL, subset_name = NULL, label = "NULL", significance = 0.05){
    # for genes that are not expressed basemean =0 and show no change in
    # expression set log2FoldChange to 0 as it is NA and causes problems
    if(! 'Gene' %in% colnames(deseq_res)){
        deseq_res = add_rownames(as.data.frame(deseq_res), var = "Gene")
    }
    # Keeping only the genes that are expressed in at least one condition so baseMean won't be 0
    deseq_res = deseq_res %>% filter(!is.na(log2FoldChange)) %>%
                mutate(padj = ifelse(is.na(padj), 1, padj))

    p = ggplot(data = deseq_res, aes(x = log2FoldChange))
    p = p + geom_line(stat = "density", aes(color = padj < significance))
    # Manually setting the xlim to -4 and +4
    p = p + coord_cartesian(xlim = c(-4.5, 4.5)) + scale_color_manual(values = unname(default_colors[c('blue', 'darkred')]))
    p = p + theme_bw() + theme(legend.position= "bottom") + labs(color = paste0('padj < ', significance)) + ggtitle(label)
    plot(p)
    
    # density = counts / sum(counts * bar width)
    if (!is.null(subsets)){
        deseq_res = deseq_res %>% mutate(Subset = Gene %in% subsets)
        q = ggplot(data = deseq_res, aes(x = log2FoldChange))
        q = q + geom_line(stat = "density", aes(color = padj < significance, linetype = Subset)) + labs(linetype = subset_name)
        # Manually setting the xlim to -4 and +4
        q = q + coord_cartesian(xlim = c(-4.5, 4.5)) + scale_color_manual(values = unname(default_colors[c('blue', 'darkred')]))
        q = q + theme_bw() + theme(legend.position= "bottom") + labs(color = paste0('padj < ', significance)) + ggtitle(label)
        
        plot(q)
    }
    
}# }}}

#still working on it
get_expression_quartiles = function(tpms, contrast){
    quartiles = melt(tpms) %>% mutate(Contrast = variable) %>%
        mutate(Contrast = ifelse(Contrast %in% names(contrast),contrast[Contrast], Contrast))

    # not sure when to use the 'valid' quartiles, which is when I remove the not expressed genes, TPM == 0
    valid = quartiles %>% filter(value != 0) %>% mutate(value = log2(1 + value)) %>% 
        group_by(Contrast) %>% mutate(quartiles_valid = ntile(value, 4)) %>% ungroup()
    quartiles = quartiles %>% mutate(value = log2(1 + value)) %>% group_by(variable) %>%
                mutate(quartiles_all = ntile(value, 4)) %>% ungroup()

    quartiles = left_join(quartiles, valid, by = c("Gene", "variable", "value", "Contrast")) %>%
        mutate(quartiles_valid = ifelse(is.na(quartiles_valid), 0, quartiles_valid))

    quartiles %>% group_by(Contrast) %>% head()
    
    return(quartiles)
}



#Make a plotting function# {{{
plotDE<-function(res, fdr){
 df<-data.frame(res, FDR=ifelse(res$padj<fdr, paste("<", fdr, sep=" "), paste(">", fdr, sep=" ")))

 p<-qplot(log(baseMean), log2FoldChange, data=df, aes(x="baseMean", y="log2FoldChange"))
 p<-p + geom_point(aes(colour=FDR))
 print (p)
}
# }}}

#Make a plotting function for introns--take into account exon de analysis# {{{
plotDEnIntron<-function(res, fdr){
 df<-data.frame(res, FDR=ifelse(res$padj<fdr, paste("<", fdr, sep=" "), paste(">", fdr, sep=" ")))

 p<-ggplot(df, aes(x=log(baseMean), y=log2FoldChange, colour=FDR, shape=DEexon, size=2))
 p<-p + geom_point(aes(colour=FDR, size=FDR, shape=DEexon), size=2)
 print (p)
}
# }}}

#Run DESeq and compare 2 conditions# {{{
NbinomTest <- function(cds, cond1, cond2, geneNames, exonDE=NULL) {
 #Run the negative binomial test on the CountDataSet(cds with estimatedDispersions) for the 2 conditions
 res <- nbinomTest( cds, cond1, cond2 )

 if(is.null(exonDE)==FALSE) res<-data.frame(res, DEexon=ifelse(res$id %in% exonDE$id, "yes", "no"))

 #will need to give it a more descriptive name and check if the conditions are printed in the plot
 pdf(paste(cond1, "vs", cond2, ".pdf", sep=""))
 ifelse(is.null(exonDE)==FALSE, plotDEnIntron(res, .05), plotDE(res, .05))
 dev.off()

 #Filter for significant genes for FDR threshold of 5% (also don't take into account genes for which padj was not defined (i.e. NA) because fold change is NaN (i.e. division of 0 by 0 (counts/size/factor-??)))
 resSig <- res[!is.na(res$padj) & res$padj < .05, ]
 #Order the genes based on p-value
 resSig<<-resSig[ order(resSig$pval), ]

 #Get the Gene Name for the e! Gene Ids of the de significant genes-REMEMBER: cannot give remco's file cause he got the genes names probably based on gene ids for MM9
 #mouse_data<-read.delim(id_to_name_file, as.is=TRUE, quote="")
 geneNames<-read.delim(geneNames, as.is=TRUE, quote="")

 #Get the strongly downregulated of the significant genes--this will just order it so the downreg come first--the list still contains the upreg
 downRegSig<-resSig[order(resSig$foldChange, -resSig$baseMean), ] 
 down_mouse_names<-geneNames$mgi_symbol[ match(downRegSig$id,geneNames$ensembl_gene_id) ]
 down_genes<-cbind(downRegSig,Name=down_mouse_names)

 write.table(down_genes, file=paste(cond1,"vs",cond2,"-downRegSig.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

 #Get the strongly upregulated of the significant genes
 upRegSig<-resSig[order(-resSig$foldChange, -resSig$baseMean), ]
 up_mouse_names<-geneNames$mgi_symbol[ match(upRegSig$id,geneNames$ensembl_gene_id) ]
 up_genes<-cbind(upRegSig,Name=up_mouse_names)

 write.table(up_genes, file=paste(cond1,"vs",cond2,"-upRegSig.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

 #Change the match arguments if you have different header--what happens for ids that don't have an mgi symbol?
 de_mouse_names<-geneNames$mgi_symbol[ match(resSig$id,geneNames$ensembl_gene_id) ]
 #DE significant genes (id to name)
 degenes<-cbind(resSig,Name=de_mouse_names)
 #Save an object of the de genes
 save(degenes, file=paste(cond1,"vs",cond2,"-de.Rda",sep=""))  

 #Write output-all de significant genes with gene name instead of gene id
 write.table(degenes, file=paste(cond1,"vs",cond2,"-de.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

 #Write output for all genes ordered by p-val
 resAll<-res[ order(res$pval), ]

 all_mouse_names<-geneNames$mgi_symbol[ match(resAll$id,geneNames$ensembl_gene_id) ]
 allgenes<-cbind(resAll,Name=all_mouse_names)
 #Save an object for all genes 
 save(allgenes, file=paste(cond1,"vs",cond2,"-de.Rda",sep=""))
 
 #Write output for all genes with gene name instead of id
 write.table(allgenes, file=paste(cond1,"vs",cond2,".txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
}# }}}

