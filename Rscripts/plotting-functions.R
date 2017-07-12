library = function (...) suppressMessages(base::library(...))

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

plot_venn = function(a, b, intersection, label, cols = NULL){# {{{
    x = c('VennDiagram', 'grid', 'RColorBrewer')
    lapply(x, suppressMessages(library), character.only=T)

    if(is.null(cols)){
        pastel_cols = default_colors[1:2]
    } else {
        pastel_cols =  unname(default_colors[cols])
    }
    venn = draw.pairwise.venn(
                               area1 = a,
                               area2 = b,
                               cross.area = intersection,
                               scaled = TRUE,
                               category = label,
                               fill = pastel_cols,
                               aplha = 0.5,
                               label.col= c("black",  "white", "black"),
                               lty = "blank",
                               cex = 2.5,
                               cat.cex = 2.5,
                               cat.pos = c(330,155),
                               cat.dist = 0.05,
                               cat.default.pos='outer',
                               ext.pos = 30,
                               ext.dist = -0.05,
                               ext.length = 0.85,
                               ext.line.lwd = 2,
                               fontface = c('plain', 'plain', 'plain'),
                               catface = c('plain', 'plain', 'plain'),
                               fontfamily = c('sans', 'sans', 'sans'),
                               catfamily = c('sans', 'sans', 'sans'),
                               ext.line.lty = "dashed"
                               #cat.just = list(c(-1, -1), c(1, 1)),
                               );
    grid.draw(venn);
}
# comment(plot_venn) will print this message
attr(plot_venn, "comment") = "plot_venn() expects the size of a, b and their intersection and a vector of length 2 with the names of the 2 categories"
# }}}

drawVenn = function(sets){# {{{
    x = c('Vennerable', 'RColorBrewer')
    lapply(x, suppressMessages(library), character.only=T)
    
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

# venn counts to venn object from vennerable # {{{
venn_counts2venn <- function(venn_cnt){
    x <- c('Vennerable', 'RColorBrewer', 'stringr')
    lapply(x, suppressMessages(library), character.only=T)

    n = which(colnames(venn_cnt) =="Counts") - 1
    set_names = colnames(venn_cnt)[1:n]
    weight = venn_cnt[,"Counts"]
    names(weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")

    vn = Venn(SetNames = set_names, Weight = weight)
    c_vn = compute.Venn(vn)
    theme_vn = VennThemes(c_vn)

    cols = c("#71E5A5", "#F1A5AD", "#E4B869", "#B7D98C",
             "#D3E164", "#6CDAE9", "#8DDAC4", "#BAC7DF", "#DBB1E6")
    shades = c("#795F0F", "#8F5B16", "#A55E17", "#8B7239", "#B65A24", 
               "#A06F4D", "#C9920E", "#A58530", "#A97046", "#B38926",
               "#CE712A", "#D38113", "#A88C5B", "#DC7733", "#A08A61", 
               "#EF9B18", "#CAA03A", "#D38655", "#F9C804", "#CF9462", 
               "#B49170", "#FABF16", "#F3A526", "#C1A463", "#FFAC4E", "#E8C36D",
               "#F9B276", "#F8CE63", "#FCC270", "#DAC090")

    for (i in 1:length(theme_vn[["Set"]])) {
        theme_vn[["Set"]][[i]]$col = cols[i]
        theme_vn[["SetText"]][[i]]$col = cols[i]
    }
    for (i in 1:length(theme_vn[["SetText"]])) {
        theme_vn[["SetText"]][[i]]$fontsize = 28
    }
    for (i in 1:length(theme_vn[["FaceText"]])) {
        theme_vn[["Face"]][[i]]$col = shades[i]
        theme_vn[["Face"]][[i]]$fill = shades[i]
        theme_vn[["FaceText"]][[i]]$col = colors()[295]
        theme_vn[["FaceText"]][[i]]$fontsize = 28
    }

    sets_faces = names(theme_vn[['Face']])[str_count(names(theme_vn[['Face']]), '1') == 1]
    sets_faces = sets_faces[order(sets_faces, decreasing = T)]
    for (x in sets_faces){
        index = grep(x, sets_faces)
        theme_vn[["Face"]][[x]]$col = theme_vn[['Set']][[index]]$col
        theme_vn[["Face"]][[x]]$fill = theme_vn[['Set']][[index]]$col
        theme_vn[["Face"]][[x]]$fontsize = 28
        theme_vn[["FaceText"]][[x]]$fontsize = 28
    }
    return(list(vn = vn, theme = theme_vn, c_vn = c_vn))
}
# }}}

flag_replicates = function(fs){# {{{
    x = c('ggplot2', 'dplyr')
    lapply(x, suppressMessages(library), character.only=T)
    
    df = read.table(fs, header=F, sep="\t", col.names=c('comparison', 'npeaks'))

    df = df %>% 
            mutate(comparison = gsub("-overlapped-peaks.txt", '', basename(as.character(comparison))),
                   group = gsub("_.vs.*", '',basename(as.character(comparison)))
                   )

    p = ggplot(df, aes(x=comparison, y=npeaks, fill=group)) + geom_bar(stat='identity', width=0.5)
    p = p + facet_grid(. ~group, scales="free")
    # removes the space betwwen 0 and x axis
    p = p + scale_y_continuous(expand=c(0,0)) + scale_x_discrete(expand=c(0.02,0))
    # clean look from http://rstudio-pubs-static.s3.amazonaws.com/3364_d1a578f521174152b46b19d0c83cbe7e.html
    # and vertical labels
    p = p + theme_classic() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
                                    axis.line.x=element_blank(),
                                    strip.background = element_rect(colour="white", fill="white"),
                                    legend.position="none")

    return(p)
}# }}}


# Plot FPKMs or TPMs for developmental genes# {{{
plot_tpm <- function(tpm, design, gtf, markers, significant, name, what = 'p2'){
    x = c('ggplot2', 'gridExtra', 'dplyr', 'reshape2', 'RColorBrewer')
    lapply(x, suppressMessages(library), character.only=T)
    source("~/source/Rscripts/ggplot-functions.R")
    select = dplyr::select
    
    add_rownames = function(df, var = 'Gene') {
        stopifnot(is.data.frame(df))
        rowname_df = setNames(data_frame(rownames(df)), var)
        cbind(rowname_df, df)
    }
    
    mark_de = function(x, gene_universe) {
        per_condition = significant %>% filter(Condition == x)
        gene_universe[[x]] = as.numeric(gene_universe[['Gene']] %in% per_condition[['Gene']])
        return(gene_universe)
    }

    tpm = add_rownames(tpm)
    # melting dataframe on gene name and conditions for which we have de info
    # Renaming variable to Library so it matches design and left_join can be used
    tpm = melt(tpm,
               id.vars = 'Gene',
#               id.vars = c('Gene', ((significant %>% select(Condition) %>% unique() %>% .[[1]]))),
               variable_name = "Library")

    tpm = full_join(tpm, design, by = "Library")
    # for each condition present in the significant dataframe check if gene de
    for(x in ((significant %>% select(Condition) %>% unique() %>% .[[1]]))) {
        tpm = tpm  %>% mutate(Group = ifelse(Gene %in% significant[significant$Condition == x, 'Gene'],
                                             'p < 0.05', 'p > 0.05'))
    }

    # Import gtf and read gene names or load genes with get_annotation(assembly)
    gtf = import(gtf)
    gene_names = as.data.frame(values(gtf)) %>% select(gene_id, gene_name) %>% unique()
    colnames(gene_names) = gsub('gene_id', 'Gene', colnames(gene_names))
    tpm = full_join(tpm, gene_names, by = "Gene")

    gg = brewer.pal(n = 6, name = "Paired")
    gg = gg[c(1:2,5:6)]
    names(gg) = rep(unique(design[['Contrast']]), each = 2)
    names(gg)[c(2,4)] = gsub("$", "-de", names(gg)[c(2,4)])


    pdf('testnew.pdf')
        p = ggplot(tpm, aes(x = Gene, y = value, color=Group)) + geom_point(size=2) + scale_color_manual(values=gg)
        p = p + geom_errorbar(stat = "hline", yintercept = "mean", width = 0.8, aes(ymax = ..y.., ymin = ..y..))
        p = p + facet_grid(Condition ~ State, scales="free")
        p = p + scale_y_continuous(breaks= c(seq(0, 100, by=20), seq(100, 800, by=200)))
        p = p + theme_bw() + xlab('') + ylab("TPM")
        p = p + theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0),
                       axis.text.y = element_text(size=7, face="plain"),
                       panel.grid.major.x = element_blank(),
                       axis.text.x=element_text(angle=90, vjust=1, face="plain")) + ggtitle(name)
        print(p)
    dev.off()

    for (x in levels(factor(tpm$Condition))){
        tpm[tpm$Condition == x & tpm$X1 %in% names(de[[x]][de[[x]]==TRUE]), 'Group'] <- paste0(tpm[tpm$Condition == x & tpm$X1 %in% names(de[[x]][de[[x]]==TRUE]), 'Group'], '_de')
    }
    tpm[['X1']] <- unlist(mclapply(as.character(tpm[['X1']]), function(x) names(markers[markers == x])))

    group_color <- tpm$Group
    group_color <- as.character(gg[match(group_color,names(gg))])
    tpm$Group <- factor(tpm$Group, levels=c('WT','KO', 'WT_de', 'KO_de'))

    if(what == 'p1'){
        p <- ggplot(tpm, aes(x=Condition, y=value, color=Group)) + geom_point() + scale_color_manual(values=gg)
        p <- p + geom_errorbar(stat = "hline", yintercept = "mean", width=0.8,aes(ymax=..y..,ymin=..y..))
        p <- p + facet_grid(~ X1,) + theme_bw() + xlab('') + ylab("FPKM")
        p <- p + theme(panel.grid.major.x= element_blank(),
                       strip.text.x = element_text(size = 10, colour = "black", angle = 90),
                       axis.text.x=element_text(angle=90, vjust=1)) + ggtitle(name)
        print(p)
   } else if (what == 'p2') {
        tpm[as.character(tpm[['X1']]) %in% names(markers)[1:7], 'State'] <- 'Pluripotency'
        tpm[as.character(tpm[['X1']]) %in% names(markers)[8:13], 'State'] <- 'Endoderm/Ectoderm'
#        tpm[as.character(tpm[['X1']]) %in% names(markers)[1:3], 'State'] <- 'Pluripotency'
#        tpm[as.character(tpm[['X1']]) %in% names(markers)[4:7], 'State'] <- 'Naive Pl.'
#        tpm[as.character(tpm[['X1']]) %in% names(markers)[8:10], 'State'] <- 'Primed Pl.'
#        tpm[as.character(tpm[['X1']]) %in% names(markers)[11:15], 'State'] <- 'Endoderm'
#        tpm[as.character(tpm[['X1']]) %in% names(markers)[c(16:17,19)], 'State'] <- 'Ectoderm'
#        tpm[as.character(tpm[['X1']]) %in% names(markers)[18], 'State'] <- 'Mesoderm'
#        tpm[as.character(tpm[['X1']]) %in% names(markers)[c(20,21)], 'State'] <- 'JAK/STAT'
        
        p <- ggplot(tpm, aes(x=X1, y=value, color=Group)) + geom_point(size=2) + scale_color_manual(values=gg)
        p <- p + geom_errorbar(stat = "hline", yintercept = "mean", width=0.8,aes(ymax=..y..,ymin=..y..))
        p <- p + facet_grid(Condition ~ State, scales="free")
        p <- p + scale_y_continuous(breaks= c(seq(0, 100, by=20), seq(100, 800, by=200)))
        p <- p + theme_bw() + xlab('') + ylab("TPM")
        p <- p + theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0),
                       axis.text.y = element_text(size=7, face="plain"),
                       panel.grid.major.x = element_blank(),
                       axis.text.x=element_text(angle=90, vjust=1, face="plain")) + ggtitle(name)
        print(p)
    }

}# }}}
