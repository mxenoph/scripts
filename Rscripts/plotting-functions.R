plot_venn <- function(a, b, intersection, label){# {{{
    x <- c('VennDiagram', 'grid', 'RColorBrewer')
    lapply(x, suppressMessages(library), character.only=T)

    pastel_col <- brewer.pal(8, "Pastel2")
    venn <- draw.pairwise.venn(
                               area1 = a,
                               area2 = b,
                               cross.area = intersection,
                               scaled = TRUE,
                               category = label,
                               fill = pastel_col[1:2],
                               aplha = 0.5,
                               lty = "blank",
                               label.col=c("black"),
                               cex = 1.5,
                               cat.cex = 1.5,
                               cat.pos = c(285, 105),
                               cat.dist = 0.09,
                               cat.just = list(c(-1, -1), c(1, 1)),
                               ext.pos = 30,
                               ext.dist = -0.05,
                               ext.length = 0.85,
                               ext.line.lwd = 2,
                               ext.line.lty = "dashed"
                               );
    grid.draw(venn);
}
# comment(plot_venn) will print this message
attr(plot_venn, "comment") <- "plot_venn() expects the size of a, b and their intersection and a vector of length 2 with the names of the 2 categories"
# }}}
