
# Get ggplot colours# {{{
gg_color_hue <- function(n) {
    x <- c('ggplot2')
    lapply(x, suppressMessages(library), character.only=T)
    
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}
# }}}

# Get ggplot legend
g_legend<-function(p){# {{{
    x <- c('ggplot2')
    lapply(x, suppressMessages(library), character.only=T)

    tmp <- ggplotGrob(p)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
# }}}

# pots in grid.arrange and shared legend plotted separately# {{{
# grid_arrange_shared_legend(p1,p2,p3)
grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
                 do.call(arrangeGrob, lapply(plots, function(x)
                                             x + theme(legend.position="none"))),
                 legend,
                 ncol = 2,
                 heights = unit.c(unit(1, "npc") - lheight, lheight))
}# }}}

