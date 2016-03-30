# Get ggplot colours# {{{
gg_color_hue <- function(n) {
    x <- c('ggplot2')
    mclapply(x, suppressMessages(library), character.only=T)
    
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}
# }}}

# Get ggplot legend
g_legend<-function(p){# {{{
    x <- c('ggplot2')
    mclapply(x, suppressMessages(library), character.only=T)

    tmp <- ggplotGrob(p)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
# }}}