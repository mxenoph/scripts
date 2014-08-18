# Get ggplot colours# {{{
gg_color_hue <- function(n) {
    library(ggplot2)
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}# }}}
