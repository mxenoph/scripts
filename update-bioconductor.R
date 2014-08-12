bioc-updates <- function(url){
    bioc= available.packages(contrib.url(url))
    inst-pkgs= installed.packages()
    common-pkgs= intersect(rownames(bioc), rownames(inst-pkgs))
    mapply(function(x, y) compareVersion(x, y), inst-pkgs[common-pkgs, "Version"], bioc[common-pkgs, "Version"])
}
