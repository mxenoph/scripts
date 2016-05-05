bioc_updates = function(url, inst_pkgs){
    # e.g. "http://www.bioconductor.org/packages/release/bioc/"
    bioc= available.packages(contrib.url(url))
    common_pkgs= intersect(rownames(bioc), rownames(inst_pkgs))
    status = mapply(function(x, y) compareVersion(x, y), inst_pkgs[common_pkgs, "Version"], bioc[common_pkgs, "Version"])
    status = status[status < 0]
    return(status)
}
