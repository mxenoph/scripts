bioc_updates = function(url, old_pkgs){
    # e.g. "http://www.bioconductor.org/packages/release/bioc/"
    bioc = available.packages(contrib.url(url))
    get_status = function(pkgs){
        common_pkgs= intersect(rownames(bioc), rownames(pkgs))
        status = mapply(function(x, y) compareVersion(x, y), pkgs[common_pkgs, "Version"], bioc[common_pkgs, "Version"])
        status = status[status < 0]
        return(status)
    }
    old_status = get_status(old_pkgs)
    new_pkgs = installed.packages()
    new_status = get_status(new_pkgs)
}
