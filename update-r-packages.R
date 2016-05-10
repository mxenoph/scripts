#!/usr/bin/env Rscript --slave

source("~/source/update-bioconductor.R")

r_version = paste(version[['major']], '.', version[['minor']], sep = '')
prev_r_version = paste(version[['major']], '.', as.character(as.numeric(version[['minor']]) - 0.1), sep = '')
search_lib = file.path(dirname(.libPaths()[grep(r_version, .libPaths())]), paste0('R-',prev_r_version))
pkgs = installed.packages(lib.loc = search_lib)

bioc_pkgs = bioc_updates('http://www.bioconductor.org/packages/release/bioc/', pkgs)
source("https://bioconductor.org/biocLite.R")
print('Installing bioconductor packages')
biocLite(names(bioc_pkgs))

pkgs = pkgs[!rownames(pkgs) %in% names(bioc_pkgs)]
for (i in 1:nrow(pkgs))
{
    pkg_name = pkgs[i, 1]
    pkg_version = pkgs[i, 'Built']
    if (pkg_version != r_version)
    {
        #print(paste('Installing', pkg_name))
        if(pkg_name == 'colorout'){
            print('Installing colorout')
            devtools::install_github("jalvesaq/colorout")
        } else if (pkg_name == 'modules'){
            print('Installing modules')
            devtools::install_github('klmr/modules')
        } else{
            install.packages(pkg_name)
        }
    }
}
