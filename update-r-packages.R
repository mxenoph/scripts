#!/usr/bin/env Rscript --slave

pkgs <- installed.packages()
r-version <- paste(version[['major']], '.', version[['minor']], sep = '')
 
for (i in 1:nrow(pkgs))
{
    pkg-name <- pkgs[i, 1]
    pkg-version <- pkgs[i, 3]
    if (pkg-version != r-version)
    {
        print(paste('Installing', pkg-name))
        install.packages(pkg-name)
    }
}
