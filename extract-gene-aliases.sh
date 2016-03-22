#!/usr/bin/env bash

#9606 for human, 10090 for mouse
#species=$1

gzip -cd gene_info.gz | awk '$1==10090{print $3"\t"$5}' > gene-aliases.txt
