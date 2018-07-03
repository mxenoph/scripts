#!/bin/env bash

gtf=$1

# Scan every 50 lines starting with line=1. If number found then print 1 and exit
check=$(awk 'NR%50==1 {if($1 ~ "^[0-9][0-9]*$" || $1 ~ "^X|Y|MT$") {print 1; exit;}}' ${gtf})

# If ${check} is undefined then it will evaluate to true
if [ -z ${check} ]
then
    echo "Chromosomes start with chr in GTF already."
else
    awk 'FS=OFS="\t" {if ($1 ~ "^[0-9][0-9]*$" || $1 ~ "^X|Y$") {$1="chr"$1}; if ($1 ~ "^MT") {$1="chrM"} print $0}' ${gtf} > ${gtf/.gtf/.tmp} && mv ${gtf/.gtf/.tmp} ${gtf}
fi
