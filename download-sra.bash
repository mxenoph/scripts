#!/usr/bin/env bash

run_info=$1

while read line
do
    srr=$(cut -d ',' -f 1 <<< $line)
    first_digits=${srr:0:6}
    bsub wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${first_digits}/${srr}/${srr}.sra
done < ${run_info}

