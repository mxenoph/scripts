#!/usr/bin/env bash 

for p in Nanog Klf4; do
    tmp=$(ls );
    bsub -M 64000 -R "rusage[mem=30000]" -n 8 "computeMatrix reference-point -S mm10/bam-compare/read-count/3KO_h0-${p}_pooled_RPKM.bw mm10/bam-compare/read-count/3KO_30m-${p}_pooled_RPKM.bw mm10/bam-compare/read-count/3KO_h4-${p}_pooled_RPKM.bw mm10/bam-compare/read-count/3KO_h24-${p}_pooled_RPKM.bw -R ${tmp} --outFileName mm10/bam-compare/read-count/matrices/${p}_pooled_RPKM.deeptools-tf-all-peaks-matrix --outFileSortedRegions mm10/bam-compare/read-count/matrices/${p}_pooled_RPKM.deeptools-tf-all-peaks-matrix-sorted-by.bed -bs 10 -p 8  --beforeRegionStartLength 2000 --afterRegionStartLength 2000";
done
