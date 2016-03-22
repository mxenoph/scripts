#!/bin/bash

#Always run with bsub and give memory and ncores to run faster

file=$1
chromSize=$2
out=$(echo $file| sed 's/\.gtf//')

cat $file | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4,$5}' > $out.exon.bed
bedtools sort -i $out.exon.bed > $out.exon.temp.bed
mv -f $out.exon.temp.bed $out.exon.bed
bedtools merge -i $out.exon.bed > $out.exon.merged.bed

cat $file | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4,$5}' > $out.gene.bed
bedtools sort -i $out.gene.bed > $out.gene.temp.bed
mv -f $out.gene.temp.bed $out.gene.bed
bedtools subtract -a $out.gene.bed -b $out.exon.merged.bed > $out.intron.bed

bedtools complement -i $out.gene.bed -g $chromSize > $out.intergenic.bed

