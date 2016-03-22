#!/usr/bin/env bash

cd /nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013

declare -A swapped
swapped=(["3KO_h24-S5P_2"]="3KO_h4-Input_3" ["3KO_h4-Input_3"]="3KO_h24-Nanog_3" ["3KO_h24-Klf4_3"]="3KO_h24-S5P_2" ["3KO_h24-Nanog_3"]="3KO_h24-Klf4_3")

mkdir -p mm10/bowtie/swapped
mkdir -p data/swapped

for x in "${!swapped[@]}"
do 
    y=${swapped[${x}]}
    echo "${x} is really ${y}"
    #samtools view -H mm10/bowtie/${x}.bam | sed "s/$x/$y/g" > ${y}.sam_header
    #samtools view mm10/bowtie/${x}.bam > ${y}.sam_main
    #cat ${y}.sam_header ${y}.sam_main > ${y}.sam
    #samtools view -bSo ${y}.bam ${y}.sam
    mv mm10/bowtie/${x}.bam mm10/bowtie/swapped/
    mv ${y}.bam mm10/bowtie/

    cp mm10/bowtie/${x}.stderr mm10/bowtie/swapped/
    cp mm10/bowtie/swapped/${x}.stderr mm10/bowtie/${y}.stderr

    # Saving the old files to a dir as I can not rename in one go
    mv data/${x}.fastq.gz data/swapped/${x}.fastq.gz
done

for x in "${!swapped[@]}"
do 
    # doing this in a second loop as if doing all in one will result in overwriting
    # files that have not been changed yet
    cp data/swapped/${x}.fastq.gz data/${y}.fastq.gz
done
