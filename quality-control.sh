#!/usr/bin/env bash

function qc_trim () {
    fastq=$1
    qc_files=$2

    encoding=$(sed -n 6p ${qc_files}/fastqc_data.txt | cut -f 2)
    if [[ ${encoding} != "Sanger / Illumina 1.9" ]]
    then
        echo "ERROR:Encoding is not Phred33"
        exit 1
    fi
    phred=33
    
    fastq_quality_trimmer -Q ${phred} -t 20 -l 20 -i ${fastq} -o ${fastq%%.*}_quality_trimmed.fastq
}

function qc () {
    fastq=$1
    path=$(dirname ${fastq})/fastqc
    mkdir -p ${path}
    target=$(basename ${fastq%%.*})
    echo $target
    fastqc -o ${path} ${fastq}

    seq_qual=$(sed -n 2p ${path}/${target}_fastqc/summary.txt | cut -f 1)
    if [[ ${seq_qual} == "FAIL" ]]
    then
        echo $seq_qual
        qc_trim ${fastq} ${path}/${target}_fastqc
        qc ${fastq%%.*}_quality_trimmed.fastq
    fi
    
}

qc $1

