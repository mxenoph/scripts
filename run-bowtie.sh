#!/bin/bash

# Example usage: bsub -n 4 -M 20000 run-bowtie.sh fastq genome_index mapper
# e.g. genome = /nfs/research2/bertone/common/genome/MM9/MM9 
# Check if arguments are supplied
if [ $# -eq 0 ]
then
    echo $0 'fastq genome mapper'
fi

# readlink converts relative to absolute path names
fastq=$(readlink -m $1)
genome=$(readlink -m $2)

# Check if an argument was supplied for the mapper to use
# -z checks if the argument is a NULL string and if so executed the body
if [ -z "$3" ]
then
    # Set bowtie as the default mapper if an argument is not provided
    mapper='bowtie'
else
    mapper=$3
fi

out=$PWD'/'$mapper'/'

if [ ! -d $out ]
then
     mkdir $out
     echo -e "Creating $out"
fi

if [ $mapper == 'bowtie' ]
then
    describer=$(sed 's/\(.fastq.*\|.fq.*\|.txt\)$//' <<< $fastq | sed 's/.*\///')

    if [[ $fastq =~ ".gz" ]]
    then
        gunzip -d $fastq
        fastq=${fastq/.gz/}
        bowtie -m 1 -S -p 4 $genome $fastq > ${out}${describer}.sam
        gunzip $fastq
    fi
    # TODO: add more options if ever get PE reads for ChIP-seq or if I align RNA-Seq 
    # with bowtie instead of gsnap (not gapped alignment)
    # -p <n> launches parallel search threads => Always bsub it with -n 4
    bowtie -m 1 -S -p 4 $genome $fastq > ${out}${describer}.sam
    # Keeping a record of how the files in the directory where generated
    echo -e `date +"%D%t%T"` "\t" "bowtie -m 1 -S $genome $fastq > ${out}${describer}.sam" >> ${out}FilesTree.log

# TODO: add more peak callers if ever need to use a different one
fi

# SAM to BAM
samtools view -bS ${out}${describer}.sam -o ${out}${describer}.uns.bam
echo -e `date +"%D%t%T"` "\t" "samtools view -bS ${out}${describer}.sam -o ${out}${describer}.uns.bam" >> ${out}FilesTree.log

# Sort BAM
samtools sort ${out}${describer}.uns.bam ${out}${describer}.sort
echo -e `date +"%D%t%T"` "\t" "samtools sort ${out}${describer}.uns.bam ${out}${describer}.sort" >> ${out}FilesTree.log

# Sort BAM
samtools index ${out}${describer}.sort.bam
echo -e `date +"%D%t%T"` "\t" "samtools index ${out}${describer}.sort.bam" >> ${out}FilesTree.log

# Remove intermediate files
#rm ${out}${describer}.sam
#rm ${out}${describer}.uns.bam
