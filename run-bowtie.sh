#!/bin/bash

#run-bowtie.sh #{{{
usage() { echo "Usage: $0 [-f <fastq>] [-g <genome index directory>] [-o <output directory>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':f:g:o:h'
while getopts $options option
do
    case $option in
        f ) fastq=(${OPTARG}) ;;
        g ) GENOME=(${OPTARG}) ;;
        o ) OUT=${OPTARG} ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

if [ -z "$fastq" ]
then
    usage
fi

if [ -z "$GENOME" ]
then
    GENOME='/nfs/research2/bertone/common/genome/MM10/MM10'
fi

if [ -z "$OUT" ]
then
    OUT=$PWD"/bowtie/"
else
    if [[ ! $OUT =~ \/bowtie\/$ ]]
    then
        OUT="$(readlink -m ${OUT})/bowtie/"
    fi
fi
# }}}

# readlink converts relative to absolute path names
fastq=$(readlink -m ${fastq})
GENOME=$(readlink -m ${GENOME})

mkdir -p ${OUT}

describer=$(sed 's/\(.fastq.*\|.fq.*\|.txt\)$//' <<< $fastq | sed 's/.*\///')

if [[ $fastq =~ ".gz" ]]
then
    gunzip -d $fastq
    fastq=${fastq/.gz/}
    #bowtie -m 1 -S -p 4 $GENOME $fastq > ${OUT}${describer}.sam
    #gunzip $fastq
fi

# TODO: add more options if ever get PE reads for ChIP-seq or if I align RNA-Seq 
# with bowtie instead of gsnap (not gapped alignment)
# -p <n> launches parallel search threads => Always bsub it with -n 4
bowtie -m 1 -S -p 4 $GENOME $fastq > ${OUT}${describer}.sam
# Keeping a record of how the files in the directory where generated
#echo -e `date +"%D%t%T"` "\t" "bowtie -m 1 -S $GENOME $fastq > ${OUT}${describer}.sam" >> ${OUT}FilesTree.log

# TODO: add more peak callers if ever need to use a different one

# SAM to BAM
samtools view -bS ${OUT}${describer}.sam -o ${OUT}${describer}.uns.bam
#echo -e `date +"%D%t%T"` "\t" "samtools view -bS ${OUT}${describer}.sam -o ${OUT}${describer}.uns.bam" >> ${OUT}FilesTree.log

# Sort BAM
samtools sort ${OUT}${describer}.uns.bam ${OUT}${describer}
#echo -e `date +"%D%t%T"` "\t" "samtools sort ${OUT}${describer}.uns.bam ${OUT}${describer}.sort" >> ${OUT}FilesTree.log

# Sort BAM
samtools index ${OUT}${describer}.bam
#echo -e `date +"%D%t%T"` "\t" "samtools index ${OUT}${describer}.sort.bam" >> ${OUT}FilesTree.log

# Regardless of whether the file was compressed to begin with, compress it
gzip $fastq

# Remove intermediate files
rm ${OUT}${describer}.sam
rm ${OUT}${describer}.uns.bam
