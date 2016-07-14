#!/bin/bash

#run-bowtie.sh #{{{
usage() { echo "Usage: $0 [-f <fastq>] [-g <genome index directory>] [-o <output directory>] [-t threads; default = 4] [-c parent folder of FASTQC generated folder]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':f:g:o:t:c:h'
while getopts $options option
do
    case $option in
        f ) fastq=(${OPTARG}) ;;
        g ) GENOME=(${OPTARG}) ;;
        o ) OUT=${OPTARG} ;;
        t ) THREADS=${OPTARG} ;;
        c ) FASTQC=${OPTARG} ;;
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

if [ -z "$THREADS" ]
then
    THREADS=4
fi

DDUP=${OUT}"ddup/"
# }}}

# readlink converts relative to absolute path names
fastq=$(readlink -m ${fastq})
GENOME=$(readlink -m ${GENOME})

mkdir -p ${OUT}
mkdir -p ${DDUP}

describer=$(sed 's/\(.fastq.*\|.fq.*\|.txt\)$//' <<< $fastq | sed 's/.*\///')

if [[ $fastq =~ ".gz" ]]
then
    gunzip -d $fastq
    fastq=${fastq/.gz/}
fi

fastqc_report=${FASTQC}/$(basename ${fastq/.f*/_fastqc})
if [ -e "${fastqc_report}/fastqc_data.txt" ]
then
    encoding=$(awk -F '\t' '$0 ~ /^Encoding/ {if($2 ~ /Sanger/) {print "--phred33-quals"; exit } else if($2 ~ /Illumina 1.5/) {print "--phred64-quals"; exit} }' ${fastqc_report}/fastqc_data.txt)
else
    # most new data are PHRED+33 so use it as default but check if mapping crashes
    encoding="--phred33-quals"
fi

# TODO: add more options if ever get PE reads for ChIP-seq or if I align RNA-Seq 
# with bowtie instead of gsnap (not gapped alignment)
# -p <n> launches parallel search threads => Always bsub it with -n 4
#bowtie -m 1 -S -p 4 $GENOME $fastq 2> ${OUT}${describer}.stderr 1> ${OUT}${describer}.sam

bowtie ${encoding} -y -m 1 --best --strata -S -p ${THREADS} --nomaqround $GENOME $fastq 2> ${OUT}${describer}.stderr 1> ${OUT}${describer}.sam
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

# Mark duplicates, $picard is an env variable
/ebi/research/software/Linux_x86_64/opt/java/jdk1.6/bin/java -Xmx2g -jar \
    /nfs/research2/bertone/software/picard-tools-1.76/MarkDuplicates.jar \
    INPUT=${OUT}${describer}.bam OUTPUT=${DDUP}${describer}_ddup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true M=${DDUP}${describer}_ddup.log

# Regardless of whether the file was compressed to begin with, compress it
gzip $fastq

# Remove intermediate files
rm ${OUT}${describer}.sam
rm ${OUT}${describer}.uns.bam
