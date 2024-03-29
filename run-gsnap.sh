#!/usr/bin/env bash

# Notes: 
# input files must have extension .fastq.gz
# example:
# for i in {0..19}; do run-gsnap.sh -f "/nfs/research2/bertone/user/mxenoph/hendrich/QC/C1UY6ACXX_2Lox-Epi_13s000232-1-1_Hendrich_lane313s000232_1_sequence.fastqc.gz /nfs/research2/bertone/user/mxenoph/hendrich/QC/C1UY6ACXX_2Lox-Epi_13s000232-1-1_Hendrich_lane313s000232_2_sequence.fastq.gz"
# -g mm10 -p 20 -x $i -s 2lox-Epi -t $PWD; done


# Parse arguments#{{{
ARGS=$(getopt -o f:g:x:p:s:t: -l "fastq:,genome:index:nparallel:sample:target:" -n "run-gsnap.sh" -- "$@")

# Bad arguments
if [ $? -ne 0 ]
then
    exit 1
fi
eval set -- "$ARGS"

while true
do
    case "$1" in
        -f | --fastq)
            #fastq="$2"; shift 2 ;;
            fastq=($2); shift 2 ;; # Read it as an array
        -g | --genome)
            genome="$2"; shift 2 ;;
        -p | --nparallel)
            N_JOBS="$2"; shift 2 ;;
        -x | --index)
            INDEX="$2"; shift 2 ;;
        -s | --sample)
            SAMPLE="$2"; shift 2 ;;
        -t | --target)
            TARGET="$2"; shift 2 ;;
        --)
            shift ; break ;;
        *) echo "Error! Invalid option provided"; exit 1 ;;
    esac
done
#}}}

TARGET="$TARGET"/gsnap
mkdir -p "$TARGET"
output="${TARGET}/${SAMPLE}.$INDEX"

# Check-point for PE data#{{{
if [[ ${#fastq[@]} -eq 2 ]]
then
    # Sanity check that the fastq provided are for mate reads
    if [[ "${fastq[0]/_1_sequence*/}" != "${fastq[1]/_2_sequence*/}" ]]
    then
        echo 'Fastq files provided are either not for mate reads or not in mate_1 mate_2 order.'
        exit 1;
    fi
fi
#}}}

# Genome directories and splice-sites paths#{{{
genome_dir=/nfs/research2/bertone/common/genome
splice_dir=/nfs/research2/bertone/user/mxenoph/genome_dir
#}}}

# GSNAP Options Explanation:#{{{
# -m: number of mismatches (for 75bp set to 4 or 5);
# -i: indel penalty (default=2)...<2 can lead to FP at the end of the read
# -N: look for novel splicing
# -w: for novel splicing event - so it doesn't consider a splice junction that is too big (runs from one end of the chromosome to the other)
# -n: number of paths to print
# -Q: print nothing if >n paths
# -O: relevant if more than one worker thread (-t) -> prints output in the same order as input#}}}

# Setting GSNAP OPTIONS#{{{
N_THREADS=10
# kmer size to use in genome database -- if not specified highest available in the db is used
#K_MER=15
MISMATCH=7
PROTOCOL='sanger'

# ${genome^^} converts string to uppercase
if [ ${genome^^} == 'MM9' ]
then
    SPLICE_SITES=${splice_dir}/M_musculus_9/MM9.maps/mm9_ensembl67_splicesites.iit
elif [ ${genome^^} == 'MM10' ]
then
    SPLICE_SITES=${splice_dir}/M_musculus_10/MM10.maps/mm10_ensembl70_splicesites.iit
else
    echo 'No splice sites provided. Exiting.'
    exit 1;
fi


OPT="-q $INDEX/$N_JOBS -t $N_THREADS"
OPT="$OPT -D ${genome_dir}/${genome^^} -d ${genome^^}"
#Option -B is batch mode to make it run faster| (-i) --indel-penalty by default is 2, what is the benefit of setting it to mismatches? --trim-mismatch-score=0 turns off trimming at ends and will give FP at the ends of reads
OPT="$OPT -B 5 -m $MISMATCH -i 1"
#Find novel splicing and limit the size of splice junctions, define distant-splice penalty (-E) for when the intron length exits the value of -w
OPT="$OPT -N 1 -w 100000 -E 100 -s $SPLICE_SITES"
#Number of paths to print and 'don't print failed alignments', name output based on the q option split up
OPT="$OPT -n 10 -Q -O --nofails -A sam --split-output=$output --gunzip --quality-protocol=$PROTOCOL"
#}}}

gsnap=/homes/kostadim/myrto/software/external/gmap-current/bin/gsnap
# Run gsnap
$gsnap $OPT ${fastq[*]}
