#!/usr/bin/env bash

usage() { echo "Usage: $0 [-a <gtf>] [-b <sorted bam>] [-f <feature (default=exon)>] [-o <output path>] [-p <parameters (default= -m union -q)>] [-s <sample name (default=basename)>] [-l <boolean; TRUE for pe library>]" 1>&2; exit 1; }
# : means takes an argument but not mandatory (if mandatory will have to check after)
LIBRARY=false

options=':a:b:f:o:p:s:l'
while getopts $options option
do
    case $option in
        a ) gtf=(${OPTARG}) ;;
        b ) bam=${OPTARG} ;;
        f ) FEATURE=${OPTARG} ;;
        o ) OUT=${OPTARG} ;;
        #need to remove the param option
        p ) PARAM=${OPTARG} ;;
        s ) SAMPLE=${OPTARG} ;;
        l ) LIBRARY=true ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

# Check that mandatory arguments were provided
# -z gives true if string is empty
if [ -z "$gtf" ] || [ -z "$bam" ]
then
    usage
fi


# Controlling output directory#{{{
if [ ! -z "$OUT" ]
then
    if [[ ! $OUT =~ \/counts$ ]]
    then 
        OUT="$(readlink -m ${OUT})/counts"
        mkdir -p $OUT
    fi
else
    OUT="$(dirname ${bam})/counts"
fi
#}}}

options='-m union -q -f bam'
# Controlling default parameters#{{{
if [ ! -z "$PARAM" ]
then
    options="$options -s $PARAM"
fi

if [ ${LIBRARY} == true ]
then
    options="$options -r name"
else
    options="$options -r pos"
fi

# -n gives true if string is not empty
if [ -n "$FEATURE" ] && [ "$FEATURE" != "exon" ]
then
    FEATURE="-t ${FEATURE} -i ${FEATURE}_id"
else
    FEATURE='-t exon -i gene_id'
fi
#}}}

#PARAM="$PARAM $FEATURE"
options="$options $FEATURE"

# Controlling output name#{{{
if [ ! -n "$SAMPLE" ]
then
    target="${OUT}/$(basename ${bam/.bam/.counts})"
else
    target="${OUT}/${SAMPLE}.counts"
fi
#}}}

bam="$(readlink -m ${bam})"

# HTSeq-count Options:#{{{
# -m <mode> either uniion, intersection strict, intersection-noempty; default:union
# -o <samout_name>: write all SAM alignment records into an output SAM file called samout_name
# -q suppress progress report and warnings
# -s whether data from strand specific assay; default: yes
# -t feature type- default is exon
# For enhancers use intersection-nonempty
# opt="-m intersection-nonempty -q -i enhancer_id"#}}}

htseq-count \
    ${options} \
    ${bam} ${gtf} \
    > ${target}

# will be printed in the log
echo "htseq-count ${options} ${bam} ${gtf} > ${target}"
