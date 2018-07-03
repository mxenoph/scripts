#!/bin/bash

#run-bowtie.sh #{{{
usage() { echo "Usage: $0 [-b <bed>] [-o <output directory>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':b:o:h'
while getopts $options option
do
    case $option in
        b ) fastq=(${OPTARG}) ;;
        o ) OUT=${OPTARG} ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

if [ -z "$bed" ]
then
    usage
fi

if [ -z "$OUT" ]
then
    OUT=$PWD"/chromHMM/"
else
    if [[ ! $OUT =~ \/chromHMM\/$ ]]
    then
        OUT="$(readlink -m ${OUT})/chromHMM/"
    fi
fi
# }}}

# readlink converts relative to absolute path names
bed=$(readlink -m ${bed})

mkdir -p ${OUT}

base_cmd="java -mx4000M -jar ChromHMM.jar"


