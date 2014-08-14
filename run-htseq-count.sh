#!/usr/bin/env bash

# Parse arguments#{{{
ARGS=$(getopt -o a:b:o:t:: -l "annotation:,bam:,options::,type::" -n "run-htseq-count.sh" -- "$@")

# Bad arguments
if [ $? -ne 0 ]
then
    exit 1
fi
eval set -- "$ARGS"

while true
do
    case "$1" in
        -a | --annotation)
            gtf="$2"; shift 2 ;;
        -b | --bam)
            bam="$2" shift 2 ;;
        -o | --options)
            case "$2" in
                "")
                    opt='-m union -q'; shift 2 ;;
                *)
                    opt="$2"; shift 2 ;;
            esac ;;
        -t | --type)
            case "$2" in
                "")
                    feature='exon'; feature_id='gene_id'; shift 2 ;;
                *)
                    feature="$2"; feature_id=$feature"_id"; shift 2 ;;
            esac ;;
        --)
            shift ; break ;;
        *) echo "Error! Invalid option provided"; exit 1 ;;
    esac
done #}}}

target_dir="$(dirname "$bam")/htseq"
target_base="$(basename "$bam")"

output="$target_dir/${target_base/.bam/.counts}"
if [ $feature != "exon" ]
then
    output="$target_dir/${target_base/.bam/.${feature}.counts}"
fi

mkdir -p "$target_dir"

# HTSeq-count Options:#{{{
# -m <mode> either uniion, intersection strict, intersection-noempty; default:union
# -o <samout_name>: write all SAM alignment records into an output SAM file called samout_name
# -q suppress progress report and warnings
# -s whether data from strand specific assay; default: yes
# -t feature type- default is exon
# For enhancers use intersection-nonempty
# opt="-m intersection-nonempty -q -i enhancer_id"#}}}

htseq-count \
    "$opt" \
    -t "$feature" -i "$feature_id" \
    "$bam" "$annotation" \
    > "$output"

