#!/usr/bin/env bash

# Parse arguments#{{{
ARGS=$(getopt -o c:t:g:o: -l "control:,treatment:,genome:,output:" -n "run-macs.sh" -- "$@")

# Bad arguments
if [ $? -ne 0 ]
then
    exit 1
fi
eval set -- "$ARGS"

while true
do
    case "$1" in
        -c | --control)
            ctrl="$2"; shift 2 ;;
        -t | --treatment)
            treat="$2"; shift 2 ;;
        -g | --genome)
            genome_esize="$2"; shift 2 ;; #e.g. mm
        -o | --output)
            target_dir="$2"/macs; shift 2 ;;
        --)
            shift ; break ;;
        *) 
            echo "Error! Invalid option provided"; exit 1 ;;
    esac
done #}}}

target_base="$(basename "$treat")"
output="${target_base/.bam/}"

mkdir -p "$target_dir"
# MACS only save output in working directory
cd "$target_dir"

if [[ $ctrl == '' ]]
then
    echo "Calling peaks without a control experiment!"
    macs14 \
        -t "$treat"\
        --gsize="$genome_esize" --name="$output" \
        --format=BAM -S --diag
else
    macs14 \
        -t "$treat" -c "$ctrl" \
        --gsize="$genome_esize" --name="$output" \
        --format=BAM -S --diag
fi

