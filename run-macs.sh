#!/usr/bin/env bash

# Parse arguments#{{{
ARGS=$(getopt -o c:t:g: -l "control:,treatment:,genome:" -n "run-macs.sh" -- "$@")

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
            genome_esize="$2"; shift 2 ;;
        *)
            echo "Invalid option!"; exit 1 ;;
    esac
done#}}}

target_dir="$(dirname "$ctrl")"/macs
target_base="$(basename "$ctrl")"
output="${target_base/.bam//}"

mkdir -p "$target_dir"
# MACS only save output in working directory
cd "$target_dir"

macs14 \
    -t "$treat" -c "$ctr" \
    --gsize="$genome_esize" --name="$output" \
    --format=BAM -S --diag \
