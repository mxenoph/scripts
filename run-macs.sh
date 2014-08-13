#!/usr/bin/env bash

#TODO: set the files to a variable but figure out how to format ls output
# parse arguments
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
            genome_size="$2"; shift 2 ;;
        *)
            echo "Invalid option!"; exit 1 ;;
    esac
done

if [ ! -d "$out/macs" ]
then
    mkdir "$out/macs"
fi
cd "$out/macs"

IFS=';' read -a arr <<< "$treat"

for t in "${arr[@]}"
do
    n=$(sed 's/\.bam//' <<< $t)
    cmd=`echo "bsub -n 4 -M 20000 -J MACS macs14 -t $t -c $ctrl --format=BAM -S --name=$n --gsize=$genome_size --diag --wig"`;
#    echo $cmd;
    $cmd;
done
