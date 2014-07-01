#!/bin/bash

#TODO: set the files to a variable but figure out how to format ls output

while getopts ":c:t:o:a:" opt; do
    case $opt in
        c)
            ctrl=${OPTARG} ;;
        t)
            treat=${OPTARG} ;;
        o)
            out=$(sed 's/\$//' <<< ${OPTARG}) ;;
        a)
            genome_size=${OPTARG} ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
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
    cmd=`echo "bsub -n 4 -M 20000 -J MACS macs14 -t $t -c $ctrl --format=BAM -S --name=$n --gsize=$genome_size --diag --wig"`;
    echo $cmd;
#    $cmd;
done
