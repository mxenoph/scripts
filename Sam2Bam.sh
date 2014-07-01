#!/bin/bash

# Get command line arguments -- list of directories
if [ "$#" -eq 0 ]
then
    dir=( $PWD )
else
    dir=( "$@" )
fi

for d in "${dir[@]}"
do
    echo "${d}"
    for f in `ls $d/*sam`
    do
        out=`sed 's/.sam/.bam/' <<< $f`;
        cmd=`echo "bsub -n 4 -M 20000 -J Sam2Bam samtools view -bS $f -o $out"`;
        $cmd;
    done
done
