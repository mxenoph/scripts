#!/bin/bash


USAGE="Usage: $0 lisOfBamFiles.txt scalefactors.txt genome(mm9 or mm10)"
if [ $# -gt 3 ]; then
    echo $USAGE
    exit 1
fi

listOfFiles=$1
scalefactors=$2
genome=$3

while read line
do
 prefix=$(sed 's/.*\///' <<< $line | sed 's/\..*//')
 JOB="$prefix.sam2bw"
 LOG="$prefix.sam2bw.log"

 if [ -z "$scalefactors" ]; then
  #Need to edit bam2bw.sh and add if statement in case I want to run the script without scaling | not so urgent though
  exit 0
 fi

 if [ $genome == "mm10" ]; then
    submit=$(bsub -N -n 6 -M 40000 -R "rusage[mem=10000]" -R "select[ncores=2]"  -J $JOB -o $LOG "bam2bw.sh $line $prefix /nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10_noHead.genome '-f $scalefactors $prefix'")
 fi
 if [ $genome == "mm9" ]; then
    #Use this for mm9
    submit=$(bsub -N -M 40000 -R "rusage[mem=10000]"  -J $JOB -o $LOG "bam2bw.sh $line $prefix  /nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/mm9.chrom.sizes '-f $scalefactors $prefix'")
 fi
 echo $submit

done < $listOfFiles
