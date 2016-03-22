#!/usr/bin/env bash
#
#Author: Par Engstrom

SAM_TO_BED_SCRIPT=sam-to-bed.pl
SCALE_SCRIPT=scale-bedgraph.pl

USAGE="Usage: bam2bigwig.sh input.bam outBaseName sizes.txt scale-args"
if [ $# -ne 4 ]; then
    echo $USAGE
    exit 1
fi

inFile=$1
baseFn=$2
sizesFile=$3
scaleArgs=$4

sortedFile=$inFile.sorted.tmp.bed
fwFile=$inFile.fw.bed
revFile=$inFile.rev.bed


# Make sorted bed file
# The '-1' option has to do with the protocol the library was made (similar to the HTSeqCount option -s reverse) | remove it or keep it accordingly
samtools view $inFile | $SAM_TO_BED_SCRIPT -m -s '-1' | sort -k1,1 -k2,2n >$sortedFile

# + strand
awk '{if($6=="+") { print $0 }}' $sortedFile > $fwFile
bedItemOverlapCount -chromSize=$sizesFile -outBounds -bed12 null $fwFile | $SCALE_SCRIPT $scaleArgs |  wigToBigWig stdin $sizesFile ${baseFn}_fwd.bw


# - strand
awk '{if($6=="-") { print $0 }}' $sortedFile > $revFile
bedItemOverlapCount -chromSize=$sizesFile -outBounds -bed12 null $revFile | $SCALE_SCRIPT $scaleArgs | wigToBigWig stdin $sizesFile ${baseFn}_rev.bw

rm $sortedFile
rm $fwFile
rm $revFile
