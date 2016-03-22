#!/bin/bash


#e.g. $0 2lox_Epi rsorted.sam queryname
prefix=$1
suffix=$2
#queryname|coordinates
sorted=$3

#LSF args
OPT='-n 6 -M 20000 -N'

#Picard tools MergeSamFiles path and options
MergeSam='java -Xmx2g -jar /nfs/research2/bertone/software/picard-tools-1.76/MergeSamFiles.jar'
INPUT_U="INPUT=$prefix.concordant_uniq.$suffix INPUT=$prefix.halfmapping_uniq.$suffix"
INPUT_M="INPUT=$prefix.concordant_mult.$suffix INPUT=$prefix.halfmapping_uniq.$suffix"
nameSuf='sorted.sam'
jobSuf='PSortHC'

if(( "$sorted" == "queryname"))
 then
  nameSuf='rsorted.sam'
  jobSuf='RSortHC'
fi

bsub $OPT -J "MergeSam$jobSuf.$prefix.uniq" -o "$prefix.uniq.MergeSam$jobSuf.log" "$MergeSam $INPUT_U SORT_ORDER=$sorted COMMENT='Product of merging gsnap halfmapping_uniq and concordant_uniq. This is sorted based on $sorted.' OUTPUT=$prefix.uniq.$nameSuf"
  
bsub $OPT -J "MergeSam$jobSuf.$prefix.mult" -o "$prefix.mult.MergeSam$jobSuf.log" "$MergeSam $INPUT_M SORT_ORDER=$sorted COMMENT='Product of merging gsnap halfmapping_mult and concordant_mult. This is sorted based on $sorted.' OUTPUT=$prefix.mult.$nameSuf"

#Add -w option of LSF 
bsub $OPT -w "done(\"MergeSam$jobSuf.$prefix.uniq\") && done(\"MergeSam$jobSuf.$prefix.mult\")" -J "MergeSam$jobSuf.$prefix.all" -o "$prefix.MergeSam$jobSuf" "$MergeSam $INPUT_M SORT_ORDER=$sorted COMMENT='Product of merging *.uniq.* and *.mult.*. This is sorted based on $sorted.' OUTPUT=$prefix.$nameSuf"
