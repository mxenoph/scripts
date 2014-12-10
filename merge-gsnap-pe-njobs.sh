#!/bin/bash

prefix=$1
njobs=$2

#LSF args
OPT='-n 6 -M 20000 -N'

#Picard tools MergeSamFiles path and options
MergeSam='java -Xmx2g -jar /nfs/research2/bertone/software/picard-tools-1.76/MergeSamFiles.jar'

for suffix in concordant_circular concordant_mult concordant_transloc concordant_uniq halfmapping_circular halfmapping_mult halfmapping_transloc halfmapping_uniq nomapping paired_mult paired_uniq_circular paired_uniq_inv paired_uniq_long paired_uniq_scr unpaired_circular unpaired_mult unpaired_transloc unpaired_uniq 
#for suffix in concordant_uniq  
 do
  INPUT=''
  for ((i=0; i<njobs; i++))
   do
    INPUT="$INPUT INPUT=$prefix.$i.$suffix "
  done
  bsub $OPT -J "MergeSamPSort.$prefix.$suffix" -o "$prefix.$suffix.MergeSamPSort.log" "$MergeSam $INPUT SORT_ORDER=coordinate COMMENT='Product of merging gsnap split output on $njobs. This is sorted based on coordinates.' OUTPUT=$prefix.$suffix.sorted.sam"
 
  #Sorting by read name is only useful for htseq hence for the halfmapping and concrdant uniq files
  if [[ "$suffix" == "concordant_uniq" || ( "$suffix" == "concordant_mult" ) || ( "$suffix" == "halfmapping_uniq" ) || ( "$suffix" == "halfmapping_mult" ) ]]
   then
  bsub $OPT -J "MergeSamRSort.$prefix.$suffix" -o "$prefix.$suffix.MergeSamRSort.log" "$MergeSam $INPUT SORT_ORDER=queryname COMMENT='Product of merging gsnap split output on $njobs. This is sorted based on query name.' OUTPUT=$prefix.$suffix.rsorted.sam"
  fi
done
