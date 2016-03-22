#!/bin/bash

#LSF args
#OPT='-M 20000 -N -R "select[ncores=2]"'

#File should be the files prefix e.g 2lox_Epi
prefix=$1
#Needs to be changed to whatever number the jobs and consequently the files are
N_JOBS=$2
COUNT=$(echo $N_JOBS| grep -oE [[:digit:]] | wc -l)
STRING=$(for(( i=1; i <= $COUNT; i++ )); do echo -n '?'; done)
#OPT="-M 30000 -N"

for suffix in nomapping unpaired_mult unpaired_uniq unpaired_circular unpaired_transloc
do
      echo -e "Running $0 on files: $prefix.$suffix"
      for (( i=0; i<$N_JOBS; i++ ))
      do

            samtools view -bS $prefix.$i.$suffix -o $prefix.$i.$suffix'.bam';
            samtools sort $prefix.$i.$suffix'.bam' $prefix.$i.$suffix'.nsort';
      done

      #If N_JOBS is single digit or triple digit change the number of ?? matching the job number
      samtools merge -f $prefix.$suffix'.bam' $prefix.$STRING.$suffix'.nsort.bam'
      #Delete the intermediate prefixs if all q prefixs sorted and combined to one
      #bsub $OPT -w "done(\"merge_$prefix.$suffix\")" "rm $prefix.??.$suffix*"
done

#For merging the $prefix.uniq and $prefix.mult bam prefixs created above
#samtools merge -f $prefix.bam $prefix.*????.bam



#Picard tools MergeSamFiles path and options
#MergeSam='java -Xmx2g -jar /nfs/research2/bertone/software/picard-tools-1.76/MergeSamFiles.jar'
#
#for suffix in  nomapping unpaired_mult unpaired_uniq unpaired_circular unpaired_transloc
##for suffix in concordant_uniq
# do
#  INPUT=''
#  for ((i=0; i<njobs; i++))
#   do
#    INPUT="$INPUT INPUT=$prefix.$i.$suffix "
#  done
#  bsub $OPT -J "MergeSamPSort.$prefix.$suffix" -o "$prefix.$suffix.MergeSamPSort.log" "$MergeSam $INPUT SORT_ORDER=sorted COMMENT='Product of merging gsnap split output on $njobs. This is sorted based on coordinates.' OUTPUT=$prefix.$suffix.sorted.sam"
#
#done
#


