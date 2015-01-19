#!/bin/bash

###########################################
###  Always run this script under bsub! ###
###########################################


#File should be the files prefix e.g 2lox_Epi
file=$1
#Needs to be changed to whatever number the jobs and consequently the files are
N_JOBS=20
OPT="-M 30000 -N"

for suffix in halfmapping_mult halfmapping_transloc halfmapping_uniq concordant_mult concordant_transloc concordant_uniq 
 do
 echo -e "Running $0 on files: $file.$suffix"
  for ((i=0; i<$N_JOBS; i++ ))
   do
   bsub $OPT -J "bam_$file.$suffix.$i" "samtools view -bS $file.$i.$suffix -o $file.$i.$suffix.bam"
   sleep 15
   bsub $OPT -w "done(\"bam_$file.$suffix.$i\")" -J "sort_$file.$suffix.$i" "samtools sort $file.$i.$suffix $file.$i.$suffix'_sorted'"  
   sleep 15
    #If N_JOBS is single digit or triple digit change the number of ?? matching the job number 
    bsub $OPT -w "done(\"sort_$file.$suffix.$i\")" -J "merge_$file.$suffix" "samtools merge -f $file.$suffix'_sorted.bam' $file.??.$suffix_sorted.bam"
    sleep 15
    #Delete the intermediate files if all q files sorted and combined to one
    bsub $OPT -w "done(\"merge_$file.$suffix\")" "rm $file.??.$suffix*"
   done
done

sleep 15
#Merge the uniq mapping to one file-- but why not using paired uniq also?
bsub $OPT -w "done(\"merge_$file.concordant_uniq\") && done(\"merge_$file.halfmapping_uniq\")" -J "merge_$file.uniq" "samtools merge -f $file.uniq_sorted.bam $file.concordant_uniq_sorted.bam $file.halfmapping_uniq_sorted.bam"
   
sleep 15
#Merge the mult mapping to one file-- but why not using paired mult also?
bsub $OPT -w "done(\"merge_$file.concordant_mult\") && done(\"merge_$file.halfmapping_mult\")" -J "merge_$file.mult" "samtools merge -f $file.mult_sorted.bam $file.concordant_mult_sorted.bam $file.halfmapping_mult_sorted.bam"
 
#For merging the $file.uniq and $file.mult bam files created above
bsub $OPT -w "done(\"merge_file.uniq\") && done(\"merge_$file.mult\")" "samtools merge -f $file.bam $file.????.bam" 
