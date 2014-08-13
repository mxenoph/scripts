#!/bin/bash


#Run this script with a text file containing the files it should run on, a gtf file and a string for the -t and -s htseq-count options if you don't want the defaults e.g. for Epi samples '-s reverse -t intron'
files=$1
gtf=$2
uOPT=$3

htseq='htseq-count_v0.5.3'

#run htseq-count on a list of files that are strand specific
while read infile
do
 echo -e "Running htseq-count on $infile"
 # Options:
 # -m <mode> either uniion, intersection strict, intersection-noempty; default:union
 # -o <samout_name>: write all SAM alignment records into an output SAM file called samout_name
 # -q suppress progress report and warnings
 # -s whether data from strand specific assay; default: yes
 # -t feature type- default is exon
 output=$(echo $infile|sed 's/$/.ss.counts/')

 #For enhancers use intersection-nonempty
 #OPT="-m intersection-nonempty -q -i enhancer_id"

 OPT="-m union -q"


 if [[ $uOPT ]]
  then
  OPT="$OPT $uOPT"

  if [[ $uOPT =~ "-s no" ]]
   then
    output=$(echo $infile|sed 's/$/.nss/')
  fi

  if [[ $uOPT =~ "-s reverse" ]]
   then
    output=$(echo $infile|sed 's/$/.sr/')
   fi


  if [[ "$uOPT" =~ -t\ [a-zA-Z]*\ ?  ]]
   then
    ftype=${BASH_REMATCH[0]}
    ftype=$(echo $ftype| sed 's/-t\s//')
    output="$output.$ftype"
  fi

  #output="$output.counts"

  #for enhancer
  if [[ $gtf =~ "plus" ]]
   then
   output="$output.plus.counts"
   elif [[ $gtf =~ "minus" ]]
   then
   output="$output.minus.counts"
   else
       output=".counts"
  fi
  ###

fi

 OPT="$OPT $infile $gtf"

 submit=$(bsub -n 6 -N -M 70000 -R "rusage[mem=10000]" -J "htseq-count" -oo "$output" $htseq $OPT)
 echo -e "$submit"
done < "$files"
