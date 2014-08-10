#!/usr/bin/env bash

while getopts ":g:o:f:" opt; do
    case $opt in
        g)
            genome=${OPTARG} ;; # genome should be in Caps
        o)
            out=$(sed 's/\$//' <<< ${OPTARG}) ;;
        f)
            files=${OPTARG} ;; # Files to run gsnap on. If paired end separate by ';'
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

if [ ! -d "$out/gsnap" ]
then
    mkdir "$out/gsnap"
fi
cd "$out/gsnap"

IFS=';' read -a arr <<< "$treat"

echo ${arr[@]}[0]
exit 1
#out_dir only for .err files --not in this version
out_dir=$1
out_dir=$(echo $out_dir|sed 's/\/$//')
files=$2
#should be in Caps
genome=$3

if [ ! -d "$out_dir" ]
then 
 mkdir "$out_dir"
 echo -e "Creating $out_dir"
fi 

#If using a gsnap version not in the /common/software change it and give full path if not in your $PATH
gsnap="/nfs/research2/bertone/user/myrto/software/external/gmap-current/bin/gsnap"

# Options:
# -m: number of mismatches (for 75bp set to 4 or 5); 
# -i: indel penalty (default=2)...<2 can lead to FP at the end of the read
# -N: look for novel splicing
# -w: for novel splicing event - so it doesn't consider a splice junction that is too big (runs from one end of the chromosome to the other)
# -n: number of paths to print
# -Q: print nothing if >n paths
# -O: relevant if more than one worker thread (-t) -> prints output in the same order as input

N_THREADS=10
#kmer size to use in genome database--if not specified highest available in the db is used
#K_MER=15
#For running with -q option
N_JOBS=20
MEM=25000
mismatch=7
PROTOCOL='sanger'
SPLICE_SITES='/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/mm10_ensembl70_splicesites.iit'

if [ "$genome" == 'MM9' ]
 then
  SPLICE_SITES='/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/MM9.maps/mm9_ensembl67_splicesites.iit'
fi

#run gsnap on a list of files 
while read infile
do
 echo -e "Running gsnap on $infile"
 #Since I'm usually giving full path for the file in the input list-file this will not work if in my previous directories I have an '_', so be CAREFULL!
 name=$(echo "$infile"|awk -F'_' '{print $2}'| sed s/-/_/)
 #For PE reads the fastq.gz listed in the input file should be mates 1 
 mate_1=$(echo "$infile")
 mate_2=$(echo "$infile"|sed s/1_sequence/2_sequence/)

 jobName="gsnap_pe.$name[1-$N_JOBS]"
 
 #for i in {0..$(($N_JOBS-1))}; do
 for ((i=0; i<$N_JOBS; i++)); do
  #Make string with gsnap options
  #Jobs and threads
  OPT="-q $i/$N_JOBS -t $N_THREADS"
  OPT="$OPT -D /nfs/research2/bertone/common/genome/$genome -d $genome"
  #Option -B is batch mode to make it run faster| (-i) --indel-penalty by default is 2, what is the benefit of setting it to mismatches? --trim-mismatch-score=0 turns off trimming at ends and will give FP at the ends of reads 
  OPT="$OPT -B 5 -m $mismatch -i 1"
  #Find novel splicing and limit the size of splice junctions, define distant-splice penalty (-E) for when the intron length exits the value of -w
  OPT="$OPT -N 1 -w 100000 -E 100 -s $SPLICE_SITES"
  #Number of paths to print and 'don't print failed alignments', name output based on the q option split up
  OPT="$OPT -n 10 -Q -O --nofails -A sam --split-output=$out_dir/$name.$i --gunzip --quality-protocol=$PROTOCOL"
  
  #echo -e
  #echo -e bsub -q research-rh6 -n $N_THREADS -N -M $MEM -R "rusage[mem=$MEM]" -oo "$out_dir/gsnap_$name.$i.log" -J "$name.$i" $gsnap $OPT $mate_1 $mate_2
  submit=$(bsub -q research-rh6 -n $N_THREADS -N -M $MEM -R "rusage[mem=$MEM]" -oo "$out_dir/gsnap_$name.$i.log" -J "$name.$i" $gsnap $OPT $mate_1 $mate_2)
  echo -e "$submit"
  done
done < "$files"
