#!/usr/bin/env bash


#out_dir only for .err file
out_dir=$1
files=$2
genome=$3

if [ ! -d "$out_dir" ]
then
 mkdir "$out_dir"
 echo -e "Creating $out_dir"
fi

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
 name=$(echo "$infile"|awk -F'_' '{print $2 "_" $3}')
 jobName="gsnap.$name[1-$N_JOBS]"

for ((i=0; i<$N_JOBS; i++)); do
   OPT="-q $i/$N_JOBS -t $N_THREADS"
   OPT="$OPT -D /nfs/research2/bertone/common/genome/$genome -d $genome"
   #Option -B is batch mode to make it run faster| (-i) --indel-penalty by default is 2, what is the benefit of setting it to mismatches? --trim-mismatch-score=0 turns off trimming at ends and will give FP at the ends of reads
   OPT="$OPT -B 5 -m $mismatch -i 1"
   #Find novel splicing and limit the size of splice junctions, define distant-splice penalty (-E) for when the intron length exits the value of -w
   OPT="$OPT -N 1 -w 100000 -E 100 -s $SPLICE_SITES"
   #Number of paths to print and 'don't print failed alignments', name output based on the q option split up
   OPT="$OPT -n 10 -Q -O --nofails -A sam --split-output=$out_dir/$name.$i --gunzip --quality-protocol=$PROTOCOL"

   submit=$(bsub -q research-rh6 -n $N_THREADS -N -M $MEM -R "rusage[mem=$MEM]" -oo "$out_dir/gsnap_$name.$i.log" -J "$name.$i" $gsnap $OPT $infile)

   #NOTE: the splicesites.iit is hard-coded in the command..be sure it's the one you want before you run it
  #submit=$(bsub -n 6 -N -M 20000 -R "rusage[mem=10000]" -J "gsnap" -oo "$out_dir/gsnap_$name.err" gsnap -d "$genome" -t 4 -m 5 -i 1 -N 1 -w 100000 -E 100 -n 10 -Q -q $i/10 -O -A sam --split-output="$name" --gunzip --quality-protocol=sanger -s /nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/MM10.maps/mm10_ensembl70_splicesites.iit "$infile")
  #submit=$(bsub -n 6 -N -M 20000 -R "rusage[mem=10000]" -J "gsnap" -oo "$out_dir/gsnap_$name.err" gsnap -d "$genome" -t 4 -m 5 -i 1 -N 1 -w 100000 -E 100 -n 10 -Q -q $i/10 -O -A sam --split-output="$name" --gunzip --quality-protocol=sanger -s /nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/MM9.maps/mm9_ensembl67_splicesites.iit "$infile")
  echo -e "$submit"
  done
done < "$files"
