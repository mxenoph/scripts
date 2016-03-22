#!/usr/bin/perl

use warnings;
use strict;

#e.g. mm10 (name only-no file)
my $genome=shift;
my $in_file=shift;

open(IN, "<$in_file") || die "Can't open $in_file\n";

while (my $file= <IN>){
 chop($file);
 $file=~ s/\n//g;
 next && print "Input file not sorted .bam\n" if ($file !~ /\.bam$/ && $file !~ /sort/g);
 my $out_file= $file;
 $out_file=~ s/\.bam$/\.wig/;
 my $name;
 if ($file=~ /^\//){
  $name=$1 if ($file=~ /.*\/(\w+)\..+/);
 }
 else{ 
  $name=$1 if ($file=~ /(\w+)\..+/);
 }

 print "Computing coverage (separately for each strand) for $name\n";
 my $command="bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" -J \"igv_$name\" -oo \"$out_file\" \"igvtools count -w 1 --strands 'first' $file $out_file $genome\"";
 print "$command\n";
 my $submit=`$command`;
 print "$submit";

 my $trunc_file= $out_file;
 $trunc_file=~ s/\.wig/noTr\.wig/;
 #This will be submitted and run when the above job has finished
 $command="bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" -w 'done(\"igv_$name\")' -J \"grep_$name\" \"grep -v 'track' $out_file > $trunc_file\"";
 $submit=`$command`; 
 print "$submit";
 
 my $neg_file= $out_file;
 $neg_file=~ s/\.wig/.sneg\.wig/;
 $command="bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" -w 'done(\"grep_$name\")' -J \"cut_neg_$name\" \"cut -f 1-2 $trunc_file > $neg_file\"";
 $submit=`$command`; 
 print "$submit";
 #print "$command\n";

 my $pos_file= $out_file;
 $pos_file=~ s/\.wig/.spos\.wig/;
 $command="bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" -w 'done(\"grep_$name\")' -J \"cut_pos_$name\" \"cut -f 1,3 $trunc_file > $pos_file\"";
 $submit=`$command`; 
 print "$submit";

 my $out=$neg_file;
 $out=~ s/\.wig/\.bw/;
 $command="bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" -w 'done(\"cut_neg_$name\")' -J \"bw_neg_$name\" \"wigToBigWig $neg_file > $out\"";
 $submit=`$command`; 
 print "$submit";

 $out=$pos_file;
 $out=~ s/\.wig/\.bw/;
 $command="bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" -w 'done(\"cut_pos_$name\")' -J \"bw_pos_$name\" \"wigToBigWig $pos_file /nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_10/mm10_noHead.genome $out\"";
 $submit=`$command`; 
 print "$submit";

 #print "$command\n";
}
close(IN);

