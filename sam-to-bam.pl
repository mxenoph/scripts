#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw (GetOptions);

my (@files, $out, $pe);
Getopt::Long::Configure("no_ignore_case", "prefix_pattern=(--|-|\/)");
# Requires at least one file
GetOptions("f=s{1,}" => \@files,
           "s:i" => \$pe,
           "o=s{1,1}" => \$out);

# If output directory not specified set to PWD
$out= $ENV{'PWD'} if (! defined $out);
# If -s option not passed set pe to 0 -- Not paired end data -> numeric sort
$pe= 0 if (! defined $pe);

my $cluster_usage= "bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" ";
$out=~ s/$/\// if ($out !~ /\/$/);

foreach my $f (@files){
(my $str = $f) =~ s/\.sam//;
$str=~ s/.+\///;

# If argument not provided I'm assuming data are single end -> numeric sort
my $bamconv= " \"samtools view -bS $f -o $out$str.bam\"";
my $cmd= $cluster_usage."-J \"sam2bam_$str\" -e \"$str.bamconv.log\"".$bamconv;
#print $cmd, "\n";
`$cmd`;

my $bamsort= " \"samtools sort $out$str.bam $out$str.sort\"";
$bamsort= " \"samtools sort $out$str.bam $out$str.sort\"" if ( $pe  != 0 );
$cmd= $cluster_usage."-w 'done(\"sam2bam_$str\")' -J \"bamsort_$str\" -e \"$str.bamsort.log\"".$bamsort;
`$cmd`;

my $bamind= " \"rm $out$str.bam; rename '.sort' '' $out$str.sort.bam; samtools index $out$str.bam\"";
$cmd= $cluster_usage."-w 'done(\"bamsort_$str\")' -J \"bamsort_$str\" -e \"$str.bamind.log\"".$bamind;
#print $cmd, "\n";
`$cmd`;
}

