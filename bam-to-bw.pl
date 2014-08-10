#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw (GetOptions);

my (@files, $genome, $out);
Getopt::Long::Configure("no_ignore_case", "prefix_pattern=(--|-|\/)");
# Requires at least one file
GetOptions("f=s{1,}" => \@files,
            "g=s{1,1}" => \$genome,
            "o=s{1,1}" => \$out);

my $cluster_usage= "bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" ";
$out=~ s/$/\// if ($out !~ /\/$/);

foreach my $f (@files){
    (my $str = $f) =~ s/\.bam//;
    $str=~ s/.+\///;

    my $genomecov= " \"bedtools genomecov -bg -ibam $f -g $genome -split > $out$str.bedGraph\"";
    my $cmd= $cluster_usage."-J \"bed_$str\" -e \"$str.genomcov.log\"".$genomecov;
#    print $cmd, "\t";
    `$cmd`;

    my $bw= " \"bedGraphToBigWig $out$str.bedGraph $genome $out$str.bw\"";
    $cmd= $cluster_usage."-w 'done(\"bed_$str\")' -J \"bw_$str\"".$bw;
    `$cmd`;
}

