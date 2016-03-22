#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw (GetOptions);
use File::Basename;
use File::Path qw (make_path);

my (@bams, $genome, $output_path);
Getopt::Long::Configure("no_ignore_case", "prefix_pattern=(--|-|\/)");
GetOptions( "b=s{1,}" => \@bams, # Requires at least one file
            "g=s{1,1}" => \$genome,
            "o=s{1,1}" => \$output_path );

my $cluster_usage= "bsub -n 2 -M 20000 -N -R \"rusage[mem=5000]\" -R \"select[ncores>2]\" ";

# Format output_path
my @tmp = fileparse($output_path);
$output_path =~ $tmp[1];
# Create output directory if does not exist
make_path( $output_path, {verbose => 1, mode => 0711} );

foreach my $bam ( @bams ){
    my ($basename, $path, $suffix) = fileparse($bam, '\.[^\.]*');

    my $genomecov= " \"bedtools genomecov -bg -ibam $bam -g $genome -split > $output_path$basename.bedGraph\"";
    my $cmd= $cluster_usage."-J \"bed_$basename\" -e \"$basename.genomcov.log\"".$genomecov;
#    print $cmd, "\t";
    `$cmd`;

    my $bw= " \"bedGraphToBigWig $output_path$basename.bedGraph $genome $output_path$basename.bw\"";
    $cmd= $cluster_usage."-w 'done(\"bed_$basename\")' -J \"bw_$basename\"".$bw;
    `$cmd`;
}

