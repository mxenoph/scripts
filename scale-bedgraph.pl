#!/usr/bin/perl
#
#Author:Par Engstrom
#
# not using /bin/env perl because it complains about IO.sobeing compiled by a different version
use warnings;
use strict;
use Getopt::Std;

sub usage()
{
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage:
$cmd <scale_factor>
$cmd -f <scale_factor_file> <key>

EndOfUsage
    exit(1);
}


my %args;
getopts("f:", \%args);

usage() unless(@ARGV == 1);

# Determine scale factor
my $scale_factor;
if($args{f}) {
    open IN, $args{f} or die "ERROR: failed to open scale factor file $args{f} for input.\n";
    while(<IN>) {
	chomp;
	my ($key, $factor) = split "\t", $_;
	if($key eq $ARGV[0]) {
	    $scale_factor = $factor;
	    last;
	}
    }
    close IN;
    defined($scale_factor) or die "ERROR: key $ARGV[0] not found in scale factor file $args{f}.\n"; 
}
else {
    $scale_factor = $ARGV[0];
}

# Process bedGraph data
while(<STDIN>) {
    chomp;
    my @f = split /\t/, $_;
    $f[3] = $f[3] / $scale_factor;
    print join("\t", @f), "\n";
}
