#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long qw (GetOptions);
use File::Basename; # For retrieving basename and extension

my (@annotation);
Getopt::Long::Configure("no_ignore_case", "prefix_pattern=(--|-|\/)");
# Requires at least one file
GetOptions("f=s{1,}" => \@annotation);

# Adds prefix 'chr' and converts 'MT' to 'chrM' for gtf and bed files #{{{
sub gtf {
    my ($gtf, $out) = @_;

    # Open filehandles
    open(IN, "<$gtf") || die "$gtf: No such file\n";
    open(OUT, ">$out") || die "$out: Can not write\n";

    while(<IN>){
        if ($_ =~ /^(\d{1,2})/g) {
            $_ =~ s/^/chr/;
            print OUT "$_";
        }
        elsif ($_=~ /^(MT){1}/) {
            $_ =~ s/^MT/chrM/;
            print OUT $_;
        }
        elsif ($_ =~ /^[X|Y]/) {
            $_ =~ s/^/chr/;
            print OUT $_;
        }
        else {
            print OUT $_;
        }
    }
    
    close(IN);
    close(OUT);
}#}}}

foreach my $file (@annotation){
    my ($basename, $path, $suffix) = fileparse($file, '\.[^\.]*');
    my $out = $path.$basename.'.with-prefix'.$suffix;

    gtf($file, $out);
}

