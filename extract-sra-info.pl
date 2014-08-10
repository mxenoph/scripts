#!/usr/bin/perl

use warnings;
use strict;

open(SRA, "< $ARGV[0]") || die "Could not open $ARGV[0]";

my (%info, %gse);

while (my $line = <SRA>){
    my @elements= split(',', $line);
    next if ($.==1);

    if($elements[11] !~ /\w/){
        if(!keys %gse){
            die "GSE* file not provided" if(!exists $ARGV[1]);
            %gse=readGSE($ARGV[1]);
            $info{$elements[0]}=$gse{$elements[28]}."\t".$elements[12]
        }
        else{
            $info{$elements[0]}=$gse{$elements[28]}."\t".$elements[12]
        }
    }
    else{
        $info{$elements[0]}=$elements[11]."\t".$elements[12]
    }
}
close(SRA);

my $out=$ARGV[0].".info";
open(OUT, "> $out") || die "Could not open $out";
foreach(keys %info){
    print OUT $_, "\t", $info{$_}, "\n";
}
close(OUT);

sub readGSE{
    my $file= shift;
    my %tmp;

    $/='^SAMPLE';
    open(GSE, "< $file") || die "Could not open $file";
    while (my $chunk= <GSE>){
        next if ($.==1);
        my ($gsm, $info);
        if ($chunk =~ /(GSM\d+).+\s+.+Sample_title\s*=\s*(\S+)/) {
            $gsm=$1; 
            $info=$2;
            $info=~ s/-\/-/_KO_/;
            $info=~ s/__/_/;
        }
        else{
            return(%tmp) if ($chunk !~ /^\s+?$/);
        }
        $tmp{$gsm}=$info;
    }
    $/="\n";
    close(GSE);
    return(%tmp)
}
