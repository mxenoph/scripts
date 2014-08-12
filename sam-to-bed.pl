#!/bin/env perl
#
#Author:Par Engstrom
#
#
use warnings;
use strict;
use Getopt::Std;
use Genoman::BlockSet;
use Genoman::GenomeFeature::Generic;
use Genoman::File::BED;

my $STRAND = 1;

main();
exit(0);


sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: $cmd [options] 

The script expects SAM on standard input and writes BED to standard output.

Options:

-m      Merge alignments of the same fragment into one entry per chromosome
        and strand. The input must be sorted by fragment ID for this to work.

-s int  Set to 1 or -1 to indicate the strand on which to count a first mate
        aligned to the forward strand. Default: 1.

EndOfUsage
    exit(1);
}


sub main {

    # Parse args
    my %args;
    getopts("fms:", \%args);
    if(defined $args{s}) {
	usage() unless($args{s} eq '1' or $args{s} eq '-1');
	$STRAND = $args{s};
    }
    my $merge = $args{m};
    usage if(@ARGV);

    # Create output object
    my $out = Genoman::File::BED->new(-fh => \*STDOUT);

    # Process alignments
    if($merge) {
	my @sam_array;
	my $cur_id = "";
	while(my $line = <STDIN>) {
	    my @sam = split "\t", $line, 7; # Split line. We only care about the first 6 fields.
	    next if($sam[1] & 0x4); # Skip non-alignments
	    if($sam[0] ne $cur_id) {
		merge_and_print_frag($out, \@sam_array) if($cur_id ne "");
		$cur_id = $sam[0];
		@sam_array = ( \@sam );
	    }
	    else {
		push @sam_array, \@sam;
	    }
	}
	merge_and_print_frag($out, \@sam_array) if($cur_id ne "");
    }
    else {
	while(my $line = <STDIN>) {
	    my @sam = split "\t", $line, 7; # Split line. We only care about the first 6 fields.
	    next if($sam[1] & 0x4); # Skip non-alignments
	    my ($chr, $strand, $blocks) = parse_sam(\@sam);
	    print_bed($out, $sam[0], $chr, $strand, $blocks);
	}
    }
}


sub merge_and_print_frag {
    my ($out, $sam_array) = @_;

    my (%fwd_aln, %rev_aln);

    foreach my $sam (@$sam_array) {
	my ($chr, $strand, $blocks) = parse_sam($sam);
	if($strand == 1) {
	    push(@{$fwd_aln{$chr}}, $blocks);
	}
	else {
	    push(@{$rev_aln{$chr}}, $blocks);
	}
    }

    while(my ($chr, $aln_array) = each %fwd_aln) {
	my $blocks =  shift @$aln_array;
	$blocks = $blocks->union(@$aln_array);
	print_bed($out, $sam_array->[0][0], $chr, 1, $blocks);
    }
    while(my ($chr, $aln_array) = each %rev_aln) {
	my $blocks =  shift @$aln_array;
	$blocks = $blocks->union(@$aln_array);
	print_bed($out, $sam_array->[0][0], $chr, -1, $blocks);
    }

}


sub parse_sam {
    my $sam = shift;
    my $strand = $STRAND;
    $strand = -$strand if($sam->[1] & 0x10);
    $strand = -$strand if(($sam->[1] & 0x81) == 0x81);
    return ($sam->[2], $strand, cigar_to_blocks($sam->[3], $sam->[5]));
}


sub cigar_to_blocks {
    my ($pos, $cigar) = @_;
    my $blocks;
    my @cig_fields = split(/([A-Z])/, $cigar);
    while(@cig_fields) {
	my $l = shift @cig_fields;
	my $op = shift @cig_fields;
	if($op eq 'M') {  # here we could also catch X and =, but I don't think they are used
	    if($l != 0) {
		push @$blocks, $pos;
		$pos += $l;
		push @$blocks, $pos-1;
	    }
	}
	elsif($op eq 'N' or $op eq 'D') {
	    $pos += $l;
	}
	elsif($op ne 'I' and $op ne 'S') {
	    die "I don't know what to do with CIGAR operation $op";
	}
    }
    return Genoman::BlockSet->new_fast($blocks);
}



sub print_bed {
     my ($out, $id, $chr, $strand, $blocks) = @_;
     my $ft = Genoman::GenomeFeature::Generic->new(
	 -name => $id, -ft_type => 'alignment',
	 -seq_id => $chr, -blocks => $blocks, -strand => $strand);
     $out->store_feature($ft);
}

# sub print_bed {
#     my ($id, $chr, $strand, $blocks) = @_;

#     $strand = strand_numeric_to_sign;

#     my (@starts, @sizes);
#     for(my $i = 0; $i < @$blocks; $i += 2) {
#         push @starts, $blocks->[$i] - $blocks[0];
#         push @sizes, $blocks->[$i+1] - $blocks->[$i] + 1;
#     }
#     push @row, (scalar(@sizes), join(',',@sizes,''), join(',',@starts,''));

#     print join("\t", $chr, $blocks->[0]-1, $blocks->[-1], $id, 0, $strand, );
# }


