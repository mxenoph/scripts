#!/bin/env perl
# from mara
use warnings;
use strict;

sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: perl $cmd chrom_sizes window_size sorted.bam

Modes:

EndOfUsage
    exit(1);
}

@ARGV == 3 or usage();

my $chr_file = $ARGV[0];
my $window_size = $ARGV[1];
my $sam_file = $ARGV[2];

my %chr = ();
my %bins = ();

open CHR, "$chr_file" or die "Could not open $chr_file\n";
while (my $line = <CHR>) {
	my @values = split("\t",$line);
	$chr{$values[0]} = $values[1];
}
close CHR;

my $current_chr = "";
my $current_end = $window_size;
my $current_length = $window_size;
my $count = 0;

open SAM, "samtools view -F 0x4 $sam_file | " or die "Could not open $sam_file\n";
#open SAM, $sam_file or die "Could not open $sam_file\n";
while (my $sam_line = <SAM>) {
	my @sam_values = split("\t",$sam_line);
#	print "$sam_values[2]\t$sam_values[3]\n";
	if ($sam_values[2] ne $current_chr) {
		if ($current_chr ne "") { $bins{"$current_chr:$current_end"} = $count; } # / $current_length; }
		while ($current_chr ne "" && $current_end < $chr{$current_chr}) {
			if ($current_end + $window_size < $chr{$current_chr}) {
        		$current_end = $current_end + $window_size;
        		$current_length = $window_size;
    		} else {
        		$current_length = $chr{$current_chr} - $current_end;
        		$current_end = $chr{$current_chr};
    		}
			$bins{"$current_chr:$current_end"} = 0;
		}
		$current_chr = $sam_values[2];
		$current_end = $window_size;
		$current_length = $window_size;
		$count = 0;
	}
	if ($sam_values[3] > $current_end) {
		while ($current_end < $sam_values[3]) {
			$bins{"$current_chr:$current_end"} = $count; # / $current_length;
			if ($current_end + $window_size < $chr{$current_chr}) {
				$current_end = $current_end + $window_size;
				$current_length = $window_size;
			} else {
				$current_length = $chr{$current_chr} - $current_end;
				$current_end = $chr{$current_chr};
			}
			$count = 0;
		}
	}
	$count++;
}
close SAM;

while ($current_end < $chr{$current_chr}) {
	$bins{"$current_chr:$current_end"} = 0; # / $current_length;
	if ($current_end + $window_size < $chr{$current_chr}) {
		$current_end = $current_end + $window_size;
		$current_length = $window_size;
	} else {
		$current_length = $chr{$current_chr} - $current_end;
		$current_end = $chr{$current_chr};
	}
}

foreach my $key (keys %bins) {
	print "$key\t$bins{$key}\n";
}
