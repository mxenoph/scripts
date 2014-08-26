#!/usr/bin/perl -w

# From Lars Juhl Jensen (EMBL) <https://www.biostars.org/p/4783/>
use strict;

my $step = 1;
my $window = 100000;

my @sigmoid = (15000, 0.00025, 0, 1);




die "Syntax: score_chip.pl part parts\n" unless scalar @ARGV == 2;
my ($part, $parts) = @ARGV;


my %segments = ();
#open GFF, "gzip -cd dmel-4-*.gff.gz |";
open GFF, "gzip -cd ../dmel-*.gff.gz |";
while (<GFF>) {
  s/\r?\n//;
  next if /^#/;
  my ($chromosome, undef, $type, $left, $right, undef, undef, undef, $name) = split /\t/;
  next unless $type eq "gene" and $left =~ /^[0-9]+$/ and $right =~ /^[0-9]+$/ and $name =~ /(FBgn[0-9]+)/;
  $name = $1;
  my $center = ($left+$right)/2;
  $segments{$chromosome}{$center}{$left}{$right}{$name} = 1;
}
close GFF;


open CHIP, "< regions_merged.tsv";
open OUT, "> score_chip.$part.out";
$_ = <CHIP>;
my $NR = 0;
while (<CHIP>) {

  next unless ($NR++)%$parts == $part;

  s/\r?\n//;
  my ($clone, $chromosome, $start, $stop, undef) = split /\t/;
  next unless $start =~ /^[0-9]+$/ and $stop =~ /^[0-9]+$/;

  my @centers = ();
  foreach my $center (keys %{$segments{$chromosome}}) {
    push @centers, $center if $center >= $start-$window and $center <= $stop+$window;
  }

  my $count = 0;
  my %name_sum = ();
  for (my $pos = $start+int(($stop-$start+1)%$step/2); $pos <= $stop; $pos += $step) {
    $count++;
    my %name_best = ();
    foreach my $center (@centers) {
      foreach my $left (keys %{$segments{$chromosome}{$center}}) {
        foreach my $right (keys %{$segments{$chromosome}{$center}{$left}}) {
          if ($pos < $left) {
            my $score = $sigmoid[2]+($sigmoid[3]-$sigmoid[2])/(1+exp($sigmoid[1]*(($left-$pos)-$sigmoid[0])));
            next unless $score > 1e-6;
            foreach my $name (keys %{$segments{$chromosome}{$center}{$left}{$right}}) {
              $name_best{$name} = $score unless exists $name_best{$name} and $name_best{$name} >= $score;
            }
          } elsif ($pos > $right) {
            my $score = $sigmoid[2]+($sigmoid[3]-$sigmoid[2])/(1+exp($sigmoid[1]*(($pos-$right)-$sigmoid[0])));
            next unless $score > 1e-6;
            foreach my $name (keys %{$segments{$chromosome}{$center}{$left}{$right}}) {
              $name_best{$name} = $score unless exists $name_best{$name} and $name_best{$name} >= $score;
            }
          } else {
            foreach my $name (keys %{$segments{$chromosome}{$center}{$left}{$right}}) {
              $name_best{$name} = 1;
            }
          }
        }
      }
    }
    foreach my $name (keys %name_best) {
      if (exists $name_sum{$name}) {
        $name_sum{$name} += $name_best{$name};
      } else {
        $name_sum{$name} = $name_best{$name};
      }
    }
  }

  foreach my $name (sort {$name_sum{$b} <=> $name_sum{$a}} keys %name_sum) {
    my $score = $name_sum{$name}/$count;
    last if $score <= 0.0005;
    printf OUT "%s\t%s\t%.3f\n", $clone, $name, $score;
  }

}
close CHIP;
close OUT;
