#!/usr/local/bin/perl

use warnings;
use strict;

use DBI;

sub usage {
    my ($cmd) = $0 =~ /([^\\\/]+)$/;
    print STDERR <<EndOfUsage;
    
$cmd

Usage: $cmd <input.genePred> <db_name>

EndOfUsage
    exit;
}

usage() unless(@ARGV == 2);

my $genePred = $ARGV[0];
my $db_name = $ARGV[1];

my $db = DBI->connect("dbi:SQLite:$db_name","","",{RaiseError => 1, AutoCommit => 1});
$db->do("CREATE TABLE GeneTable (chrom TEXT, name TEXT PRIMARY KEY, strand TEXT, txStart INTEGER, txEnd INTEGER, cdsStart INTEGER, cdsEnd INTEGER, exonCount INTEGER, exonStarts TEXT, exonEnds TEXT, name2 TEXT)");

open IN, "$genePred" or die "Could not open file $genePred\n";
while (my $line = <IN>) {
	chomp $line;
	my @values = split("\t",$line);
	$db->do("INSERT INTO GeneTable VALUES (\'$values[1]\',\'$values[0]\',\'$values[2]\',\'$values[3]\',\'$values[4]\',\'$values[5]\',\'$values[6]\',\'$values[7]\',\'$values[8]\',\'$values[9]\',\'$values[11]\')");
}

