#!/user/bin/env perl -w

#------------------------------------------------------------------------------
#
$usage = qq(
Usage:  $0 <GENES-FILE> <BEDGRAPH-FILE> <OUT-PREFIX>\n);

die($usage) if (@ARGV < 3);
#-----------------------------------------------------------------------------
use List::Util qw[min max];
use File::Basename;

my $GENES = $ARGV[0];
my $BEDGRAPH = $ARGV[1];

$genes = read_genes_into_hash($GENES);
$bedgraph = read_bedgraph_into_hash($BEDGRAPH);

my $OUT = $ARGV[2]."_gene_signal.txt";
open(OUT, ">$OUT") or die "Can't open $OUT\n";

# compute scores for genes
foreach my $gene (keys %$genes) {
  $genes->{$gene}->{'prom_max'}=compute_max_bin($genes->{$gene}->{'chr'},([$genes->{$gene}->{'prom_start'},$genes->{$gene}->{'prom_end'}]));
  $genes->{$gene}->{'body_max'}=compute_max_bin($genes->{$gene}->{'chr'},([$genes->{$gene}->{'body_start'},$genes->{$gene}->{'body_end'}]));
}


# output results
print OUT "gene.id\tgene.name\tprom.max\tbody.max\n";
foreach my $gene (keys %$genes) {
  print OUT "$gene\t$genes->{$gene}->{'name'}\t$genes->{$gene}->{'prom_max'}\t$genes->{$gene}->{'body_max'}\n";
}

close(OUT);



sub compute_max_bin {
  my ($chr,@coords) = @_;

  my $max=0;
    
  foreach (@coords){
    my ($start,$end)=@{$_};
    my $pos=$start-200;
   # print "$chr\t$pos\t$start\t$end\n";
    while ($pos < $end) {
      if (exists $bedgraph->{$chr}->{$pos} && $bedgraph->{$chr}->{$pos}->{'end'} > $start) {
	$max=max($max,$bedgraph->{$chr}->{$pos}->{'count'});
      }
      $pos++;
    }
  }
  return($max);
}

sub read_bedgraph_into_hash {
  my ($filename) = @_;
  my $hash_ref ={};
  open(FILE, $filename) or die("Could not open file!");
  while(my $line = <FILE>) {
    chomp($line);
    my @data = split('\t', $line);
    # adapt for if chrs are in EnsEMBL format
    my $chr= $data[0]=~/chr/ ? $data[0] : "chr".$data[0];
    $chr=~ s/\s//g;
    if ($data[3] > 5) {
      $hash_ref->{$chr}->{$data[1]}->{'end'}=$data[2];
      $hash_ref->{$chr}->{$data[1]}->{'count'}=$data[3];
    }
  }
  close(FILE);
  return ($hash_ref);
}

sub read_genes_into_hash {
  my ($filename) = @_;
  my $hash_ref ={};
  open(FILE, $filename) or die("Could not open file!");
  while(my $line = <FILE>) {
    chomp($line);
    my @data = split('\t', $line);
    my $gene_id= $data[0];
    if ($gene_id=~/ENSMUSG/) {
      # adapt for if chrs are in EnsEMBL format
      my $chr= $data[3]=~/chr/ ? $data[3] : "chr".$data[3];
      $chr=~ s/\s//g;
      $hash_ref->{$gene_id}->{'chr'}=$chr;
      $hash_ref->{$gene_id}->{'start'}=$data[4];
      $hash_ref->{$gene_id}->{'end'}=$data[5];

      $strand=$data[8];

      if ($strand eq "1") {
	$hash_ref->{$gene_id}->{'tss'}=$data[6];
	#$hash_ref->{$gene_id}->{'prom_start'}=$data[6]-1000;
	$hash_ref->{$gene_id}->{'prom_start'}=$data[6]-2000;
	$hash_ref->{$gene_id}->{'prom_end'}=$data[6]+500;
	$hash_ref->{$gene_id}->{'body_start'}=$data[6]+501;
	$hash_ref->{$gene_id}->{'body_end'}=$data[5];
      }
      else {
	$hash_ref->{$gene_id}->{'tss'}=$data[7];
	$hash_ref->{$gene_id}->{'prom_start'}=$data[7]-500;
	#$hash_ref->{$gene_id}->{'prom_end'}=$data[7]+1000;
	$hash_ref->{$gene_id}->{'prom_end'}=$data[7]+2000;
	$hash_ref->{$gene_id}->{'body_start'}=$data[4];
	$hash_ref->{$gene_id}->{'body_end'}=$data[7]-501;
      }

      $hash_ref->{$gene_id}->{'name'}=$data[9];
      $hash_ref->{$gene_id}->{'type'}=$data[18];
      
    }
  }
  close(FILE);
  return ($hash_ref);
}








