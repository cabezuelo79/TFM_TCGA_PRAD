#!usr/bin/perl
use Data::Dumper;
open (IN, "prueba.tsv");

my %hash;
my %patient;
my %geneid;

while(my $linea=<IN>){
	chomp($linea);
	if($linea =~/(^DO\w+)\t.+\t.+\t.+\t.+\t.+\t(.+)\t.+\t(.+)\tGRCh37/) {
		my $id=$1;
		my $gene=$2;
		my $raw=$3;
		$geneid{$gene}=$raw;
		$patient{$id}=$raw;
		$hash{$gene}{$id}=$raw;
	}
}
close (IN);

foreach $id (sort keys %patient) {
               print $id."\t";
}
print "\n";

foreach $gene (keys %geneid) {
                print $gene."\t"; 
                       foreach $id (sort keys %patient) {
                               foreach ($hash{$gene}{$id}) {
                                          print $raw;
                                          print "\n";
                                   }
                  }
}