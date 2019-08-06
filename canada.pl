#!usr/bin/perl
use Data::Dumper;
open (IN, "prueba.tsv");

use strict;

my %hash;
my %patient;
my %geneid;

while(<IN>){
	chomp;
	if($_ !~/^icgc_donor_id/){
		#Esto seria más facil con un split de tabuladores porque el fichero es tsv
		my @data=split(/\t/);
		$geneid{$data[7]}++;
		$patient{$data[0]}++;
		$hash{$data[7]}{$data[0]}=$data[9];
	}
		
}
close (IN);

####### Lo más fácil #####################################
#foreach $id (sort {$a cmp $b} keys %patient) {
#	print $id."\t";
#}
#print "\n";
#
#foreach my $genename (sort {$a cmp $b} keys %geneid) {
#	print $genename."\t"; 
#	foreach my $paciente_id (sort {$a cmp $b} keys %patient) {
#		print $hash{$genename}{$paciente_id} ."\t";
#	}
#	print "\n";
#}
#########################################################


#Lo mas correcto (lo anterior incluye un tab al final)###

#Cabeceras unidas por un tab. Incluyo un tab al principio para q se vean bien los genes (creo que a R esto no le gusta asi
# que despues de comprobar que esta bien, ese tab habria que quitarlo)

print "GeneNames\t". join("\t",sort {$a cmp $b} keys %patient ) ."\n";
#Recorro los genes
foreach my $genename (sort {$a cmp $b} keys %geneid) {
	#Lo imprimo
	print $genename."\t";
	my @values;
	foreach my $paciente_id (sort {$a cmp $b} keys %patient) {
		#almaceno los valores en un array (que es ordenado) y los imprimo debajo con un join
		push(@values,$hash{$genename}{$paciente_id});
	}
	print STDERR scalar(@values) ."\n";
	print join("\t",@values) ."\n";
}

#########################################################
