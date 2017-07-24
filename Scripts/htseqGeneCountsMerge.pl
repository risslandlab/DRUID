#!/usr/bin/env perl


use v5.20.0;

use strict;
use warnings;

use autodie;

opendir my $pwd, '.'; 
my @tsv_files = grep { /tsv$/ } readdir $pwd;
my %indices;
my @types;
my $table_ref = {};
while ( my $file = shift @tsv_files ) {
	say $file;
	my ( $index ) = $file =~ /^(.*)_[^_]+\.tsv$/;
	my ( $type ) = $file =~ /_([^_]+)\.tsv$/;
	$indices{$index} = 1;
	push @types, $type;
	open my $file, '<', $file;
	while ( <$file> ) {
		chomp;
		my @values = split /\t/;
		next if $values[0] eq "";
		next if $values[0] =~ /^__/;
		my $gene = $values[0];
		my $FPKM = $values[1];
		$table_ref->{$type}{$gene}{$index} = $FPKM;
	}
	close $file;
}

my @indices = sort { $a cmp $b } keys %indices;
for my $type (@types){
	open my $out, '>', "combined_${type}.txt";
	say $out join( "\t", ("Gene", "Organism", @indices ) );
	my $type_table_ref = $table_ref->{$type};
	for my $gene ( sort( keys( %$type_table_ref ) ) ) {
		my $geneName = substr($gene, 0, -1);
		my $organism = substr($gene, -1);

		say $out join( "\t", ($geneName, $organism, map{ $type_table_ref->{$gene}{$_} // 0 } @indices ))
	}
	close $out;

}
