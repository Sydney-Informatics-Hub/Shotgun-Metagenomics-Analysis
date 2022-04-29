#!/usr/bin/env perl

#########################################################
#
# Platform: NCI Gadi HPC
# Description: see https://github.com/Sydney-Informatics-Hub/Shotgun-Metagenomics-Analysis
#
# Author/s: Cali Willet
# cali.willet@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the scientific and technical assistance
# <or e.g. bioinformatics assistance of <PERSON>> of Sydney Informatics
# Hub and resources and services from the National Computational
# Infrastructure (NCI), which is supported by the Australian Government
# with access facilitated by the University of Sydney.
#
#########################################################

use warnings;
use strict;


# per sample:
my $cohort = '<cohort>';
my $config = "./Inputs/$cohort\.config";
my @samples = split(' ', `awk 'NR>1 {print \$2}' $config`); 

my $all_cat = ''; 
my $specieshash = {};

my $indir = './Insertion_sequences'; 
my $outdir = "$indir\/Filtered_IS_with_species"; 
`mkdir -p $outdir`; 

foreach my $sample (@samples) {
	# Collect species ID of contigs:
	print "Collecting species for $sample...\n"; 
	my $kraken = "./Speciation_contigs/$sample\/$sample\.kraken2.standard.out";
	undef %{$specieshash}; 
	open (I, $kraken) || die "$! $kraken\n";
	while (my $line = <I>) {
		chomp $line;
		my ($category, $contig, $species, @rest) = split('\t', $line);  
		if (! $specieshash->{$contig}) {
			$specieshash->{$contig}->{species} = $species; 
		} 
		else {
			print "Duplicate entries for $contig - terminating, please fix this script\n"; die; 
		}
	} close I;	
	
	# Now gather IS sequences from Prokka:
	print "Writing collated IS and species data for $sample...\n"; 
	my $prokka = "$indir\/$sample\/$sample\.gff";
	open (I, $prokka) || die "$! $prokka\n";
	my $out = "$outdir\/$sample\.IS.txt";
	$all_cat .= " $out"; 
	open (OUT, ">$out") || "$! write $out\n"; 
	print OUT "#Sample\tContig\tSpecies\tStart\tEnd\tName\tInference\tProduct\n";  
	while (my $line = <I>) {
		chomp $line;
		#Contig source  feature start   end     score   strand  frame   attribute
		if ( ($line!~m/^\#/) && ($line =~m/Passenger|Transposase/) ) {	
			my @cols = split('\t', $line); 
			my $contig = $cols[0];
			my $start = $cols[3];
			my $end = $cols[4];
			my $species = $specieshash->{$contig}->{species}; 
			my @attributes = split(';', $cols[8]);
			my ($type, $name, $inference, $product);  
			foreach my $attr (@attributes) {
				if ($attr=~m/^Name=/) {
					($type, $name) = split('=', $attr);  
				}
				elsif ($attr=~m/^inference=/) {
					($type, $inference) = split('=', $attr);
				}
				elsif ($attr=~m/^product=/) {
					($type, $product) = split('=', $attr);
				}
			} 
			print OUT "$sample\t$contig\t$species\t$start\t$end\t$name\t$inference\t$product\n"; 
		
		}
	} close I; close OUT; 
}

# per cohort:	
my $all_out = "$outdir\/IS_$cohort\_allSamples.txt";
open (A, ">$all_out") || die "$! write $all_out\n"; 
print A "#Sample\tContig\tSpecies\tStart\tEnd\tName\tInference\tProduct\n"; 
close A; 
`cat $all_cat | sed '/^\#/d' | sort | uniq >> $all_out`; 
