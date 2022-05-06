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

my $cohort = '<cohort>'; 

#1) Collect ARG data from curated list: 
my $curated = './Inputs/curated_ARGs_list.txt';
print "Collecting curated ARGs from $curated\n"; 

my $genehash = {}; 
my $dup = 0;
open (I, $curated) || die "$! $curated\n"; 
while (my $line = <I>) {
	chomp $line;
	if ($line!~/^\#/) {
		my ($gene, $family, $class, @vars) = split(' ', $line);
		if (! $genehash->{$gene}) {
			$genehash->{$gene}->{family} = $family; 
			$genehash->{$gene}->{class} = $class;
			$genehash->{$gene}->{id} = $gene; # the gene name to use for all variant names of this gene
			foreach my $var (@vars) {
				$genehash->{$var}->{family} = $family; 
				$genehash->{$var}->{class} = $class;
				$genehash->{$var}->{id} = $gene;			
			}
		} 
		else {
			print "Duplicate entries for $gene\n"; 
			$dup++; 
		}
	}	
} close I;

if ( $dup) {
	print "Please correct gene duplicate entries and resubmit. Terminating.\n"; die; 
} 
 

#2) Report curated ARGs from abricate output and link with species from kraken contig output:
# Print per sample:
my $list = "./Inputs/$cohort\_samples.list";
open (L, $list) || die "$! $list\n";
chomp (my @samples = <L>);
close L;
my $all_cat = ''; 
my $specieshash = {};
my $countshash = {}; 
foreach my $sample (@samples) {
	# Collect species ID of contigs:
	print "Collecting species for $sample...\n"; 
	my $kraken = "./Speciation_contigs/$sample\/$sample\.kraken2.standard.out";
	undef %{$specieshash}; 
	undef %{$countshash}; 
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
	
	# Collect counts
	print "Collecting normalised and raw ARG read counts for $sample...\n";
	my $counts = "./ARGs/ARG_read_counts/$sample\.curated_ARGs.counts.norm";
	open (C, $counts) || die "$! $counts\n"; 
	chomp (my $header = <C>);
	while (my $line = <C>) {
		chomp $line;
		my @cols = split ('\t', $line); 
		my $gene = $cols[1];
		my $contig = $cols[2];
		my $start = $cols[3];
		my $end = $cols[4];
		my $strand = $cols[5];
		my $count = $cols[7];
		my $tpm = $cols[11];
		my $rpkm = $cols[14];
		$countshash->{"$gene\_$contig\_$start\_$end\_$strand"}->{count} = $count;
		$countshash->{"$gene\_$contig\_$start\_$end\_$strand"}->{tpm} = $tpm;
		$countshash->{"$gene\_$contig\_$start\_$end\_$strand"}->{rpkm} = $rpkm;	
	} close C; 
	
	print "Writing collated ARG and species data for $sample...\n"; 
	my $args = "./ARGs/Abricate/$sample\/$sample\.ARGs.txt"; # all databases collated previosuly into one TSV per sample
	open (I, $args) || die "$! $args\n";
	my $out = "./ARGs/Curated_ARGs/$sample\.curated_ARGs.txt";
	$all_cat .= " $out"; 
	open (OUT, ">$out") || "$! write $out\n"; 
	print OUT "#Sample\tGene\tResistance_mechanism\tDrug_class\tSpecies\tContig\tStart\tEnd\tStrand\tRaw_count\tTPM\tRPKM\tCoverage\tCoverage_map\tGaps\t%Coverage\t%Identity\tDatabase\tAccession\tProduct\tResistance\n"; 	
	chomp ($header = <I>); 
	while (my $line = <I>) {
		chomp $line;
		#FILE   SEQUENCE(contig)        START   END     STRAND  GENE    COVERAGE        COVERAGE_MAP    GAPS    %COVERAGE       %IDENTITY       DATABASE        ACCESSION       PRODUCT RESISTANCE	
		my @cols = split('\t', $line); 
		my $contig = $cols[1];
		my $species = $specieshash->{$contig}->{species};
		my $gene = $cols[5];
		if ($genehash->{$gene}) {
			my $id = $genehash->{$gene}->{id};
			# Keep all the columns from ARGs input, BUT replace the gene column with the curated gene name
			# and put the gene in column 2, as well as the resistance mechanism (family) and drug class from args list and counts from hash 
			my ($file, $contig, $start, $end, $strand, $gene, @rest) = @cols; 
			my $count = $countshash->{"$id\_$contig\_$start\_$end\_$strand"}->{count};
			my $tpm = $countshash->{"$id\_$contig\_$start\_$end\_$strand"}->{tpm};
			my $rpkm = $countshash->{"$id\_$contig\_$start\_$end\_$strand"}->{rpkm};
			print OUT "$sample\t$id\t$genehash->{$gene}->{family}\t$genehash->{$gene}->{class}\t$species\t$contig\t$start\t$end\t$strand\t$count\t$tpm\t$rpkm";
			
			my $columns = @rest; 
			for ( my $i = 0; $i < $columns; $i ++ ) {
					print OUT "\t$rest[$i]";
			} print OUT "\n"; 
		}
	} close I; close OUT; 
}

# Print per cohort:	
my $all_out = "./ARGs/Curated_ARGs/$cohort\_allSamples.curated_ARGs.txt";
open (A, ">$all_out") || die "$! write $all_out\n"; 
print A "#Sample\tGene\tResistance_mechanism\tDrug_class\tSpecies\tContig\tStart\tEnd\tStrand\tRaw_count\tTMP\tRPKM\tCoverage\tCoverage_map\tGaps\t%Coverage\t%Identity\tDatabase\tAccession\tProduct\tResistance\n"; 	
close A; 
`cat $all_cat | sed '/^\#/d' | sort | uniq >> $all_out`; 
