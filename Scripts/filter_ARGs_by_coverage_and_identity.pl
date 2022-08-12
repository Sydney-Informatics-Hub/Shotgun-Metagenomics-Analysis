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

my $cohort='<cohort>';

# User-specified parameters for filtering:
my $cover = 70; # Percentage coverage threshold for filtering - adjust as required
my $identity = 80; # Percentage identity threshold for filtering - adjust as required

# Gather curated genes
my $curated = './Inputs/curated_ARGs_list.txt';
print "Collecting curated ARGs from $curated\n"; 
my @genes = (); 
open (C, $curated) || die "$! $curated\n";
chomp (my $header = <C>);
while (my $line = <C>) {
	my ($gene, @rest) = split('\t', $line); 
	push @genes, $gene; 
} close C; 

# Collect data
my $in = "./ARGs/Curated_ARGs/$cohort\_allSamples.curated_ARGs.txt";
my $storehash = {}; 

open (I, $in) || die "$! $in\n";
chomp ($header = <I>); 
while (my $line = <I>) {
	my @cols = split('\t', $line); 
	my $sample = $cols[0];
	my $gene = $cols[1];
	my $raw_count = $cols[9];
	my $tpm = $cols[10];
	my $cover = $cols[15];
	my $identity = $cols[16];
	
	if ( ($cover >= $cover) && ($identity >= $identity) ) { # Filter for coverage and identity percent
		if ($storehash->{$sample}->{$gene}) { # gene has been seen before - sum them		
			my $sum_raw = $raw_count + $storehash->{$sample}->{$gene}->{raw_count};
			my $sum_tpm = $tpm + $storehash->{$sample}->{$gene}->{tpm}; 
			$storehash->{$sample}->{$gene}->{raw_count} = $sum_raw;
			$storehash->{$sample}->{$gene}->{tpm} = $sum_tpm; 
		}
		else {
			$storehash->{$sample}->{$gene}->{raw_count} = $raw_count;
			$storehash->{$sample}->{$gene}->{tpm} = $tpm;
		}
	}
} close I;

# Add headers to outfiles
my $out_raw = "./ARGs/Curated_ARGs/$cohort\_allSamples_curated_ARGs_rawCount_Rdataframe.txt";
my $out_tpm = "./ARGs/Curated_ARGs/$cohort\_allSamples_curated_ARGs_TPM_Rdataframe.txt";
print "Printing raw counts to $out_raw\n";
print "Printing TPM-normalised counts to $out_tpm\n";
open (OR, ">$out_raw") || die "$! write $out_raw\n"; 
open (OT, ">$out_tpm") || die "$! write $out_tpm\n"; 
print OR "Sample_ID";
print OT "Sample_ID";
foreach my $gene (@genes) {
	print OR "\t$gene";
	print OT "\t$gene";
} 
print OR "\n"; 
print OT "\n";	

# Print counts to outfiles	
foreach my $sample (sort keys %{$storehash}) {
	print OR "$sample"; 
	print OT "$sample"; 
	foreach my $gene (@genes) { #use array format to ensure sort order in output
		if ($storehash->{$sample}->{$gene}) {
			print OR "\t$storehash->{$sample}->{$gene}->{raw_count}";
			print OT "\t$storehash->{$sample}->{$gene}->{tpm}";
		}
		else {
			print OR "\tNA";
			print OT "\tNA"; 
		}	
	}
	print OR "\n";
	print OT "\n";	
} close OR; close OT; 			
