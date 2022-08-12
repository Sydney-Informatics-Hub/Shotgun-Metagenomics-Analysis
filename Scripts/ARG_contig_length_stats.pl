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

# Collect ARG-containing contigs into hash structure: 
print "Collecting ARG-containing contigs\n"; 
my $args = "./ARGs/Curated_ARGs/$cohort\_allSamples.curated_ARGs.txt";
my $arg_contig_storehash = {};
open (S, $args) || die "$! $args\n"; 
chomp (my $header = <S>); 
while (my $line = <S>) {
	chomp $line; 
	my @cols = split('\t', $line);
	my $sample = $cols[0];
	my $contig = $cols[5]; 
	$arg_contig_storehash->{$sample}->{$contig} = 1; # Do not care how many times the contig has been seen so store arbitrary value for the contig key 	
} close S; 

# Collect contig lengths from assembly and compute descriptive stats:
print "Collect contig length data from assemblies\n"; 
my $samples = "./Inputs/$cohort\_samples.list";
my $samplehash= {}; 
open (S, $samples) || die "$! $samples\n"; 
while (my $sample = <S>) {
	chomp $sample;
	print "\tComputing descriptive contig length stats for $sample\n"; 
	my $assembly = "./Assembly/$sample\/$sample\.filteredContigs.fa"; 
	open (A, $assembly) || die "$! $assembly\n"; 
	
	my @all_contig_lengths = ();
	my $all_contig_count = 0;
	my $all_sum_contig_length = 0;
	
	my @arg_contig_lengths = ();
	my $arg_contig_count = 0;
	my $arg_sum_contig_length = 0;	
  
	while (my $line = <A>) {
		chomp $line;
		if ($line=~m/^\>/) {
			my @cols = split('=', $line); 
			my $length = $cols[-1]; # last field on the contig id line
			my $contig = $cols[0];
			$contig=~s/^\>//;
			$contig=~s/ flag$//;
			$all_contig_count++;
			$all_sum_contig_length += $length; 
			push @all_contig_lengths, $length; 
			
			if ($arg_contig_storehash->{$sample}->{$contig}) {
				$arg_contig_count++;
				$arg_sum_contig_length += $length; 
				push @arg_contig_lengths, $length; 				
			}
		}
	} close A; 
	
	# mean, min, max:
	my $all_mean = $all_sum_contig_length / $all_contig_count;
	my @all_sorted = sort { $a <=> $b } @all_contig_lengths; 
	my $all_min = $all_sorted[0];
	my $all_max = $all_sorted[-1]; 
	
	my $arg_mean = $arg_sum_contig_length / $arg_contig_count;
	my @arg_sorted = sort { $a <=> $b } @arg_contig_lengths; 
	my $arg_min = $arg_sorted[0];
	my $arg_max = $arg_sorted[-1]; 	
	
	# stdev:
	my $all_sum_a = 0;
	foreach my $len (@all_contig_lengths) {
		my $a = ($len - $all_mean)**2; 
		$all_sum_a += $a; 
	}
	my $all_b = $all_sum_a / $all_contig_count;
	my $all_stdev = sqrt($all_b); 
	
	my $arg_sum_a = 0;
	foreach my $len (@arg_contig_lengths) {
		my $a = ($len - $arg_mean)**2; 
		$arg_sum_a += $a; 
	}
	my $arg_b = $arg_sum_a / $arg_contig_count;
	my $arg_stdev = sqrt($arg_b); 	
	
	# save results
	$samplehash->{$sample}->{all_mean} = $all_mean;
	$samplehash->{$sample}->{all_min} = $all_min;	  
	$samplehash->{$sample}->{all_max} = $all_max;
	$samplehash->{$sample}->{all_stdev} = $all_stdev;
	$samplehash->{$sample}->{all_count} = $all_contig_count;
	
	$samplehash->{$sample}->{arg_mean} = $arg_mean;
	$samplehash->{$sample}->{arg_min} = $arg_min;	  
	$samplehash->{$sample}->{arg_max} = $arg_max;
	$samplehash->{$sample}->{arg_stdev} = $arg_stdev;
	$samplehash->{$sample}->{arg_count} = $arg_contig_count;		
}

# Print results for all samples 
my $out = "./ARGs/Curated_ARGs/$cohort\_allSamples_curated_ARGs_contig_length_stats.txt";
print "Printing results for all samples to $out\n"; 
open (O, ">$out") || die "$! write $out\n"; 
print O "#Sample\tAll_contigs_count\tAll_contigs_length_mean\tAll_contigs_length_stdev\tAll_contigs_length_min\tAll_contigs_length_max\tARG_contigs_count\tARG_contigs_length_mean\tARG_contigs_length_stdev\tARG_contigs_length_min\tARG_contigs_length_max\n"; 
foreach my $sample (sort keys %{$samplehash}) {
	print O "$sample\t$samplehash->{$sample}->{all_count}\t$samplehash->{$sample}->{all_mean}\t$samplehash->{$sample}->{all_stdev}\t$samplehash->{$sample}->{all_min}\t$samplehash->{$sample}->{all_max}";
	print O "\t$samplehash->{$sample}->{arg_count}\t$samplehash->{$sample}->{arg_mean}\t$samplehash->{$sample}->{arg_stdev}\t$samplehash->{$sample}->{arg_min}\t$samplehash->{$sample}->{arg_max}\n";  
} close O; 



