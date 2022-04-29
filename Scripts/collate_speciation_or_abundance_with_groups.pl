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


# Inputs:
my $config = './Inputs/<cohort>.config';

my $input = $ARGV[0]; # The speciation or abundance all-sample collated file, output of collate_speciation.pl or collate_abundance.pl
chomp $input; 


# Collect samples by group IDs
my $grouphash = {};
my @groups = '';
open (S, $config) || die "$! $config\n"; 
while (my $line = <S>) {
	my ($id, $sample, $platform, $centre, $group) = split(' ', $line); 
	if (!$grouphash->{$group}) {
		push @groups, $group; 	
	}
	$grouphash->{$group}->{$sample} = 1; 
} close S; 


# Header is needed to get the column numbers for each sample
open (O, $input)|| die "$! $input\n"; 
chomp (my $header = <O>); 
close O;
my @header = split('\t', $header); 


# Collect column numbers pertaining to samples for each group
foreach my $group (sort keys %{$grouphash}) {
	print "Collecting sample columns for group $group ";
	my $cut_list = '';
	foreach my $sample (sort keys %{$grouphash->{$group}}) {
		my ($index) = grep { $header[$_] eq $sample } (0 .. @header-1);
		$index++; 
		$cut_list.=",$index";  	
	} 
	$cut_list=~s/^\,/1,/; # make sure the first column with species names is printed to all groups outputs
 	# Parse the column numbers, as a comma delimited string list, to cut
	my $group_out = $input; 
	$group_out =~ s/allSamples/$group/;
	print "and printing to $group_out\n"; 
	`cut -f $cut_list $input > $group_out`; 
} 



