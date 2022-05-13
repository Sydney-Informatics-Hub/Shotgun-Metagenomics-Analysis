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
my $config = "./Inputs/$cohort\.config";

my $indir = './ARGs/Curated_ARGs';
my $raw_input = "$indir\/$cohort\_allSamples_curated_ARGs_rawCount_Rdataframe.txt";
my $norm_input = "$indir\/$cohort\_allSamples_curated_ARGs_TPM_Rdataframe.txt";

# Collect samples by group IDs
my $grouphash = {};
my @groups = '';
open (S, $config) || die "$! $config\n"; 
chomp (my $header = <S>); 
while (my $line = <S>) {
	chomp $line;
	my ($id, $sample, $platform, $centre, $group) = split(' ', $line); 
	if (!$grouphash->{$group}) {
		push @groups, $group; 	
	}
	$grouphash->{$group}->{$sample} = 1; 
} close S; 

# Print per-group output files
foreach my $group (sort keys %{$grouphash}) {
        print "Printing per-group output for group $group\:\n";
        my $group_temp = "./$group\_ID_list.temp";
	open (G, ">$group_temp") || die "$! $group_temp\n"; 
        foreach my $sample (sort keys %{$grouphash->{$group}}) {
                print G "$sample\n"; 	
        } close G; 

	# TPM: 	
        my $group_out = $norm_input;
        $group_out =~ s/allSamples/$group/;
	print "\t$group_out\n";	
	`head -1 $norm_input > $group_out`;
	`join -t \$'\t' $group_temp $norm_input >> $group_out`;
	
	# Raw counts:         
	$group_out = $raw_input;
        $group_out =~ s/allSamples/$group/;
	print "\t$group_out\n";	
	`head -1 $raw_input > $group_out`;
	`join -t \$'\t' $group_temp $raw_input >> $group_out`;	
	
	# Delete temp file:
	unlink($group_temp);
}
	
