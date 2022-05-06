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

my $set = $ARGV[0]; chomp $set; # do reads and contigs separately

my $dir = "./Abundance_$set";
my $cohort = '<cohort>'; 
my $config = "./Inputs/$cohort\.config";
my @samples = split(' ', `awk 'NR>1 {print \$2}' $config`); 

my $out_all = "$dir\/Bracken2_$cohort\_$set\_allSamples.txt";

my $brakhash = {};
my $sciformathash = {}; 
my $nameformathash = {};
my $first = 0; 
my $c = 0; 

# Collect abundance data, store in hash so that the species order can be preserved when printing out
foreach my $sample (@samples) { 
	my $report = "$dir/$sample\/$sample\.bracken";
	open (K, $report) || die "$! $report\n"; 
	chomp (my $header = <K>); 
	while (my $line = <K>) {
		chomp $line; 
		my @cols = split('\t', $line);
		my $sci_name = "$cols[0]\_$cols[1]";
		$sci_name =~ s/^\s+|\s+$//g;
		$sci_name =~ s/\s+/\_/g;
		my $ab = $cols[-1]; 	 
		$brakhash->{$sci_name}->{$sample}->{abundance} = $ab; 	
		$nameformathash->{$first}->{sample_name} = $sample; 
		
	} close K; $first++; 	
}

# Print headers
open (O, ">$out_all") || die "$! write $out_all\n";
print O "#Species";
foreach my $sample (sort by_number keys %{$nameformathash}) {
	print O "\t$nameformathash->{$sample}->{sample_name}";		 		
} print O "\n"; 

foreach my $sci_name (sort keys %{$brakhash}) { 
	print O "$sci_name";
	foreach my $sample_num (sort by_number keys %{$nameformathash}) {
		my $sample = $nameformathash->{$sample_num}->{sample_name};
		my $ab = 0;
		if (	$brakhash->{$sci_name}->{$sample} ) {
			$ab = $brakhash->{$sci_name}->{$sample}->{abundance};
		}
		if ($sci_name =~m/sticklandii/) {
			print "$sample\t$sci_name\t$ab\n"; 
		}
		
		print O "\t$ab"; 						
	} print O "\n"; 		
} close O; 

sub by_number {
	$a<=>$b;
}
