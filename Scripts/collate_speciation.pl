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

# Parse the name of the $set variable ('reads' or 'contigs') as first and only CLI argument

use warnings;
use strict;


my $set = $ARGV[0]; chomp $set; # do for contigs and reads separately
my $cohort = '<cohort>'; 
my $config = "./Inputs/$cohort\.config"; 
my $kdir = "./Speciation_$set";

my @samples = split(' ', `awk 'NR>1 {print \$2}' $config`); 

my $out_all = "$kdir\/Kraken2_$cohort\_$set\_allSamples.txt";

# Collect data to hash
my $krakhash = {};
my $sciformathash = {}; 
my $nameformathash = {};
my $first = 0; 
my $c = 0; 
foreach my $sample (@samples) {

	my $kreport = "$kdir/$sample\/$sample\.kraken2.standard.report";
	open (K, $kreport) || die "$! $kreport\n"; 
	while (my $line = <K>) {
		chomp $line; 
		my @cols = split('\t', $line);
		my $percent = $cols[0];
		my $sci_name = $cols[-1];
		$sci_name =~ s/^\s+|\s+$//g;
		$sci_name =~ s/\s+/\_/g;
		if (!$first) { #use the species order of the first sample to order the remainder
			$c++; 
			$sciformathash->{$c}->{sci_name} = $sci_name; 
		}
		$krakhash->{$sci_name}->{$sample}->{percent} = $percent; 	
		$nameformathash->{$first}->{sample_name} = $sample; 
		
	} close K; 
	$first++; 	
}

# Print headers:
open (O, ">$out_all") || die "$! write $out_all\n";
print O "#Sci_name";
foreach my $sample (sort by_number keys %{$nameformathash}) {
	print O "\t$nameformathash->{$sample}->{sample_name}";		 		
} print O "\n";

# Print output in order:
foreach my $c (sort by_number keys %{$sciformathash}) {
	my $sci_name = $sciformathash->{$c}->{sci_name}; 
	print O "$sci_name";
	foreach my $sample_num (sort by_number keys %{$nameformathash}) {
		my $sample = $nameformathash->{$sample_num}->{sample_name};
		my $percent = $krakhash->{$sci_name}->{$sample}->{percent};
		print O "\t$percent"; 						
	} print O "\n"; 		
} close O; 

sub by_number {
	$a<=>$b;
}
