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

# Collect the data into structure: 
print "Collecting counts into memory\n"; 
my $list = "./Inputs/$cohort\_samples.list";
open (L, $list) || die "$! $list\n"; 
my $tpmhash = {}; 
while (my $sample = <L>) {
	chomp $sample;
	my $counts = "./ARGs/ARG_read_counts/$sample\.curated_ARGs.reformat.counts";
	open (C, $counts) || die "$! $counts\n"; 
	chomp (my $header = <C>); #did not use # header prefix for compatibility with R
	my $counter = 0; 
	while (my $line = <C>) {
		chomp $line;
		if ($line!~m/^\_/) { 
			$counter++; 
			my ($id, $length, $count) = split (' ', $line); 
			my $rpk = ($count/($length/1000));
			$tpmhash->{$sample}->{$counter}->{$id}->{rpk} = $rpk; # use counter to retain print order
			$tpmhash->{$sample}->{$counter}->{$id}->{length}= $length;
			$tpmhash->{$sample}->{$counter}->{$id}->{count}= $count;
		}
	} close C; 
} close L; 

# Get the per sample RPK: 
my $all_cat = '';
my $info = "./$cohort\_target_reads_and_assembly_summary.txt";
foreach my $sample (keys %{$tpmhash}) {
	print "\tCalculating TPM for sample $sample\n"; 
	my $out = "./ARGs/ARG_read_counts/$sample\.curated_ARGs.counts.norm"; 
	$all_cat .= " $out"; 
	my $rpk_sum = 0;
	foreach my $counter (sort by_number keys %{$tpmhash->{$sample}}) {
		foreach my $id (keys %{$tpmhash->{$sample}->{$counter}}) {
			$rpk_sum += $tpmhash->{$sample}->{$counter}->{$id}->{rpk};
		}
	}
	my $tpm_scale = $rpk_sum / 1000000;

	# Get the per gene tpm:	
	foreach my $counter (sort by_number keys %{$tpmhash->{$sample}}) {
		foreach my $id (keys %{$tpmhash->{$sample}->{$counter}}) {
			my $tpm = $tpmhash->{$sample}->{$counter}->{$id}->{rpk} / $tpm_scale; 
			$tpmhash->{$sample}->{$counter}->{$id}->{tpm} = $tpm; 
		}
	}
	
	# Now do for RPKM/FPKM: 
	#need total fragments
	print "\tCalculating RPKM for sample $sample\n";
	open (O, ">$out") || die "$! write $out\n";  
	print O "#Sample\tGene\tContig\tStart\tEnd\tStrand\tLength\tRaw_counts\tRPK\tRPK_sum\tTPM_scale\tTPM\tLibrarySize\tRPKM_scale\tRPKM\n"; 
	my $reads = `grep $sample $info | awk '{print \$3}'`;
	chomp $reads;
	my $fragments = $reads / 2; 
	my $rpkm_scale = $fragments / 1000000;
	foreach my $counter (sort by_number keys %{$tpmhash->{$sample}}) {
		foreach my $id (keys %{$tpmhash->{$sample}->{$counter}}) {
			my ($gene, $contig, $start, $end, $strand) = split ('\:', $id); 
			my $count = $tpmhash->{$sample}->{$counter}->{$id}->{count};
			my $length = $tpmhash->{$sample}->{$counter}->{$id}->{length}; 			
			my $rpkm = $count / ( ($length / 1000) * ($fragments / 1000000) );  
			print O "$sample\t$gene\t$contig\t$start\t$end\t$strand\t$length\t$count\t$tpmhash->{$sample}->{$counter}->{$id}->{rpk}\t$rpk_sum\t$tpm_scale\t$tpmhash->{$sample}->{$counter}->{$id}->{tpm}\t$fragments\t$rpkm_scale\t$rpkm\n";			
		}
	}	
}

# Collate all in cohort: 	
my $all_out = "./ARGs/ARG_read_counts/$cohort\.curated_ARGs.counts.norm";
print "Collating all in cohort $cohort to $all_out\n";
open (A, ">$all_out") || die "$! write $all_out\n";  
print A "#Sample\tGene\tContig\tStart\tEnd\tStrand\tLength\tRaw_counts\tRPK\tRPK_sum\tTPM_scale\tTPM\tLibrarySize\tRPKM_scale\tRPKM\n"; 
close A; 
`cat $all_cat | sed '/^\#/d' | sort | uniq >> $all_out`; 

#-----------------------------------------------
sub by_number {
	$a<=>$b; 
}
