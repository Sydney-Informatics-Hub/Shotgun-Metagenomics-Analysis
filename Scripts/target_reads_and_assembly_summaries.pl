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

my $sample=$ARGV[0];

my $out = "./Assembly/$sample\/$sample\.summary.txt";
open (O, ">$out") || die "$! write $out\n";
 
#Get raw reads from FastQC:
my @pairs = split(' ', `ls ./Fastq/$sample*_R1*.fastq.gz`); 
my $raw_reads = 0; 
foreach my $pair (@pairs) {
	chomp $pair; 
	$pair=`basename $pair`; 
	chomp $pair;  
	my $qc_dir = "./FastQC/$pair"; 
	$qc_dir =~s/\.fastq.gz/\_fastqc/; 
	my $qc_data = "$qc_dir\/fastqc_data.txt";   
	my $r1 = `grep "Total Sequences" $qc_dir\/fastqc_data.txt | awk '{print \$(NF)}'`;
	chomp $r1;
	my $raw_reads_pair = 2 * $r1;
	$raw_reads+=$raw_reads_pair;  
}
 
#Get target reads from zcat | wc -l (does not seem clear cut from megahit log)
my $lines = `zcat ./Target_reads/*$sample\* | wc -l`;
chomp $lines;
my $target_reads = $lines/4;
	
#Calculate % host contamination:
my $host_reads = $raw_reads - $target_reads; 
my $host = sprintf("%.3f",($host_reads/$raw_reads*100)); 

#Get megahit summary:
my $log = "./Assembly\/$sample\/$sample\.log"; 
my $data = `tail -2 $log | head -1`;
chomp $data;
#2019-10-16 12:38:07 - 97037 contigs, total 161134802 bp, min 200 bp, max 467257 bp, avg 1660 bp, N50 6086 bp	
my @cols = split(' ', $data); 
my $contigs = $cols[3];
my $tot_bp = $cols[6];
my $max_bp = $cols[12];
my $av_bp = $cols[15];
my $N50 = $cols[18];
	
#Filtered contigs count:
my $filt = `grep ">" ./Assembly\/$sample\/$sample\.filteredContigs.fa | wc -l`;
chomp $filt;
	
print O "$sample\t$raw_reads\t$target_reads\t$host\t$contigs\t$tot_bp\t$max_bp\t$av_bp\t$N50\t$filt\n"; 	 
close O; 
