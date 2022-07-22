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

# GFF format:
#Contig source  feature start   end     score   strand  frame   attribute
#k141_105473     Prodigal:002006 CDS     202     921     .       -       0       ID=ANKLIJNP_00066;Parent=ANKLIJNP_00066_gene;Name=IS4811_IS1031_IS5_ORF_1;gene=IS4811_IS1031_IS5_ORF_1;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:IS.faa:IS4811;locus_tag=ANKLIJNP_00066;product=Passenger Gene

# abricate tab format:
#FILE   SEQUENCE        START   END     STRAND  GENE    COVERAGE        COVERAGE_MAP    GAPS    %COVERAGE       %IDENTITY       DATABASE        ACCESSION       PRODUCT RESISTANCE
#T1121B_21092006_D1_D2.filteredContigs.fa        k141_10276      685     812     -       tet(M)_1        1793-1920/1920  ..............= 0/0     6.67    100.00  resfinder       X92947  tet(M)  Doxycycline;Tetracycline;Minocycline

# columns to print:
# 1 contig = col 2 from abricate
# 2 source = abricate
# 3 feature = gene
# 4 start = col 3
# 5 end = col 4
# 6 score = .
# 7 strand = col 5
# 8 frame = ? ( 0 or . )
# 9 attribute = Name=<gene>;gene=<gene>

`mkdir -p ./ARGs/Curated_GFF`;

#) Collect ARG data from curated list:
my $curated = './Inputs/curated_ARGs_list.txt';
my $genehash = {}; 
open (I, $curated) || die "$! $curated\n"; 
my $dup = 0; 
while (my $line = <I>) {
	chomp $line;
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
		$dup ++; 
	}	
} close I; 

if ( $dup) {
	print "Please correct gene duplicate entries and resubmit. Terminating.\n"; die; 
} 


#) Open the abricate output, and for any entry that hits the curated genes or their variant names, output as GFF:
my $list = './Inputs/test_samples.list';
open L, '<', $list;
chomp (my @samples = <L>);
close L;

foreach my $sample (@samples) {
	my $args = "./ARGs/Abricate/$sample\/$sample\.ARGs.txt"; 
	my $out = "./ARGs/Curated_GFF/$sample\.curated_ARGs.gff";
	print "Writing $out\n"; 
	open (I, $args) || die "$! $args\n";
	my @lines = (); 
	chomp (my $header = <I>); 
	while (my $line = <I>) {
		chomp $line;
		#FILE   SEQUENCE(contig)        START   END     STRAND  GENE    COVERAGE        COVERAGE_MAP    GAPS    %COVERAGE       %IDENTITY       DATABASE        ACCESSION       PRODUCT RESISTANCE	
		my @cols = split('\t', $line); ; 
		my $gene = $cols[5];
		if ($genehash->{$gene}) { 
			my $id = $genehash->{$gene}->{id}; 
			my $unique_id = "$id\:$cols[1]\:$cols[2]\:$cols[3]\:$cols[4]"; 
			my $string_to_print = "$cols[1]\tabricate\tgene\t$cols[2]\t$cols[3]\t.\t$cols[4]\t.\tID=$unique_id\;Name=$id\;gene=$id\n";
			push @lines, $string_to_print; 
		}
		else {
			print "WARN: Could not find $gene from $args in ../Curated_ARGs/curated_ARGs_list.txt\n";
		}
	} close I;

	# ) Print out unique entries
	open (OUT, ">$out") || "$! write $out\n"; 
	my @unique = uniq( @lines );
	foreach my $line (@unique) {
		print OUT "$line"; 
	} close OUT; 
}	

#------- 
sub uniq {
	my %seen;
	return grep { !$seen{$_}++ } @_;
}
#
