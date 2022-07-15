#! /usr/bin/perl
#
#Sorry this script could've been much more efficient..

set -e

#########################################################
#
# Platform: NCI Gadi HPC
#
# Author: Tracy Chew
# tracy.chew@sydney.edu.au
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

use strict;
use warnings;

# Change path to your sample.list file. Used to obtain sampleIDs
my $inputs="./Inputs/samples.list";

# Save outputs. Change to your desired output names
my $outdir="../Diamond_NCBI_ARGs"; # also location of per sample output files
my $summary="$outdir\/Allsamples_resistome.txt";
`mkdir -p $outdir`;
`rm -rf $summary`;

# % RESISTOME AS: TOTAL ARGS / TOTAL GENES 
# Write headers
open(S,">>",$summary)||die "Could not write to $!\n";
print S "#SampleID\tTimepoint\tARG_matches\tNCBINR_matches\tFailed_length\t%Resistome\n";

open(I,"$inputs")||die "$!\n";
while(<I>){
	
	# Get files required for each sample
	chomp(my $sampleid=$_);
	my $diamond="../Diamond_NCBI/$sampleid\.DIAMOND.tsv";

	# THIS NEEDS TO BE CHANGED TO SUIT YOUR PROJECT 
	chomp(my $analysis_set=`grep $sampleid ./Inputs/ADL*/samples.list | awk -F'/' '{print \$3}'`);
	chomp(my $timepoint=`grep -l $sampleid ./Inputs/$analysis_set\/*_samples.list | grep -o "T."`);
	if(!$timepoint){ $timepoint = "NA"; }
	my $ARGs="../../$analysis_set\_analysis_Aug21/Curated_ARGs/$sampleid\.curated_ARGs.txt";
	my $prodigal="../Prodigal_CDS/$sampleid\.CDS.gff";

	open(S,">>",$summary)||die "Could not write to $!\n";

	# Save outputs
	my $out="$outdir\/$sampleid\.DIAMOND_ARGs.txt";

	print "Reading $sampleid curated ARGs with start and end positions derived from ABRICATE\n";
	# Save sample curated ARGs to memory
	my $arghash;
	my $sppwithARGhash;
	open(A,$ARGs) || die "Could not open $ARGs: $!\n";
	while(<A>){
		chomp $_;
		my(@info)=split("\t",$_);
		
		my $contig = $info[5];
		my $genename = $info[1];
		my $startpos = $info[6];
		my $endpos = $info[7];
		my $species = $info[4];
		(my $taxid) = ($species =~ m/taxid ([0-9]+)/g);
		#print "$contig\t$genename\t$startpos\t$endpos\t$species\t$taxid\n";

		$arghash->{$contig}->{gene} = $genename;
		$arghash->{$contig}->{startpos} = $startpos;
		$arghash->{$contig}->{endpos} = $endpos;
		$arghash->{$contig}->{species} = $species;
		
		# Get a unique record of species with an ARG
		$sppwithARGhash->{$taxid}->{species} = $species;
	}
	close A;

	# Save positions for each CDS on MEGAHIT contig from Prodigal output
	my $cdshash;

	print "Reading $sampleid Prodigal predicted CDS sequences with start and end positions\n";
	open(P,$prodigal)||die"$!\n";
	while(<P>){
		chomp $_;
		if($_!~/^#/){
			my(@info)=split(" ",$_);
			
			my $contig = $info[0];
			my $cdsstart = $info[3];
			my $cdsend = $info[4];
			my @vars = split(";",$info[8]);
			my $genenum = $vars[0];
			$genenum =~ s/^ID=[0-9+]*_//;
			
			my $geneid = "$contig\_$genenum";
	
			#print "$contig\t$geneid\t$cdsstart\t$cdsend\n";
	
			$cdshash->{$geneid}->{cdsstart} = $cdsstart;
			$cdshash->{$geneid}->{cdsend} = $cdsend;

		}
	}

	# Keep counts
	my $ARG_match=0;
	my $total_diamond=0;
	my $failed_length=0;
	
	print "Matching $sampleid top hit DIAMOND result per CDS to curated ARGs using contig and gene positions\n";
	# Print per sample results to file
	open(O,">$out")||die "Could not write to $out: $!\n";
	
	# Print header
	print O "#Contig\tProdigal_CDS\tSpecies\tGene\tARG_start\tARG_end\tProdigal_start\tProdigal_end\n";
	
	# Match NCBI NR genes identified by DIAMOND to ARG using contig
	# DIAMOND does not contain start and end positions - obtain this from Prodigal data
	open(D,$diamond)||die "Could not open $diamond: $!\n";
	while(<D>){
		chomp $_;
		$total_diamond++;
		my(@info)=split("\t",$_);
		
		my $length = $info[3];
		my $qstart = $info[6];	
		my $qend = $info[7];
		my $contig_gene = $info[0];
		
		# Can contain multiple ids
		my $diamond_taxid = $info[12];
		
		(my $query = $contig_gene) =~ s/_[0-9+]*$//g;
		
		# Filter by length
		if ( $length > 25 ){
			# Get gene positions on MEGAHIT contig from prodigal output
			# Start or end positions match
			if ($arghash->{$query} && ($arghash->{$query}->{startpos} == $cdshash->{$contig_gene}->{cdsstart} || $arghash->{$query}->{endpos} == $cdshash->{$contig_gene}->{cdsend}) ) {
				print O "$query\t$contig_gene\t$arghash->{$query}->{species}\t$arghash->{$query}->{gene}\t$arghash->{$query}->{startpos}\t$arghash->{$query}->{endpos}\t$cdshash->{$contig_gene}->{cdsstart}\t$cdshash->{$contig_gene}->{cdsend}\n";
				$ARG_match++;
				
			}
		}
		else{
			$failed_length++;
		}
	}

	my $passed_diamond=$total_diamond-$failed_length;
	my $resistome=$ARG_match/$passed_diamond*100;
	print S "$sampleid\t$timepoint\t$ARG_match\t$total_diamond\t$failed_length\t$resistome\n";
	close O;
}


