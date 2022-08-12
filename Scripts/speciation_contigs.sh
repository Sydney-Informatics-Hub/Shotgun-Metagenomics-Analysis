#!/bin/bash

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


set -e

sample=`echo $1 | cut -d ' ' -f 1`

database=./kraken2_standard_db_build_25-06-2022 # UPDATE THIS PATH TO REFLECT YOUR DATABASE

input=./Assembly/${sample}/${sample}.filteredContigs.fa
final_outdir=./Speciation_contigs/${sample}
outdir=${PBS_JOBFS}/${sample}
mkdir -p ${outdir}
mkdir -p ${final_outdir}

out=${outdir}/${sample}.kraken2.standard.out
report=${outdir}/${sample}.kraken2.standard.report
krona=${outdir}/${sample}.krona.html


kraken2 \
	--db $database \
	--threads $NCPUS \
	--output $out \
	--report $report \
	--report-zero-counts \
	--use-names \
	$input


ktImportTaxonomy \
	-m 3 -t 5 \
	$report \
	-o $krona

 
tar -zcf ${final_outdir}/${sample}.krona.html.files.tar.gz ${outdir}/${sample}.krona.html.files #tens of Ks of these per sample, pushing iNode over quota

cp $out $final_outdir
cp $report $final_outdir
cp $krona $final_outdir


