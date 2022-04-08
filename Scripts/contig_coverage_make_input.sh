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

config=<config>
inputs=./Inputs/contig_coverage.inputs

rm -f $inputs

echo Minimum mapping quality to include for coverage calculation:
read mapq
echo

echo Minimum base quality to include for coverage calculation:
read baseQ
echo

echo Using minimum mapping quality $mapq and minimum base quality $baseQ for coverage calculation
echo

awk 'NR>1 {print $2}' $config > $inputs	

sed -i "s/$/,${mapq},${baseQ}/" $inputs

tasks=`wc -l < $inputs`
printf "Number of metagenomic assembly tasks to run: ${tasks}\n"

