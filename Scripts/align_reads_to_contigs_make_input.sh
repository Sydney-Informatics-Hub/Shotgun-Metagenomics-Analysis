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

config=./Inputs/<cohort>.config
inputs=./Inputs/align_reads_to_contigs.inputs

rm -f $inputs

awk 'NR>1' ${config} | while read LINE
do 
        sample=`echo $LINE | cut -d ' ' -f 1`
        labSampleID=`echo $LINE | cut -d ' ' -f 2`
        platform=`echo $LINE | cut -d ' ' -f 3`
	centre=`echo $LINE | cut -d ' ' -f 4`
        lib=1
	
	printf "${labSampleID},${platform},${centre},${lib}\n" >> $inputs
done

tasks=`wc -l < $inputs`
printf "Number of alignment tasks to run: ${tasks}\n"

