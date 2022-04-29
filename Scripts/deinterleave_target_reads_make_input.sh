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
inputs=./Inputs/deinterleave_target_reads.inputs

rm -rf ${inputs} ${inputs}-unsorted

awk 'NR>1' ${config} | while read LINE
do 
	sample=`echo $LINE | cut -d ' ' -f 1`
	bytes=$(find ./Target_reads -type f -name "*${sample}*.fq.gz" -ls | awk '{total += $1} END {print total}')
	printf "${sample}\t${bytes}\n" >> ${inputs}-unsorted
done

# Reverse numeric sort on bytes
sort -rnk 2  ${inputs}-unsorted > ${inputs}
rm -rf ${inputs}-unsorted

tasks=`wc -l < $inputs`
printf "Number of deinterleave tasks to run: ${tasks}\n"

