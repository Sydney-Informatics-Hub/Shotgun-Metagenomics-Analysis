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

cohort=<cohort>
config=./Inputs/${cohort}.config
out=./${cohort}_target_reads_and_assembly_summary.txt

printf "#Sample\tRaw_reads\tTarget_reads\t%%_host\tNum_contgs\tTot_bp\tMax_bp\tMean_bp\tN50\tFiltered_contigs\n" > $out

samples=$(awk 'NR>1 {print $2}' $config)
samples=($samples)

for ((i = 0; i < ${#samples[@]}; i++ ))
do 
	sample=${samples[$i]}
	summary=./Assembly/${sample}/${sample}.summary.txt
	cat $summary >> $out
	rm -rf $summary
done


