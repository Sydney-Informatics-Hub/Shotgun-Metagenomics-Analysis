#!/bin/bash

#########################################################
#
# Platform: NCI Gadi HPC
# Description: see https://github.com/Sydney-Informatics-Hub/Shotgun-Metagenomics-Analysis
#
# Author/s: Cali Willet; Tracy Chew
# cali.willet@sydney.edu.au; tracy.chew@sydney.edu.au
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

config=./Inputs/<cohort>.config)

logdir=./Logs/humann2
mkdir -p ${logdir}
mkdir -p Functional_profiling

awk 'NR>1' ${config} | while read LINE
do 
	sample=`echo $LINE | cut -d ' ' -f 2`
	echo Submitting humann2 for $sample
	qsub -N h2-${sample} -o ${logdir}/${sample}.o -e ${logdir}/${sample}.e -v sample="$sample" ./Scripts/functional_profiling.pbs
	sleep 1
	echo
done
