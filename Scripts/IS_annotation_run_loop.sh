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

logdir=./Logs/IS
mkdir -p ${logdir} ./Insertion_sequences

awk 'NR>1' ${config} | while read LINE
do 
        sample=`echo $LINE | cut -d ' ' -f 2`
        echo Submitting IS annotation for $sample
        qsub -N IS-${sample} -o ${logdir}/${sample}.o -e ${logdir}/${sample}.e -v sample="$sample" ./Scripts/IS_annotation.pbs
        sleep 1
        echo
done


