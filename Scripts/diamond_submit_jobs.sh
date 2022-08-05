#!/bin/bash

#########################################################
#
# Platform: NCI Gadi HPC
# Description: see https://github.com/Sydney-Informatics-Hub/Shotgun-Metagenomics-Analysis
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

set -e

inputs=./Inputs/<cohort>_samples.list

logdir=./Logs/Diamond_NCBI
script=./Scripts/diamond_taxon.pbs

mkdir -p ${logdir}
mkdir -p ./Diamond_NCBI

while read sample
do
        echo Submitting diamond for $sample
        qsub -N d-${sample} -o ${logdir}/${sample}.diamond.o -e ${logdir}/${sample}.diamond.e -v sample="$sample" ${script}
        sleep 1
        echo
done < $inputs
