#! /bin/bash

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

inputs=./Inputs/remove_host.inputs

find ./Fastq -name "*_R1*.f*q.gz" -exec stat -Lc '%s %n' {} \; | sort -rnk 1 | awk '{print $2}'| sed 's/_R1.*//' > $inputs

tasks=`wc -l < $inputs`

printf "Number of remove host tasks to run: ${tasks}\n"

