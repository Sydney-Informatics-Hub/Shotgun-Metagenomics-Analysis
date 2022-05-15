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

fastq=`echo $1 | cut -d ',' -f 1`

prefix=$(basename $fastq | sed  's/.f.*q.gz$//g')
out=./FastQC/${prefix}
log=./Logs/FastQC/${prefix}.log


fastqc --extract -o ${out} ${fastq} >> ${log} 2>&1 
