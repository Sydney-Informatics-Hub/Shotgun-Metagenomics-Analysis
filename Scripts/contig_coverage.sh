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

sample=`echo $1 | cut -d ',' -f 1` 
mapq=`echo $1 | cut -d ',' -f 2` # minimum mapping quality, -q
baseQ=`echo $1 | cut -d ',' -f 3` # minimum base quality, -Q 

assdir=./Assembly
outdir=./Align_to_assembly

contigs=${assdir}/${sample}/${sample}.contigs.fa
sampledir=${outdir}/${sample}
bam=${sampledir}/${sample}.sort.bam

samtools coverage \
	-q ${mapq} -Q ${baseQ} \
	--reference ${contigs} \
	-o ${sampledir}/${sample}.q${mapq}Q${baseQ}.cov \
	${final}




