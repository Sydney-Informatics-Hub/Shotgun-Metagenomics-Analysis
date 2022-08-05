#!/bin/bash

#########################################################
#
# Platform: NCI Gadi HPC
# Description: see https://github.com/Sydney-Informatics-Hub/Shotgun-Metagenomics-Analysis
# 
# Author/s: Tracy Chew
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

sample=$1

contigs=./Assembly/${sample}/${sample}.filteredContigs.fa
outdir=./Prodigal_CDS

prodigal -p meta \
        -i ${contigs} \
        -f gff \
	-q \
        -a ${outdir}/${sample}.CDS.prot.fa \
        -o ${outdir}/${sample}.CDS.gff

