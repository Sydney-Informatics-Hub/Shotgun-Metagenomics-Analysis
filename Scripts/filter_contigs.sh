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

sample=`echo $1 | cut -d ',' -f 1` 
mapq=`echo $1 | cut -d ',' -f 2`
baseQ=`echo $1 | cut -d ',' -f 3`

cov_dir=./Align_to_assembly
ass_dir=./Assembly 

cov=${cov_dir}/${sample}/${sample}.q${mapq}Q${baseQ}.cov
list=${cov_dir}/${sample}/${sample}.filtIDs.list
contigs=${ass_dir}/${sample}/${sample}.contigs.fa
filtered_contigs=${ass_dir}/${sample}/${sample}.filteredContigs.fa

# Apply filtering criteria to the SAMtools 'coverage' file using awk (adjust the below to suit your desired filtering thresholds):
awk '$7>=1' $cov | awk '{print $1}' > $list

# Extract matching IDs to a new contigs fasta:
seqtk subseq $contigs $list > $filtered_contigs

# Clean up:
rm -rf $list
