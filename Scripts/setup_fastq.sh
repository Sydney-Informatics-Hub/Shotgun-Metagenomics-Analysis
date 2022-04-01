#! /bin/bash

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

path=$1

fastqs=$(find $path -type f -name "*f*q.gz")
fastqs=($fastqs)

echo "$(date): Found ${#fastqs[@]} *fastq.gz files. Creating symlinks in ./Fastq directory"

for fastqpath in ${fastqs[@]}
do
	fastq=$(basename ${fastqpath})
	ln -f -s $fastqpath ./Fastq/$fastq
done


