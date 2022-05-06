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

sample=$1

bam=./Align_to_assembly/${sample}/${sample}.sort.dedup.bam
gff=./ARGs/Curated_GFF/${sample}.curated_ARGs.gff
logfile=./Logs/HTSeq/${sample}.log
out=./ARGs/ARG_read_counts/${sample}.curated_ARGs.counts

rm -rf ${logfile} ${out}

$HOME/.local/bin/htseq-count \
        --order pos \
        --stranded no \
        --minaqual 20 \
        --type gene \
        --idattr ID \
        --mode union \
        --nonunique none \
        --secondary-alignments ignore \
        --supplementary-alignments score \
        --counts_output ${out} \
        ${bam} ${gff} 2>>${logfile}
