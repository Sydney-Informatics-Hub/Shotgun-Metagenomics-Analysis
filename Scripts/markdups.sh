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

sample=$1

assdir=./Align_to_assembly/${sample}
bam=${assdir}/${sample}.sort.bam
dedup=./${assdir}/${sample}.sort.dedup.bam
metrics=${assdir}/${sample}.dedup_metrics.txt 
log=./Logs/Mark_dups/${sample}.dedup.log

rm -rf $dedup $metrics $log

# Max java heap space, 30 GB per CPU per task (ie, Gadi hugemem nodes)
mem=$(expr $NCPUS \* 30)

gatk MarkDuplicates \
	--java-options "-Xmx${mem}G -Dsamjdk.compression_level=5 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	-I ${bam} \
	-O ${dedup} \
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
	-M ${metrics} > ${log} 2>&1 

samtools index -@ ${NCPUS} ${dedup}
