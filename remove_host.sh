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

prefix=$1

fq1=$(ls ${prefix}*_R1*.f*q.gz)
fq2=$(ls ${prefix}*_R2*.f*q.gz)

prefix=$(basename $prefix)

out=./Target_reads/${prefix}.interleaved.extracted.fq.gz
log=./Logs/Remove_host/${prefix}.bbmap.log

# BBmap will look in ./ref for the reference index (created by bbmap_prep.pbs)
# Script tool is  used to force log to end up in the right place

script -c "bbmap.sh \
	-Xmx42g \
	in=${fq1} \
	in2=${fq2} \
	outu=${out} \
	threads=$NCPUS \
	overwrite=t \	
	usejni=t \
	unpigz=t" ${log}
