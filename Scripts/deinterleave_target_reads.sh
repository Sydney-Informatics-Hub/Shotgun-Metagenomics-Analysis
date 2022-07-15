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

sample=`echo $1 | cut -d ' ' -f 1`

indir=./Target_reads
outdir=./Target_reads_paired

fq=$(ls ./${indir}/${sample}*fq.gz)
fq=($fq)

for ((i = 0; i < ${#fq[@]}; i++ ))
do
	file=$(basename ${fq[$i]})
	r1=${file/.interleaved.extracted.fq.gz/_R1_paired.extracted.fq}
	r2=${file/.interleaved.extracted.fq.gz/_R2_paired.extracted.fq}

	reformat.sh \
		in=${fq[$i]} \
		ow=t \
		out1=${outdir}/${r1} \
		out2=${outdir}/${r2}
	
	\rm -rf ${outdir}/${r1}.gz ${outdir}/${r2}.gz
	pigz -p $NCPUS ${outdir}/${r1}
	pigz -p $NCPUS ${outdir}/${r2}

done


