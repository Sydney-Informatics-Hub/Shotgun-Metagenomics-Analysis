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

contigs=./Assembly/${sample}/${sample}.filteredContigs.fa
outdir=./ARGs/Abricate/${sample}
mkdir -p ./ARGs/Abricate ${outdir}

abricate --threads $NCPUS --nopath \
	--db ncbi ${contigs} > ${outdir}/${sample}.ncbi.tab
	
abricate --threads $NCPUS --nopath \
	--db resfinder ${contigs} > ${outdir}/${sample}.resfinder.tab
		
abricate --threads $NCPUS --nopath \
	--db card ${contigs} > ${outdir}/${sample}.card.tab

#########################################################
# Additional databases: 		
#abricate --threads $NCPUS --nopath \
#	--db abricate ${contigs} > ${outdir}/${sample}.abricate.tab

#abricate --threads $NCPUS --nopath \
#	--db argannot ${contigs} > ${outdir}/${sample}.argannot.tab

#abricate --threads $NCPUS --nopath \
#	--db ecoh ${contigs} > ${outdir}/${sample}.ecoh.tab

#abricate --threads $NCPUS --nopath \
#	--db plasmidfinder ${contigs} > ${outdir}/${sample}.plasmidfinder.tab

#abricate --threads $NCPUS --nopath \
#	--db vfdb ${contigs} > ${outdir}/${sample}.vfdb.tab

#abricate --threads $NCPUS --nopath \
#	--db ecoli_vf ${contigs} > ${outdir}/${sample}.ecoli_vf.tab
