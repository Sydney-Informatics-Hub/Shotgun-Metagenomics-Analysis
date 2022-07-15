#!/bin/bash

# BLAST Sample CDS contigs (Prodigal) to
# NCBI NR protein DB using DIAMOND
# Due to extremley variable walltimes (that cannot be correlated
# with any attribute of the data eg size) not using parallel mode

set -e

#########################################################
#
# Platform: NCI Gadi HPC
#
# Author: Tracy Chew
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

samples=$(awk '{print $1}' ./Inputs/samples.list)
samples=(${samples})
logdir=./Logs/Diamond_NCBI

mkdir -p ${logdir}
mkdir -p ../Diamond_NCBI

for (( i = 0; i < ${#samples[@]}; i++ ))
do
        sample=${samples[$i]}
        echo Submitting diamond for $sample
        qsub -N d-${sample} -o ${logdir}/${sample}.o -e ${logdir}/${sample}.e -v sample="$sample" ./diamond_taxon.pbs
        sleep 1
        echo
done
