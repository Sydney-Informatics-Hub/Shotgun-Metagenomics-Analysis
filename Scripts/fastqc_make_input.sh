#! /bin/bash

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

outdir=./FastQC
logdir=./Logs/FastQC
inputs=./Inputs/fastqc.inputs

mkdir -p ${outdir}
mkdir -p ${logdir}

rm -rf ${inputs} ${inputs}-unsorted

for fastq in ./Fastq/*f*q.gz
do
        basename=$(basename "$fastq" | cut -d. -f1)
        out=${outdir}/${basename}
        logfile=${logdir}/${basename}.log
	bytes=$(ls -s "${fastq}" | awk '{print $1}')
        mkdir -p ${out}
        echo "${fastq},${out},${logfile},${bytes}" >> ${inputs}-unsorted
done

# Reverse numeric sort bytes, comma delimited unsorted file
sort -rnk 4 -t ',' ${inputs}-unsorted > ${inputs}
rm -rf ${inputs}-unsorted

num_tasks=`wc -l ${inputs}`
echo "Number of fastQC tasks to run: $num_tasks"
