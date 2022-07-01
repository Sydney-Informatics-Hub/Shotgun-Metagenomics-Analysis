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

dt=`date '+%d-%m-%Y'`

dbname=kraken2_standard_db_build_${dt}

printf "Using kraken2 database name: ${dbname}\n"
echo

printf "Submitting NCBI download for bacterial library (8-9 hours run time)\n" 
qsub -N k2-bacteria -l walltime=10:00:00 -l mem=60GB -o ./Logs/kraken2_download_bacteria.o  -e ./Logs/kraken2_download_bacteria.e -v library="bacteria",dbname="$dbname" ./Scripts/kraken2_download.pbs
echo
sleep 1

printf "Submitting NCBI download for taxonomy and human, viral, archaea, and Univec_Core libraries (< 1 hour run time)\n" 
qsub -N k2-download -l walltime=02:00:00 -l mem=12GB -o ./Logs/kraken2_download_other.o  -e ./Logs/kraken2_download_other.e -v library="archaea viral human UniVec_Core",dbname="$dbname" ./Scripts/kraken2_download.pbs
echo

printf "Once all libaries are downloaded, please submit ./Scripts/kraken2_build.pbs\n"


