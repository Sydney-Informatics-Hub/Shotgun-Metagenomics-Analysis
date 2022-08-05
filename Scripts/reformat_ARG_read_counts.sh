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

list=./Inputs/<cohort>_samples.list

while read sample
do
	echo Reformatting sample $sample
	counts=./ARGs/ARG_read_counts/${sample}.curated_ARGs.counts
	cc=$(wc -l < ${counts})
	cc=$(( $cc - 5 )) # the last 5 lines of the htseq output are summaries not ARGs
	
	gc=$(wc -l < ./ARGs/Curated_GFF/${sample}.curated_ARGs.gff)
        
	if [[ $gc -ne $cc ]]
        then
                echo ERROR: $sample GFF has $gc lines counts has $cc lines
		echo Exiting - please fix and resubmit
		exit
        fi	
	
	rf=./ARGs/ARG_read_counts/${sample}.curated_ARGs.reformat.counts
	printf "ID\tLENGTH\tRAW_COUNT\n" > ${rf}
	while read LINE
	do 
		id=$(echo $LINE | cut -d ' ' -f 1)
		count=$(echo $LINE | cut -d ' ' -f 2)
		start=$(echo $id | cut -d ':' -f 3)
		end=$(echo $id | cut -d ':' -f 4) 
		length=$(( $end - $start))
		printf "${id}\t${length}\t${count}\n" >> ${rf}
	done < ${counts}
done < ${list}
