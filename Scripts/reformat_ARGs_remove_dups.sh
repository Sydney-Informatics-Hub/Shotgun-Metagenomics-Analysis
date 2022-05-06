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
all_out=./ARGs/Abricate/<cohort>.ARGs.txt
rm -rf ${all_out}

while read SAMPLE
do
	echo Reformatting $SAMPLE
	
	dir=./ARGs/Abricate/${SAMPLE}
	final=${dir}/${SAMPLE}.ARGs.txt
	
	# Combine all database output:
	cat $dir/${SAMPLE}*.tab > ${final} 
	
	# Remove extra headers:
	sed -i '1!{/^\#/d;}' ${final}
	
	# Remove exact duplicates:	
	no_dups=${dir}/${SAMPLE}.ARGs.txt-no-dups
	head -1 ${final} > ${no_dups}
        lines=$(wc -l < $final)
        ((lines--))
        tail -n ${lines} ${final} | sort | uniq >> ${no_dups}
        mv ${no_dups} ${final}
	
	# Add to all-out file:
	cat ${final} >> $all_out
		
done < $list

# Remove extra headers:
sed -i '1!{/^\#/d;}' ${all_out}
