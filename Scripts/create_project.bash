#!/bin/bash

#########################################################
#
# Platform: NCI Gadi HPC
# Description: see https://github.com/Sydney-Informatics-Hub/Shotgun-Metagenomics-Analysis
#
# Author/s: Cali Willet; Tracy Chew
# cali.willet@sydney.edu.au; tracy.chew@sydney.edu.au
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


#### functions####
function storage {
echo Do you require read/write to any Gadi storage other than /scratch/${project}? If yes, please enter all separated by space \[enter for no\]:

read more_storage
IFS=' ' read -r -a array <<< "$more_storage"

lstorage=scratch/${project}
for i in "${!array[@]}"
do
    path=$(echo ${array[i]} | sed 's/^\///g')
    lstorage+="+scratch/${path}"
done

echo
echo PBS lstorage directive will be: $lstorage
echo Is this correct? Enter y or n

read answer

if [[ $answer != y ]]
then 
	storage
else
	echo Using lstorage $lstorage
	echo
	return 0	
fi

}
###################

# Config file
echo Enter your cohort name / the basename of your <cohort>.config file:
read cohort
config=./Inputs/${cohort}.config

if [ -f $config ] 
then
	echo Using config $config
	sed -i "s/^cohort=.*/cohort=${cohort}/" ./Scripts/*.sh
	sed -i "s/^cohort=.*/cohort=${cohort}/" ./Scripts/*.pbs 	
	echo
else
	echo $config does not exist - please fix. Aborting.
	exit
fi

# NCI project
echo Enter the name of your NCI project:
read project

echo Using NCI project $project for accounting and /scratch/${project} for read/write
sed -i "s/#PBS -P.*/#PBS -P ${project}/" ./Scripts/*.pbs
echo


# Call storage function as many times as needed
storage
sed -i "s|#PBS -lstorage=.*|#PBS -lstorage=${lstorage}|" ./Scripts/*.pbs


###################

echo The scripts in ./Scripts directory have now been updated to include the following:
printf "\tNCI accounting project: ${project}\n \
\tPBS lstorage directive: ${lstorage}\n \
\tCohort config file: ./Inputs/${cohort}.config\n \







