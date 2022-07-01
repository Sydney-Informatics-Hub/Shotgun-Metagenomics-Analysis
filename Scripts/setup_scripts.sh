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


#### functions####
function storage {
echo Do you require read/write to any Gadi storage other than /scratch/${project}? If yes, please enter all paths separated by space \[enter for no\]:

read more_storage
IFS=' ' read -r -a array <<< "$more_storage"

lstorage=scratch/${project}
for i in "${!array[@]}"
do
    path=$(echo ${array[i]} | sed 's/^\///g')
    lstorage+="+${path}"
done

echo
printf "PBS -l storage directive will be: $lstorage\n"
echo Is this correct? Enter y or n

read answer

if [[ $answer != y ]]
then 
	storage
else
	echo Using storage $lstorage
	echo
	return 0	
fi

}
###################

# Config file
printf "Enter your cohort name / the basename of your <cohort>.config file:\n"
read cohort
config=./Inputs/${cohort}.config

if [ -f $config ] 
then
	echo Using config $config
	sed -i "s/<cohort>/${cohort}/" ./Scripts/*.sh ./Scripts/*.pbs ./Scripts/*.pl
	echo
else
	echo $config does not exist - please fix. Aborting.
	exit
fi

echo Making sample list ./Inputs/${cohort}_samples.list
awk 'NR>1 {print $2}' ${config} > ./Inputs/${cohort}_samples.list
echo

# NCI project
echo Enter the name of your NCI project:
read project

echo Using NCI project $project for accounting and /scratch/${project} for read/write
sed -i "s/#PBS -P.*/#PBS -P ${project}/" ./Scripts/*.pbs
echo

###################

# Call storage function as many times as needed
storage
sed -i "s|#PBS -l storage=.*|#PBS -l storage=${lstorage}|" ./Scripts/*.pbs

###################

# Reference genome

printf "Is your metagenomics data extracted from a host? (eg tissue, saliva) - enter 'y' or 'n'\n"
read host

if [ $host == "y" ] || [ $host == "Y" ] 
then 
	printf  "Host-extracted data requires a host removal step against a reference genome.\nEnter the name of your host reference genome sequence fasta file, including full path:\n"
	read refpath
	
	if [ -f $refpath ] 
	then
		echo Using host reference genome sequence $refpath
		sed -i "s|<reference>|${refpath}|" ./Scripts/bbmap_prep.pbs
		echo
	else
		echo $refpath does not exist - please fix. Aborting.
		exit
	fi
else 
	echo No host reference genome sequence specified.
	echo
fi


###################

echo The scripts in ./Scripts directory have now been updated to include the following:
printf "\tNCI accounting project: ${project}\n \
\tPBS lstorage directive: ${lstorage}\n \
\tCohort config file: ./Inputs/${cohort}.config\n"
if [ -f $refpath ]
then
	printf "\tHost reference genome sequence: ${refpath}\n"
fi
echo





