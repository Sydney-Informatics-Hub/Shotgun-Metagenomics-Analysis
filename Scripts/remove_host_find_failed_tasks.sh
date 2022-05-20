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


grep "exited with status 0" ./Logs/remove_host.e | awk '{print $11}' | sort  > ./Inputs/remove_host_passed

sort ./Inputs/remove_host.inputs > ./Inputs/remove_host.inputs-sorted

grep -v -f  ./Inputs/remove_host_passed ./Inputs/remove_host.inputs-sorted > ./Inputs/remove_host.inputs-failed

\rm -rf ./Inputs/remove_host.inputs-sorted ./Inputs/remove_host_passed

rerun=$(wc -l < ./Inputs/remove_host.inputs-failed )

echo There are $rerun remove host tasks to rerun
echo These have been written to ./Inputs/remove_host.inputs-failed
echo Please edit the resource requests and run ./Scripts/remove_host_failed_run_parallel.pbs
