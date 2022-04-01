# Shotgun Metagenomics Analysis
Analysis of metagenomic shotgun sequences including assembly, speciation, ARG discovery and more

## Description
The input for this analysis is paired end next generation sequencing data from metagenomic samples. The workflow is designed to be modular, so that individual modules can be run depending on the nature of the metagenomics project at hand. More modules will be added as we develop them - this repo is a work in progress!

These scripts have been written specifically for NCI Gadi HPC, wich runs PBS Pro, however feel free to use and modify for anothre system if you are not a Gadi user. 

### Part 1. Setup and QC
Download the repo. You will see directories for `Fastq`, `Inputs`, `Reference` and `Logs`. You will need to copy or symlink your fastq to `Fastq`, sample configuration file (see below) to `Inputs` and the reference genome sequence of your host species (if applicable) to `Reference` for host contamination removal.
 

#### Fastq inputs
The scripts assume all fastq files are paired, gzipped, and all in the one directory named 'Fastq'. If your fastq are within a convoluted directory structure (eg per-sample directories) or you would simply like to link them from an alternate location, please use the script `setup_fastq.sh`.

To use this script, parse the path name of your fastq as first argument on the command line, and run the script from the base working directory (<your_path>/Shotgun-Metagenomics-Analysis) which will from here on be referred to as `workdir`. Note that this script looks for `f*q.gz` files (ie fastq.gz or fq.gz) - if yours differ in suffix, please adjust the script accordingly.

```
bash ./Scripts/setup_fastq.sh </path/to/your/parent/fastq/directory>
```

#### Configuration/sample info
The only required input configuration file should be named <cohort>.config, where <cohort> is the name of the current batch of samples you are processing, or some other meaningful name to your project; it will be used to name output files. The config file should be placed inside the $workdir/Inputs directory, and include the following columns, in this order:

```
1. Sample ID - used to identify the sample, eg if you have 3 lanes of sequencing per sample, erach of those 6 fastq files should contain this ID that si in column 1
2. Lab Sample ID - can be the same as column 1, or different if you have reason to change the IDs eg if the seq centre applies an in-house ID. Please make sure IDs are unique within column 1 and unique within column 2
3. Group - eg different time points or treatment groups. If no specific group structure is relevant, please set this to 1 (do not leave blank!) 
3. Platform - should be Illumina; other sequencing platforms are not tested on this workflow
4. Sequencing centre name
5. Library - eg if you have 2 sequencing libraries for the same sample. Can be left blank, or assigned to 1. Blank will be assigned libray ID of 1 during processing.
```

Please do not have spaces in any of the values for the config file. 


#### General setup

All scripts will need to be edited to reflect your NCI project code at the `-P <project>` and `-l <storage> directive. Please run the script create_project.sh and follow the prompts to complete some of the setup for you. 

Note that you will need to manually edit the PDS resource requests for each PBS script; guidelines/example resources will be given at each step to help you do this. As the 'sed' commands within this script operate on .sh and .pbs files, this setup script has been intentionally named .bash (easiest solution).

Remember to submit all scripts from your `workdir`. 

`bash ./Scripts/create_project.sh`

For jobs that execute in parallel, there are 3 scripts: one to make the 'inputs' file listing hte details of each parallel task, one job execution shell script that is run over each task in parallel, and one PBS launcher script. The process is to submit the make input script, check it to make sure your job details are correct, edit the resources directives depending on the number and size of your parallel tasks, then submit the PBS launcher script with `qsub`. 

#### QC

Run fastQC over each fastq file in parallel. Adjust the resources as per your project. To run all files in parallel, set the number of NCPUS requested equal to the number of fastq files (remember that Gadi can only request <1 node or multiples of whole nodes). The make input script sorts the fastq files largest to smallest, so if you have a discrpeancy in file size, optimal efficiency can be achieved by requested less nodes than the total required to run all your fastq in parallel.

FastQC does not multithread on a single file, so CPUs per parallel task is set to 1. Example walltimes on Gadi 'normal' queue:  one 1.8 GB fastq = 4 minutes; one 52 GB fastq file = 69.5 minutes.

Make the fastqc parallel inputs file by running (from `workdir`):
`bash ./Scripts/fastqc_make_inputs.sh`

Edit the resource requests in `fastqc_run_parallel.pbs` according to your number of fastq files and their size, then submit:
`qsub fastqc_run_parallel.pbs`

To ease manual inspection of the fastQC output, running `multiqc` is recommended. This will collate the individual fastQC reports into one report. This can be done on the login node for small sample numbers, or using the below script for larger cohorts. Edit the PBS directives, then run:

`qsub multiqc.pbs`

Save a copy of ./MultiQC/multiqc_report.html to your local disk then open in a web browser to inspect the results. 

#### Quality filtering and trimming

Will be added at a later date. This is highly dependent on the quality of your data and your individual project needs so will be a guide only. 

### Part 2. Removal of host contamination. 

If you have metagenomic data extracted from a host, you will need a copy of the host reference genome sequence in order to remove any DNA sequences belonging to the host. Even if your wetlab protocol included a host removal step, it is still important to run bioinformatic host removal.


#### Prepare the reference
Ensure you have a copy of the reference genome (or symlink) in ./Fasta. This workflow requires BBtools(tested with version 37.98). As of writing, BBtools is not available as a global app on Gadi. Please install locally and make "module loadable", or else edit the scripts to point directly to your local BBtools installation.

BBtools repeat masking will use all available threads on machine and 85% of available mem by default. For a mammalian genome, 2 hours on one Gadi 'normal' node is sufficient for repeat masking. 

Update the name of your reference fastq in the `bbmap_prep.pbs` script (and BBtools, see note above), then run:
`qsub ./Scripts/bbmap_prep.pbs`

#### Host contamination removal

TBC 1/4/22... 
