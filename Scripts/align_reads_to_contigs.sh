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

sample=`echo $1 | cut -d ',' -f 1` 
platform=`echo $1 | cut -d ',' -f 2`
centre=`echo $1 | cut -d ',' -f 3`
lib=`echo $1 | cut -d ',' -f 4` 

assdir=./Assembly
fastqdir=./Target_reads
outdir=./Align_to_assembly
log=./Logs/Align_to_assembly/${sample}.log
rm -rf $log

contigs=${assdir}/${sample}/${sample}.contigs.fa
sampledir=${outdir}/${sample}
mkdir -p $sampledir

#1) 
dt=$(date)
printf "${dt}: Indexing the assmebled contigs with BWA\n" >> $log
bwa index $contigs 2>>${log}
echo

#2) 
dt=$(date)
printf "${dt}: Indexing the assembled contigs with SAMtools\n" >> $log
samtools faidx $contigs 2>>${log}
echo

#3) Align the extracted target reads to the assembled reads
# For samples with multiple pairs of fastq, need to map pairs separately then merge with sambamba
fq=$(ls ${fastqdir}/*${sample}*) #interleaved
fq=($fq)

for (( i = 0; i < ${#fq[@]}; i++ ))
do
	fqname=${fq[$i]}
	fqname=$(basename ${fq[$i]})
	fqname=$(echo $fqname | cut -d '.' -f 1)
	out=${sampledir}/${fqname}.bam
	rm -rf ${out}*
	
	# Get the metadata for the BAM. Sequencing centre and platform are derived from config. Lane and flowcell are derived from the fastq. 
	# Assumes standard illumina read ID format (flowcell field 3 and lane field 4 of : delim string)
	flowcell=$(zcat ${fq[$i]} | head -1 | cut -d ':' -f 3) # standard illumina format
	lane=$(zcat ${fq[$i]} | head -1 | cut -d ':' -f 4)	
	
	dt=$(date)
	printf "${dt}: Aligning ${fq[$i]} to ${contigs}, output is ${out}\n"
	bwa mem -M -t $NCPUS $contigs \
		-R "@RG\tID:${sample}_${lib}\tPL:${platform}\tSM:${sample}\tLB:${sample}_${lib}\tCN:${centre}" \
		-p ${fq[$i]} 2>>${log} \
		| samtools sort -n -@ $NCPUS \
		-o ${out} - 2>>${log}
	echo
done

#4) Merge the alignments
merged=${sampledir}/${sample}.merged.nameSorted.bam
rm -rf $merged
bams=$(find ${sampledir} -name "*${sample}*.bam" | xargs echo)

if [ ${#fq[@]} -gt 1 ]
then 
	# Merge multiple BAMs per sample:
	dt=$(date)
	printf "${dt}: Merging ${bams} to ${merged}\n"	
	sambamba merge -t $NCPUS ${merged} ${bams} 2>>${log}
else 
	# Rename single BAM per sample for downstream compatibility:
	dt=$(date)
	printf "${dt}: Rename ${bams[0]} to ${merged} for downstream compatibility\n"
	mv ${bams[0]} ${merged}
fi
echo

#5) 
dt=$(date)
printf "${dt}: Sorting and indexing ${merged} to ${final}\n"
final=${sampledir}/${sample}.sort.bam
rm -rf $final

samtools sort -@ $NCPUS \
	-o ${final} ${merged} 2>>${log}

samtools index $final 2>>${log}

#6) 
dt=$(date)
printf "${dt}: Computing idxstats\n"
samtools idxstats -@ $NCPUS ${final} > ${final}.idxstats 2>>${log}




