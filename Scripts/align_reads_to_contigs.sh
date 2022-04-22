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



contigs=${assdir}/${sample}/${sample}.contigs.fa
sampledir=${outdir}/${sample}
mkdir -p $sampledir

#1) Index the assmebled contigs with BWA
bwa index $contigs


#2) Index the assembled contigs with SAMtools
samtools faidx $contigs


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
	
	bwa mem -M -t $NCPUS $contigs \
		-p ${fq[$i]} \
		| samtools sort -n -@ $NCPUS \
		-o ${out}
done

#4) Merge the alignments
bams=$(find ${sampledir} -name "*${sample}*.bam" | xargs echo)
merged=${sampledir}/${sample}.merged.nameSorted.bam
rm -rf $merged

if [ ${#fq[@]} -gt 1 ]
then 
	# Merge multiple BAMs per sample:	
	sambamba merge -t $NCPUS ${merged} ${bams}
else 
	# Rename single BAM per sample for downstream compatibility:
	mv ${fq[0]} ${merged}
fi


#5) Sort and index the merged alignment
final=${sampledir}/${sample}.sort.bam
rm -rf $final

samtools sort -@ $NCPUS \
	-o ${final} ${merged}

samtools index $final

#6) Compute idxstats
samtools idxstats -@ $NCPUS ${final} > ${final}.idxstats




