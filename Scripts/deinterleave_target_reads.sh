#!/bin/bash

set -e

sample=`echo $1 | cut -d ' ' -f 1`

indir=./Target_reads
outdir=./Target_reads_paired

fq=$(ls ./${indir}/${sample}*fq.gz)
fq=($fq)

for ((i = 0; i < ${#fq[@]}; i++ ))
do
	file=$(basename ${fq[$i]})
	r1=${file/.interleaved.extracted.fq.gz/_R1_paired.extracted.fq}
	r2=${file/.interleaved.extracted.fq.gz/_R2_paired.extracted.fq}

	reformat.sh \
		in=${fq[$i]} \
		ow=t \
		out1=${outdir}/${r1} \
		out2=${outdir}/${r2}
	
	\rm -rf ${outdir}/${r1}.gz ${outdir}/${r2}.gz
	pigz -p $NCPUS ${outdir}/${r1}
	pigz -p $NCPUS ${outdir}/${r2}

done


