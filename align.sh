#!/bin/bash
#PBS -S /bin/sh
#PBS -o logs/
#PBS -M ammar.aziz@menzies.edu.au
#PBS -m n
#PBS -j oe
#PBS -l ncpus=4
#PBS -l walltime=96:00:00

cd $wkdir
echo $wkdir
echo "reference: $ref"

# Align to ref
if [ ! -s "${seq}_SDR.bam" ]; then
	echo "Aligning to reference................................."
	$BWA mem -t 2 -R '@RG\tID:'$org'\tSM:'$seq'\tPL:ILLUMINA\' $ref ${seq}_1_sequence.fastq.gz ${seq}_2_sequence.fastq.gz > ${seq}.sam
	echo "Convert to bam, sort and index................................."
	$SAMTOOLS view -@ 2 -h -b -q 1 ${seq}.sam > ${seq}.bam 
	$SAMTOOLS sort -@ 2 -o ${seq}_sorted.bam ${seq}.bam 
	echo "Removing duplicates................................."
	$JAVA $SET_VAR $GATK MarkDuplicates -I ${seq}_sorted.bam -O ${seq}_SDR.bam --REMOVE_DUPLICATES true -M $wkdir/logs/${seq}_dups_metrics.txt --VALIDATION_STRINGENCY LENIENT
	$SAMTOOLS index ${seq}_SDR.bam
fi

#clean up
echo removing ${seq}.sam ${seq}.bam  ${seq}_sorted.bam
rm ${seq}.sam ${seq}.bam  ${seq}_sorted.bam

# variant calling using HaplotypeCaller
if [ ! -s "${seq}_raw.g.vcf" ]; then
	echo "calling variants using haplotypeCaller in gvcf mode................................."
	$JAVA $SET_VAR $GATK HaplotypeCaller -ERC GVCF -R $ref -ploidy $ploidy -I $wkdir/${seq}_SDR.bam -O $wkdir/variants/${seq/_SDR.bam/}_raw.g.vcf 
fi

sleep 2
exit 0