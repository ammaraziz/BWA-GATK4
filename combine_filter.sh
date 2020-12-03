#!/bin/bash
#PBS -S /bin/sh
#PBS -o logs/
#PBS -M ammar.aziz@menzies.edu.au
#PBS -m n
#PBS -j oe
#PBS -l ncpus=2
#PBS -l walltime=96:00:00

# generic variables
cd $wkdir
#needed to melt VCF files
which python
python -V


#GATK Variables
CLUSTER_SNP=3
CLUSTER_WINDOW_SNP=10
#MLEAF_SNP=0.95
QD_SNP=10.0
MQ_SNP=30.0
FS_SNP=10.0
HAPLO_SNP=20.0
QUAL_SNP=30.0
LOW_DEPTH=2 # A value of 2 means that regions with less than half of the average coverage of the entire genome will fail
HIGH_DEPTH=3 # A value of 5 means that regions with 5x depth greater than the average coverage will fail

# combine gvcf files
if [ ! -s "$wkdir/variants/cohort.g.vcf" ]; then
	echo -e "Combining all gvcf files using CombineGVCFs............................................................................................"
	find $PWD/variants/ -name "*_raw.g.vcf" > $PWD/variants/input.list
	$JAVA $SET_VAR $GATK CombineGVCFs -R $ref --variant $PWD/variants/input.list -O $PWD/variants/cohort.g.vcf
fi

# joint genotyping using GenotypeGVCFs
if [ ! -s "$wkdir/variants/cohort_jc.g.vcf" ]; then
	echo "joint calling using GenotypeGVCFs............................................................................................"
	$JAVA $SET_VAR $GATK GenotypeGVCFs --annotations-to-exclude InbreedingCoeff -R $ref -V $PWD/variants/cohort.g.vcf -O $PWD/variants/cohort_jc.g.vcf
fi

# Variant filtering 
if [ ! -s "variants/cohort_jcf.g.vcf" ]; then
	echo "filtering............................................................................................"
	echo $JAVA $SET_VAR $GATK
	$JAVA $SET_VAR $GATK VariantFiltration -V variants/cohort_jc.g.vcf -O variants/cohort_jcf.g.vcf --cluster-size $CLUSTER_SNP --cluster-window-size $CLUSTER_WINDOW_SNP --filter-expression "QD < $QD_SNP" --filter-name "QDFilter" --filter-expression "MQ < $MQ_SNP" --filter-name "MQFilter" --filter-expression "FS > $FS_SNP" --filter-name "FSFilter" --filter-expression "QUAL < $QUAL_SNP" --filter-name "StandardFilters"
fi

# subset variants, keep only SNPs and INDEls
# melt for more R friendly format.
if [ ! -s "variants/filtered_SNPs_clean.g.vcf" -a  ! -s "variants/filtered_INDELs_clean.g.vcf" ]; then
	echo "selecting variants............................................................................................"
	$JAVA $SET_VAR $GATK SelectVariants -R $ref -V variants/cohort_jcf.g.vcf -O variants/filtered_SNPs.g.vcf --exclude-filtered --select-type-to-include SNP
	$MELT variants/filtered_SNPs.g.vcf > variants/filtered_SNPs_melt.g.vcf
	$JAVA $SET_VAR $GATK SelectVariants -R $ref -V variants/cohort_jcf.g.vcf -O variants/filtered_INDELs.g.vcf --exclude-filtered --select-type-to-include INDEL
	$MELT variants/filtered_INDELs.g.vcf > variants/filtered_INDELs_melt.g.vcf
fi

sleep 2
exit 0