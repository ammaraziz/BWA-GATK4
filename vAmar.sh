#!/bin/bash
#PBS -S /bin/sh
#PBS -M ammar.aziz@menzies.edu.au
#PBS -m n
#PBS -j oe
#PBS -o logs/
#PBS -l ncpus=2
#PBS -l walltime=96:00:00
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cat $DIR/ascii.txt
echo -e "\n"
sleep 2
##############################################################
##############################################################
##############################################################
USAGE="Usage: vAmar.sh -r reference.fasta"

if [ "$#" == "0" ]; then                      # If zero arguments were supplied,
  echo "Error: no reference provided."
  echo "$USAGE"                               # display a help message
  exit 1                                      # and return an error.
fi
# read the options
# r: = mandatory
# p:: = optional 
TEMP=`getopt -o r:p: --long reference:ploidy:: -n 'vAmar.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
	#####
        -r|--reference)
            case "$2" in
                *) ref=$2 ; shift 2 ;;
            esac ;;
        -p|--ploidy)
            case "$2" in
                *) ploidy=$2 ; shift 2 ;;
            esac ;;
        -a|--arga)
            case "$2" in
                "") ARG_A='some default value' ; shift 2 ;;
                *) ARG_A=$2 ; shift 2 ;;
            esac ;;
	#####		
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done
##############################################################
##############################################################
##############################################################
echo -e "You have selected $ref as the reference."
echo -e "Finding reads.....\r"
echo -e "Ploidy set to: $ploidy"
sleep 1


# Program Variables
BWA="/home/aaziz/analysis/bp/mixture/gatk4/bwa"
SAMTOOLS="/home/aaziz/bin/miniconda3/bin/samtools"
GATK="/home/aaziz/bin/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"
JAVA="/home/aaziz/bin/miniconda3/bin/java"
MELT="$DIR/vcf_melt.py"
# Universal variables
cpus=2
#SET_VAR="-jar -XX:+UseSerialGC -Xmx4G"
SET_VAR="-jar -XX:ParallelGCThreads=4 -Xmx8G"
wkdir=$PWD
WALL_T=96:00:00
org=haploid

# Creates a list of .fastq files
sequences_tmp=($(find $PWD/*_1_sequence.fastq.gz -printf "%f "))
sequences=("${sequences_tmp[@]/_[12]_sequence.fastq.gz/}")
n=${#sequences[@]}
echo -e "detected samples: ${sequences[@]} \n"
sleep 1
# Dependancy function for alignment
depend_align() {
  ids=`cat align_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
  depend="-W depend=afterok${ids}"
  echo $depend
}

depend_ids() {
	ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
	depend="-W depend=afterok${ids}"
	echo $depend
}

#cleanup/creation
if [ ! -d "logs" ]; then
	mkdir logs
fi

if [ ! -d "variants" ]; then
	mkdir variants
fi

if [ -s align_ids.txt ]; then
    rm align_ids.txt
fi

if [ -s qsub_ids.txt ]; then
    rm qsub_ids.txt
fi


#reference: 
#	index with samtools and BWA
#	dictionaries with GATK and picard

#bwa index

	if [ ! -s "$ref.bwt" ]; then
	echo -e "BWA reference indexing"
	cmd="$BWA index -a is -p $ref $wkdir/$ref"
	qsub_id=`qsub -N bwa_index -j oe -m n -l ncpus=$cpus,walltime=$WALL_T -v command="$cmd" "$DIR"/Header.pbs`
	echo -e "index\t$qsub_id" >> qsub_ids.txt
	fi
	
#samtools index
	if [ ! -s "$ref.fai" ]; then
	echo -e "SAMtools reference indexing"
	cmd="$SAMTOOLS faidx $wkdir/$ref"
	qsub_id=`qsub -N sam_index -j oe -m n -M ammar.aziz@menzies.edu.au -l ncpus=$cpus,walltime=$WALL_T -v command="$cmd" "$DIR"/Header.pbs`
	echo -e "SAM_index\t$qsub_id" >> qsub_ids.txt
	fi

#picard dictionary
	if [ ! -s "${ref/.fasta/}.dict" ]; then
		echo -e "Picard dictionary creation"
		cmd="$JAVA $SET_VAR $GATK CreateSequenceDictionary -R $wkdir/$ref -O ${ref/.fasta/}.dict"
		qsub_id=`qsub -N dictionary -j oe -m n -M ammar.aziz@menzies.edu.au -l ncpus=$cpus,walltime=$WALL_T -v command="$cmd" "$DIR"/Header.pbs`
		echo -e "PICARD_dict\t$qsub_id" >> qsub_ids.txt

	fi

#if [ -s qsub_ids.txt ]; then
	depend_idx=$(depend_ids)
#fi

for (( i=0; i<n; i++ )); do
	if [ ! -s variants/${sequences[i]}_raw.g.vcf ]; then
		echo "qsub: Alignment and variant calling of ${sequences[i]}"
		var="seq=${sequences[i]}, ref=$ref, SAMTOOLS=$SAMTOOLS, BWA=$BWA, org=$org, wkdir=$PWD, GATK=$GATK, JAVA=$JAVA, SET_VAR=$SET_VAR, ploidy=$ploidy"
   		align_ids=`qsub -N ${sequences[i]}.align -j oe -m n -M ammar.aziz@menzies.edu.au -l ncpus=$cpus,walltime=$WALL_T "$depend_idx" -v "$var" "$DIR"/align.sh`
		echo -e "align_${sequences[$i]}\t$align_ids" >> align_ids.txt;
		if [ -s align_ids.txt ]; then
			depend_align_ids=$(depend_align)
		fi
	fi
done

	if [ ! -s "variants/cohort_jcf.g.vcf" ]; then
		echo "qsub: Combined variant calling and filtering"
		var="ref=$ref, wkdir=$PWD, GATK=$GATK, JAVA=$JAVA, SET_VAR=$SET_VAR, MELT=$MELT"
   		#combine_filter=`qsub -N combine_filter -j oe -m n -M ammar.aziz@menzies.edu.au -l ncpus=$cpus,walltime=$WALL_T "$depend_align_ids" -v "$var" "$DIR"/combine_filter.sh`
		combine_filter=`qsub -N combine_filter -j oe -m n -M ammar.aziz@menzies.edu.au -l ncpus=$cpus,walltime=$WALL_T $depend_align_ids -v "$var" "$DIR"/combine_filter.sh`
	fi

exit 0
