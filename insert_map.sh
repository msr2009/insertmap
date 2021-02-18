#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#requires bowtie2, samtools, blast, bioawk

#a loop to get arguments
#while getopts t:x:1:2:d: flag
#do
#    case "${flag}" in
#        t) threads=${OPTARG};;
#		x) _name=${OPTARG};;
#		1) read1=${OPTARG};;
#		2) read2=${OPTARG};;	
#		d) WORKING_DIR=${OPTARG};;
#	esac
#done

# Default values for blank parameters
INIT=0
CLEAN=0
GENOME=
PREFIX=
THREADS=1
EVALUE=1e-70

#display help
HELP(){
	echo "PLACEHOLDER FOR PROGRAM DESCRIPTION"
	echo
	echo "FOR MAPPING"
	echo "syntax: -d WORKING_DIRECTORY -1 READ1 -2 READ2 -g GENOME_FASTA -x PREFIX"
	echo "required options:"
	echo "-d, --dir	working directory"
	echo "-x, --prefix	prefix for files (index, results, etc)"
	echo "-g, --genome	location of genome FASTA file"
	echo "-1		location of forward read FASTQ"
	echo "-2		location of reverse read FASTQ"
	echo
	echo "optional parameters:"
	echo "-t, --threads	number of threads to use for alignment, samtools"
	echo "-e		BLAST e-value (default=1e-70)"
	echo "-h, --help	show this help"
	echo 
	echo "MAINTENANCE OPTIONS"
	echo "--init		inititalize directory structure (requires -d)"
	echo "--clean		clean working directory (requires -d)"	
}

if [ $# -eq 0 ]
then
	HELP
	exit
fi	

#PARSE OPTIONS
while [ $# -gt 0 ]; do
    case "$1" in
		-d|--dir)
			WORKING_DIR="$2"
			shift 2
			;;
		-1)
            READ1="$2"
            shift 2
            ;;
        -2)
            READ2="$2"
            shift 2
            ;;
		-x|--prefix)
	    	PREFIX="$2"
	    	shift 2
	    	;;
		-t|--threads)
			THREADS="$2"
			shift 2
			;;
		-g|--genome)
			GENOME="$2"
			shift 2
			;;
		-e)
			EVALUE="$2"
			shift 2
			;;
        --init) 
            INIT=1
            shift 1
            ;;
		--clean) 
            CLEAN=1
            shift 1
            ;;
		-h|--help) #HELP ME PLEASE!
			HELP
			exit;;
		*)
			HELP
            exit;;
    esac
done

ARRAY_SEQS_DIR=${WORKING_DIR}/ARRAY_SEQS
BED_DIR=${WORKING_DIR}/BED_FILES
_name=${WORKING_DIR}/${PREFIX}

###MAINTENANCE STUFF

#initialize directories
if (( ${INIT} == 1 )); then
	if [ -d ${WORKING_DIR} ]; then
		echo "${WORKING_DIR} exists."
	else		
		mkdir ${WORKING_DIR}
	fi
	echo "creating directory structure at ${WORKING_DIR}"
	mkdir ${ARRAY_SEQS_DIR}
	mkdir ${BED_DIR}
	mkdir ${WORKING_DIR}/READS/
	exit 1
fi

#clean directories
#so, rm .bed's, .bt2's, .bam's
if (( ${CLEAN} == 1 )); then
	echo "cleaning ${WORKING_DIR}"
	echo "removing BED, BAM, BT2"
	rm ${WORKING_DIR}/*.bed
	rm ${WORKING_DIR}/*.bam*
	rm ${WORKING_DIR}/*.bt2
	rm ${ARRAY_SEQS_DIR}/*.bed
	rm ${BED_DIR}/*
	exit 1
fi

#check that required parameters are set
if [ -z $WORKING_DIR ]
then
	echo "ERROR: must set working directory path with -d"
	echo
	HELP
	exit
fi	

if [ -z $GENOME ]
then
	echo "ERROR: must set path to genome with -g"
	echo
	HELP
	exit
fi	

if [ -z $PREFIX ]
then
	echo "ERROR: must set prefix with -x"
	echo 
	HELP
	exit
fi	

if [ -z $READ1 ]
then
	echo "ERROR: must set read1 path with -1"
	echo
	HELP
	exit
fi

if [ -z $READ2 ]
then
	echo "ERROR: must set read2 path with -2"
	echo
	HELP
	exit
fi

##############
###MAIN SCRIPT	
#############
echo
echo "######################################"
echo NAME: ${PREFIX}
echo GENOME: ${GENOME}
echo READ1: ${READ1}
echo READ2: ${READ2}
echo WORKING DIRECTORY: ${WORKING_DIR}
echo ARRAY SEQS DIRECTORY: ${ARRAY_SEQS_DIR}
echo "######################################"

echo 
echo "######################################"
echo "CONVERTING ARRAY SEQUENCE FILES"
echo "######################################"

#first convert files that need to be converted to FASTA
for seq in ${ARRAY_SEQS_DIR}/*; do
	python convert_to_fasta.py ${seq}
done

echo
echo "######################################"
echo "IDENTIFY ARRAY SEQ MAPPINGS IN GENOME"
echo "######################################"

### blast against genome db and convert top hits to bed file
### blastn finds all "good" matches to genome
### bioawk converts blast output to bed file (and orders coords appropriately)
### bedtools merges overlapping hits (and sorts)
for seq in ${ARRAY_SEQS_DIR}/*.fa; do	
	seqname=`basename ${seq} .fa`
	sh identify_genomic_mappings.sh -g ${GENOME} -s ${seq} -e ${EVALUE} > ${BED_DIR}/${seqname}.genomic.bed
done

#merge all bed files
#cat ${BED_DIR}/*.genomic.bed > ${BED_DIR}/tmp_all.bed
#bedtools sort tmp_all.bed | bedtools merge - > ARRAY.genomic.bed
#rm tmp_all.bed

echo
echo "######################################"
echo "ALIGNING TO READS TO GENOME"
echo "######################################"

sh align_to_genome_bt2.sh -g ${GENOME} -x ${_name} -t ${THREADS} -1 ${READ1} -2 ${READ2} -d ${ARRAY_SEQS_DIR}

### create bowtie2-index
#IDX=${GENOME}
#for f in ${ARRAY_SEQS_DIR}/*.fa; do 
#		IDX=${IDX},${f}
#done

#bowtie2-build $IDX ${_name}

#echo -e "\n"
#echo "######################################"
#echo "ALIGNING READS"
#echo "######################################"
#echo -e "\n"

###align reads and convert to bam
### -a finds all alignments, rather than the default (greedy) approach
#bowtie2 -t -k 5 -p ${THREADS} -x ${_name} -1 ${READ1} -2 ${READ2} | samtools view -bS > ${_name}.bam

echo 
echo "######################################"
echo "PROCESSING ALIGNED READS"
echo "######################################"

sh process_alignment.sh ${_name}.bam

###sort, rmdup, and index bam
#samtools collate -O ${_name}.bam | samtools fixmate -m - - | samtools sort - | samtools markdup -s - ${_name}.srt.rmdup.bam
#samtools sort ${_name}.bam | samtools rmdup - ${_name}.srt.rmdup.bam
#samtools index ${_name}.srt.rmdup.bam ${_name}.srt.rmdup.bam.bai

echo
echo "######################################"
echo "IDENTIFYING DISCORDANT READS"
echo "######################################"

sh discordant_alignments.sh -b ${_name}.srt.rmdup.bam --dir ${WORKING_DIR}

###for each chromosome (that's not the worm chromosomes I-V, X, MtDNA)
###make a bed file merging all the discordant mappings
#for chr in $(samtools idxstats ${_name}.srt.rmdup.bam | sed \$d | cut -f1); do
#	if [ ${chr} != "I" ] && [ ${chr} != "II" ] && [ ${chr} != "III" ] && [ ${chr} != "IV" ] && [ ${chr} != "V" ] && [ ${chr} != "X" ] && [ ${chr} != "MtDNA" ]; then
#        	echo "finding discordant reads for ${chr}"
#			#idenfity discordant reads, make bed file for positions in genome, merge regions 
#        	samtools view -q 10 ${_name}.srt.rmdup.bam ${chr} | bioawk -csam '{if($rnext != "="){print $rnext, $pnext, $pnext+1}}' | bedtools sort | bedtools merge -d 500 > ${WORKING_DIR}/${chr}.discordant.bed
#			#remove regions for each that map to genome (e.g., Psnt-1, tbb-2UTR, etc...)
#			bedtools subtract -A -a ${WORKING_DIR}/${chr}.discordant.bed -b ${ARRAY_SEQS_DIR}/ARRAY.blast.bed > ${WORKING_DIR}/${chr}.discordant.unique.bed
#	fi
#done
