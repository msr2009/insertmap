#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#requires bowtie2, samtools

HELP(){
	echo "align_to_genome_bt2.sh"
	echo "wrapper script to align array reads to a genome"
	echo "creates alignment genome+array index, then aligns"
	echo "does _not_ perform post-alignment processsing"
	echo
	echo "syntax: align_to_genome.sh -g GENOME -x PREFIX -1 READ1 -2 READ2 --array_dir"
	echo "-g, --genome	location of genome FASTA"
	echo "-x, --prefix	prefix for naming output (should contain full path)"
	echo "-1	location of forward read FASTQ"
	echo "-2	location of reverse read FASTQ"
	echo "--array_dir	directory with array sequences"
	echo "-t, --threads	number of threads to use for alignment"
	echo "-h, --help	show this help"
}

#default parameter values
THREADS=1

#PARSE OPTIONS
while [ $# -gt 0 ]; do
	case "$1" in
		-g|--genome)
			GENOME="$2"
			shift 2
			;;
		-x|--prefix)
			PREFIX="$2"
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
		-t|--threads)
			THREADS="$2"
			shift 2
			;;
		-d|--array_dir)
			ARRAY_SEQS_DIR="$2"
			shift 2
			;;
		-h|--help) #HELP!
			HELP
			exit;;
		*)
			HELP
			exit;;
	esac
done

#create bowtie2 index
IDX=${GENOME}
for f in ${ARRAY_SEQS_DIR}/*.fa; do
	   IDX=${IDX},${f}	
done

bowtie2-build --threads ${THREADS} $IDX ${PREFIX}

###align reads and convert to bam
bowtie2 -t -k 5 -p ${THREADS} -x ${PREFIX} -1 ${READ1} -2 ${READ2} | samtools view -bS > ${PREFIX}.bam

