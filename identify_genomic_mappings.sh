#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#requires bowtie2, bedtools, blast, bioawk

#display help
HELP(){
	echo "identify_genomic_mappings.sh"
	echo "performs initial mapping of array sequences to genome"
	echo "and compiles mappings into BED file."
	echo
	echo "syntax: identify_genomics_mappings.sh -g GENOME_FASTA -s ARRAY_FASTA"
	echo "options:"
	echo "-g, --genome	location of genome FASTA file"
	echo "-s, --seq	array sequence FASTA file"
	echo "-e	BLAST e-value (default 1e-70)"
}

#PARSE OPTIONS
while [$# -gt 0]; do
	case "$1" in
		-g|--genome)
			GENOME="$2"
			shift 2
			;;
		-s|--seq)
			SEQ="$2"
			shift 2
			;;
		-e)
			EVALUE="$2"
			shift 2
			;;
		-h|--help) #HELP!
			HELP
			exit;;
	esac
done

seqname=`basename ${SEQ} .fa`

echo MAPPING ${seqname} against genome

blastn -db ${GENOME} -query ${SEQ} -evalue 1e-70 -outfmt 6 | \
bioawk -t 'if($9>$10){print $2, $10, $9} else{print $2, $9, $10}}' | \
bedtools sort | bedtools merge | bioawk -t -v sn=${seqname} '{print $0, sn}'
