#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#REQUIRES: samtools, bedtools

HELP(){
	echo "discordant_alignments.sh"
	echo "finds reads that map between an array chromosome and the genome"
	echo
	echo "syntax: discordant_alignments.sh -b BAM --dir WORKING_DIRECTORY [--exclude]"
	echo "NOTE: --exclude currently isn't implemented"
	echo
	echo "-b, --bam	sorted, deduplicated BAM file"
	echo "--dir		path to working directory"
}

#default parameter values

#PARSE OPTIONS
while [ $# -gt 0 ]; do
	case "$1" in
		-b|--bam)
			BAM="$2"
			shift 2	
			;;
		--dir)
			WORKING_DIR="$2"
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

###for each chromosome (that's not the worm chromosomes I-V, X, MtDNA)
###make a bed file merging all the discordant mappings
for chr in $(samtools idxstats ${BAM} | sed \$d | cut -f1); do
        if [ ${chr} != "I" ] && [ ${chr} != "II" ] && [ ${chr} != "III" ] && [ ${chr} != "IV" ] && [ ${chr} != "V" ] && [ ${chr} != "X" ] && [ ${chr} != "MtDNA" ]; then
                echo "finding discordant reads for ${chr}"
				#idenfity discordant reads, make bed file for positions in genome, merge regions 
                samtools view -q 10 ${BAM} ${chr} | bioawk -csam '{if($rnext != "="){print $rnext, $pnext, $pnext+1}}' | bedtools sort | bedtools merge -d 500 > ${WORKING_DIR}/BED_FILES/${chr}.discordant.bed
				#remove regions for each that map to genome (e.g., Psnt-1, tbb-2UTR, etc...)
				bedtools subtract -A -a ${WORKING_DIR}/BED_FILES/${chr}.discordant.bed -b ${WORKING_DIRECTORY}/BED_FILES/${chr}.genomic.bed > ${WORKING_DIR}/BED_FILES/${chr}.discordant.unique.bed
        fi
done

