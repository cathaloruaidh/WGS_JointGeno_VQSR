#!/bin/bash

############################################################
#
# Assembly of NGS data from FASTQ to GVCF 
# Cathal Ormond 2017. 
#
# Alignment of reads from FASTQ file to reference genome
# Conversion to BAM file
#
# Based heavily on scripts found here:
# https://github.com/genepi-freiburg/gwas
#
############################################################




### Call the variants
# Run the picard command


#read  -n 1 -p "Input Selection:" mainmenuinput
log "Apply VQSR to indels at tranche filter level ${INDEL_TRANCHE_LEVEL}" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T ApplyRecalibration \
	-R ${REF_FASTA} \
	-input ${RESULTS_DIR}/${INDEL_VCF_GZ} \
	-mode INDEL \
	--ts_filter_level ${INDEL_TRANCHE_LEVEL} \
	-recalFile ${RESULTS_DIR}/${INDEL_RECAL_TXT} \
	-tranchesFile ${RESULTS_DIR}/${INDEL_TRANCHES} \
	-o ${RESULTS_DIR}/${INDEL_RECAL_VCF} \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log 

	
ARIND_RET=$?


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${INDEL_RECAL_VCF}

BGZ_RET=$?


log "TABIX gzip'd file"

tabix -f ${RESULTS_DIR}/${INDEL_RECAL_VCF_GZ}

TBX_RET=$?




# If the mapping gave an error, return, else cleanup
if [ $ARIND_RET -ne 0 ] 
then
	exit 15
elif [ $BGZ_RET -ne 0 ]
then
	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi

