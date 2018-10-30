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


log "Apply VQSR on SNPs at tranche filter level ${SNP_TRANCHE_LEVEL}" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T ApplyRecalibration \
	-R ${REF_FASTA} \
	-input ${RESULTS_DIR}/${SNP_VCF_GZ} \
	-mode SNP \
	--ts_filter_level ${SNP_TRANCHE_LEVEL} \
	-recalFile ${RESULTS_DIR}/${SNP_RECAL_TXT} \
	-tranchesFile ${RESULTS_DIR}/${SNP_TRANCHES} \
	-o ${RESULTS_DIR}/${SNP_RECAL_VCF} \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log 

	
ARSNP_RET=$?


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${SNP_RECAL_VCF}

BGZ_RET=$?


log "TABIX gzip'd file" 3

tabix -f ${RESULTS_DIR}/${SNP_RECAL_VCF_GZ}

TBX_RET=$?


mv ${RESULTS_DIR}/${SNP_TRANCHES}.pdf ${GRAPHICS_DIR}



# If the mapping gave an error, return, else cleanup
if [ $ARSNP_RET -ne 0 ] 
then
	exit 15
elif [ $BGZ_RET -ne 0 ]
then
	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi

