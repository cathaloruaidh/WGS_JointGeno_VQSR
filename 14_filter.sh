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




### Filter the cohort file based on specified thresholds
# Run the vcftools command

vcftools \
	--gzvcf ${RESULTS_DIR}/${RENAME_VCF_GZ} \
	--out ${RESULTS_DIR}/${FILTER_VCF} \
	--temp ${TEMP_DIR} \
	--maf ${MAF_FILTER} \
	--hwe ${HWE_FILTER} \
	--max-missing ${LMISS_FILTER} \
	--remove-filtered-all \
	--recode \
	--recode-INFO-all \
	2> >(tee ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log >&2)

VCFT_RET=$?

log "Filtered Variants" 3

mv ${RESULTS_DIR}/${FILTER_VCF}.recode.vcf ${RESULTS_DIR}/${FILTER_VCF}


bgzip -f ${RESULTS_DIR}/${FILTER_VCF}

BGZ_RET=$?

log "BGZIP output file done" 3


tabix -f ${RESULTS_DIR}/${FILTER_VCF_GZ}

TBX_RET=$?

log "TABIX gzip'd file done" 3


# If the mapping gave an error, return, else cleanup
if [ $VCFT_RET -ne 0 ] 
then
	exit 15
elif [ $BGZ_RET -ne 0 ]
then
	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi

