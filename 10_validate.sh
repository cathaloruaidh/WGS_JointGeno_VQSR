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


log "Validate the merged file." 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T ValidateVariants \
	-R ${REF_FASTA} \
	--dbsnp ${DBSNP} \
	--reference_window_stop 300 \
	-V ${RESULTS_DIR}/${RECAL_VCF_GZ} \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log

VV_RET=$?




# If the mapping gave an error, return, else cleanup
if [ ${VV_RET} -ne 0 ]
then
	exit 15
fi

