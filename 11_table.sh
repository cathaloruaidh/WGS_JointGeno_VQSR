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


log "Convert variants to Table" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T VariantsToTable \
	-R ${REF_FASTA} \
	-V ${RESULTS_DIR}/${RECAL_VCF_GZ} \
	-F CHROM -F POS -F ID -F QUAL -F TYPE -F FILTER \
	-F DP -F FS -F InbreedingCoeff -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR \
	-raw \
	-o ${RESULTS_DIR}/${RECAL_TABLE} \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log
	
VTT_RET=$?



# If the mapping gave an error, return, else cleanup
if [ $VTT_RET -ne 0 ] 
then
	exit 15
fi

