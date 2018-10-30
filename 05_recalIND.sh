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


log "Calculate VQSR on indels" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T VariantRecalibrator \
	-R ${REF_FASTA} \
	-input ${RESULTS_DIR}/${INDEL_VCF_GZ} \
	-resource:mills,known=false,training=true,truth=true,prior=12.0 ${INDELS} \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
	-an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
	-mode INDEL \
	-tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 \
	-tranche 99.95 -tranche 99.94 -tranche 99.93 -tranche 99.92 -tranche 99.91 \
	-tranche 99.90 -tranche 99.80 -tranche 99.70 -tranche 99.60 -tranche 99.50 \
	-tranche 99.00 -tranche 98.00 -tranche 97.00 -tranche 96.00 -tranche 95.00 \
	-tranche 94.00 -tranche 93.00 -tranche 92.00 -tranche 91.00 -tranche 90.00 \
	--maxGaussians 4 \
	-recalFile ${RESULTS_DIR}/${INDEL_RECAL_TXT} \
	-tranchesFile ${RESULTS_DIR}/${INDEL_TRANCHES} \
	-rscriptFile ${GRAPHICS_DIR}/${INDEL_PLOTS} \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log 
	
VRIND_RET=$?




# If the mapping gave an error, return, else cleanup
if [ $VRIND_RET -ne 0 ] 
then
	exit 15
fi

