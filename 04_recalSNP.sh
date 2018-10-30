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


log "Calculate VQSR on SNPs" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T VariantRecalibrator \
	-R ${REF_FASTA} \
	-input ${RESULTS_DIR}/${SNP_VCF_GZ} \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 ${OMNI} \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${THOUSANDG} \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
	-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.99 -tranche 99.98 -tranche 99.97 -tranche 99.96 \
	-tranche 99.95 -tranche 99.94 -tranche 99.93 -tranche 99.92 -tranche 99.91 \
	-tranche 99.90 -tranche 99.80 -tranche 99.70 -tranche 99.60 -tranche 99.50 \
	-tranche 99.00 -tranche 98.00 -tranche 97.00 -tranche 96.00 -tranche 95.00 \
	-tranche 94.00 -tranche 93.00 -tranche 92.00 -tranche 91.00 -tranche 90.00 \
	-recalFile ${RESULTS_DIR}/${SNP_RECAL_TXT} \
	-tranchesFile ${RESULTS_DIR}/${SNP_TRANCHES} \
	-rscriptFile ${GRAPHICS_DIR}/${SNP_PLOTS} \
	--maxGaussians 4 \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log 
	
VRSNP_RET=$?




# If the mapping gave an error, return, else cleanup
if [ $VRSNP_RET -ne 0 ] 
then
	exit 15
fi

