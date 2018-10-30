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


log "Merge SNP, indel and other variants files" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T CombineVariants \
	-R ${REF_FASTA} \
	-V:snp ${RESULTS_DIR}/${SNP_RECAL_VCF_GZ} \
	-V:indel ${RESULTS_DIR}/${INDEL_RECAL_VCF_GZ} \
	-V:other ${RESULTS_DIR}/${OTHER_HARD_VCF_GZ} \
	-o ${RESULTS_DIR}/${RECAL_VCF} \
	-assumeIdenticalSamples \
	-genotypeMergeOptions PRIORITIZE \
	-priority snp,indel,other \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log
	
CV_RET=$?


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${RECAL_VCF}

BGZ_RET=$?


log "TABIX gzip'd file" 3

tabix -f ${RESULTS_DIR}/${RECAL_VCF_GZ}

TBX_RET=$?


mv ${RESULTS_DIR}*pdf ${GRAPHICS_DIR}



# If the mapping gave an error, return, else cleanup
if [ $CV_RET -ne 0 ] 
then
	exit 15
elif [ $BGZ_RET -ne 0 ]
then
	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi

