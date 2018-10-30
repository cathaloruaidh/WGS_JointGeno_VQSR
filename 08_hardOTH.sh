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


log "Hard filter on other variants" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T VariantFiltration \
	-R ${REF_FASTA} \
	-V ${RESULTS_DIR}/${OTHER_VCF_GZ} \
	--filterExpression "QD < 2.0" \
	--filterName "OtherHardQD" \
	--filterExpression "FS > 200.0" \
	--filterName "OtherHardFS" \
	--filterExpression "SOR > 10.0" \
	--filterName "OtherHardSOR" \
	--filterExpression "ReadPosRankSum < -20.0" \
	--filterName "OtherHardReadPosRankSum" \
	-o ${RESULTS_DIR}/${OTHER_HARD_VCF} \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log 


VFOTH_RET=$?


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${OTHER_HARD_VCF}

BGZ_RET=$?


log "TABIX gzip'd file" 3

tabix -f ${RESULTS_DIR}/${OTHER_HARD_VCF_GZ}

TBX_RET=$?






# If the mapping gave an error, return, else cleanup
if [ $VFOTH_RET -ne 0 ] 
then
	exit 15
elif [ $BGZ_RET -ne 0 ]
then
	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi

