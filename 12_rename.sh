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


log "Rename samples in VCF to FAM names" 3



if [[ ! -f ${OUTPUT_FILE}.samples.txt ]]
then
	for LINE in `cat ${SOURCE_DIR}/${SOURCE_FILE}`
	do 
		#BASE=${LINE##*/}
		#echo ${BASE%.g.vcf.gz}
		echo -e "${LINE##*/} $(bcftools query -l ${LINE})" 
	done  \
		| LC_ALL=C sort -k2,2 \
		| awk '{print $1}' \
		| sed -e 's/.g.vcf.gz//' \
		> ${OUTPUT_FILE}.samples.txt
fi

bcftools reheader -s ${OUTPUT_FILE}.samples.txt -o ${RESULTS_DIR}/${RENAME_VCF_GZ} ${RESULTS_DIR}/${RECAL_VCF_GZ}
	
REHEAD_RET=$?

exit 0

#log "BGZIP output file" 3

#bgzip -f ${RESULTS_DIR}/${RENAME_VCF}

#BGZ_RET=$?


log "TABIX gzip'd file" 3

tabix -f ${RESULTS_DIR}/${RENAME_VCF_GZ}

TBX_RET=$?




# If the mapping gave an error, return, else cleanup
if [ $REHEAD_RET -ne 0 ] 
then
	exit 15
#elif [ $BGZ_RET -ne 0 ]
#then
#	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi

