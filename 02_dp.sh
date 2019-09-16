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

log "Getting coverage on autosomes and sex-chromosomes" 3
bcftools query -r $(for CHR in `seq 1 22` ; do echo -n "chr${CHR}," ; done)chrX,chrY -f '%DP\n' ${RESULTS_DIR}/${GENO_VCF_GZ} > ${RESULTS_DIR}/${GENO_DP_TXT}
BCF_RET=$?

if [ $BCF_RET -ne 0 ]
then
	exit 15
fi

log "Calculating max coverge threshold" 3

DP_THRESH=$(Rscript --vanilla -e "library(data.table) ; data <-fread('${RESULTS_DIR}/${GENO_DP_TXT}') ; data <- as.numeric(unlist(data)) ; cat(mean(data) + 5*sd(data))")
R_RET=$?


if [ $R_RET -ne 0 ]
then
	exit 15
fi


log "Hard filter on Coverage (> ${DP_THRESH})" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T VariantFiltration \
	-R ${REF_FASTA} \
	-V ${RESULTS_DIR}/${GENO_VCF_GZ} \
	--filterExpression "DP > ${DP_THRESH}" \
	--filterName "HighDP" \
	-o ${RESULTS_DIR}/${DP_VCF} \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log 


VFOTH_RET=$?

NUM_REMOVED=$(Rscript --vanilla -e "library(data.table) ; data <- fread('${RESULTS_DIR}/${GENO_DP_TXT}') ; sum(data > ${DP_THRESH}) ")

log "${NUM_REMOVED} variants removed." 4


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${DP_VCF}

BGZ_RET=$?


log "TABIX gzip'd file" 3

tabix -f ${RESULTS_DIR}/${DP_VCF_GZ}

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

