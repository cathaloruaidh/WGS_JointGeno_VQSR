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


# Note that we will not be using the mitochondrial chromosomes, so we will
# extract it now and save it in case it is used again.

vcftools --gzvcf ${RESULTS_DIR}/${GENO_VCF_GZ} \
	--chr chrM \
	--recode \
	--recode-INFO-all \
	--out ${RESULTS_DIR}/${GENO_VCF_GZ%vcf.gz}chrM

mv ${RESULTS_DIR}/${GENO_VCF_GZ%vcf.gz}chrM.recode.vcf ${RESULTS_DIR}/${GENO_VCF_GZ%vcf.gz}chrM.vcf
bgzip -f ${RESULTS_DIR}/${GENO_VCF_GZ%vcf.gz}chrM.vcf

log "Extracted chrM for future reference." 4


### Extract all SNPs from variant file
# Run the GATK command


log "Selecting SNPs" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T SelectVariants \
	-R ${REF_FASTA} \
	--variant ${RESULTS_DIR}/${DP_VCF_GZ} \
	-o ${RESULTS_DIR}/${SNP_VCF} \
	-selectType SNP \
	$(for CHR in `seq 1 22` X Y ; do echo -n "-L chr${CHR} " ; done) \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.SNP.log 
	
GENO_RET=$?


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${SNP_VCF}

BGZ_RET=$?


log "TABIX the gzip'd file" 3

tabix -f ${RESULTS_DIR}/${SNP_VCF_GZ}

TBX_RET=$?



# If the extraction gave an error, exit, else continue
if [ $GENO_RET -ne 0 ] 
then
	exit 15
elif [ $BGZ_RET -ne 0 ]
then
	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi





### Extract all SNPs from variant file
# Run the GATK command


log " " 3
log "Selecting indels" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T SelectVariants \
	-R ${REF_FASTA} \
	--variant ${RESULTS_DIR}/${DP_VCF_GZ} \
	-o ${RESULTS_DIR}/${INDEL_VCF} \
	-selectType INDEL \
	$(for CHR in `seq 1 22` X Y ; do echo -n "-L chr${CHR} " ; done) \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.INDEL.log 
	
GENO_RET=$?


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${INDEL_VCF}

BGZ_RET=$?


log "TABIX gzip'd file" 3

tabix -f ${RESULTS_DIR}/${INDEL_VCF_GZ}

TBX_RET=$?



# If the extraction gave an error, exit, else continue
if [ $GENO_RET -ne 0 ] 
then
	exit 15
elif [ $BGZ_RET -ne 0 ]
then
	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi





### Extract all SNPs from variant file
# Run the GATK command


log " " 3
log "Selecting other variants" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T SelectVariants \
	-R ${REF_FASTA} \
	--variant ${RESULTS_DIR}/${DP_VCF_GZ} \
	-o ${RESULTS_DIR}/${OTHER_VCF} \
	-xlSelectType SNP \
	-xlSelectType INDEL\
	$(for CHR in `seq 1 22` X Y ; do echo -n "-L chr${CHR} " ; done) \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.OTHER.log \
#	-nt ${NPROCS}
	
GENO_RET=$?


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${OTHER_VCF}

BGZ_RET=$?


log "TABIX gzip'd file" 3

tabix -f ${RESULTS_DIR}/${OTHER_VCF_GZ}

TBX_RET=$?



# If the extraction gave an error, exit, else continue
if [ $GENO_RET -ne 0 ] 
then
	exit 15
elif [ $BGZ_RET -ne 0 ]
then
	exit 15
elif [ $TBX_RET -ne 0 ]
then
	exit 15
fi

