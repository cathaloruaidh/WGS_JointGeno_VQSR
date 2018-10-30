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




### Joint Genotyping
# Run the GATK command

VARIANTS=""

for f in `cat ${SOURCE_DIR}/${SOURCE_FILE}`
do
	VARIANTS="${VARIANTS} --variant ${f}"
	log "Using file ${f}" 3
done


log "Joing Genotyping" 3

java ${JAVA_OPTIONS} -jar ${GATK_FILE} -T GenotypeGVCFs \
	-R ${REF_FASTA} \
	${VARIANTS} \
	--annotation InbreedingCoeff \
	--annotation FisherStrand \
	--annotation QualByDepth \
	--annotation ChromosomeCounts \
	--annotation StrandOddsRatio \
	--dbsnp ${DBSNP} \
	-o ${RESULTS_DIR}/${GENO_VCF} \
	--standard_min_confidence_threshold_for_calling 10.0 \
	--downsample_to_coverage 1000 \
	--downsampling_type BY_SAMPLE \
	-log ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.log 
	
GENO_RET=$?


log "BGZIP output file" 3

bgzip -f ${RESULTS_DIR}/${GENO_VCF}

BGZ_RET=$?


log "TABIX gzip'd file" 3

tabix -f ${RESULTS_DIR}/${GENO_VCF_GZ}

TBX_RET=$?



# If the mapping gave an error, return, else cleanup
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

