#!/bin/bash

############################################################
#
# Assembly of NGS data from FASTQ to GVCF 
# Cathal Ormond 2017. 
#
# Initial set up and creation of parameters. 
#
# Based heavily on scripts found here:
# https://github.com/genepi-freiburg/gwas
#
############################################################






### Text variables
PASS_TEST_LIGHT="[\e[102mPASSED\e[0m]"
PASS_TEST="[\e[42mPASSED\e[0m]"
FAIL_TEST_LIGHT="[\e[101mFAILED\e[0m]"
FAIL_TEST="[\e[41mFAILED\e[0m]"




### Processors
TPROCS=$(grep -c ^processor /proc/cpuinfo)
#TPROCS=14

if [ -z "${NPROCS}" ] 
then
	NPROCS=16
fi

APROCS=$(( ${NPROCS} - 1 ))


log "Using ${NPROCS} processors out of a total ${TPROCS}" 3
log " " 3




# Max Memory used by one process
# set to > 1/5 of Vishnu's max
JAVA_TEMP="-Djava.io.tmpdir=${TEMP_DIR}"
JAVA_OPTIONS=" -Xms${JAVA_MIN} -Xmx${JAVA_MAX} ${JAVA_TEMP}"





### Script directory
# the QC scripts live here. 
# Test if the scripts directory exists
if [ ! -d "${CURR_SCRIPT_DIR}" ]
then
	log "Problem with script directory" 1
	exit 3
fi




### Global Directories and Files
# Test if the source directory exists
if [ ! -d ${SOURCE_DIR} ]
then
	log "Problem with source directory: ${SOURCE_DIR}" 1
	exit 8
fi

log  "Currently directory: ${SOURCE_DIR}" 4
log " " 4 




### Create required files and directories
# Results directory for all results files
if [ ! -d ${RESULTS_DIR} ]
then
	log "Creating results directory: ${RESULTS_DIR##*/}" 3
	mkdir -p ${RESULTS_DIR}

else
	log "Results directory already exists at: ${RESULTS_DIR##*/}" 4
fi



# Temporary directory for any intermediate steps
if [ ! -d ${TEMP_DIR} ]
then
	log "Creating temp directory: ${TEMP_DIR##*/}" 3
	mkdir -p ${TEMP_DIR}

else
	log "Temp directory already exists at: ${TEMP_DIR##*/}" 4 
fi



# Graphics directory for all images for further analysis
if [ ! -d ${GRAPHICS_DIR} ]
then
	log "Creating graphics directory: ${GRAPHICS_DIR##*/}" 3
	mkdir -p ${GRAPHICS_DIR}

else
	log "Graphics directory already exists at: ${GRAPHICS_DIR##*/}" 4
fi



log " " 3



### Set the program paths

GATK34_FILE=${TOOL_DIR}/gatk/GenomeAnalysisTK_3.4.jar
GATK37_FILE=${TOOL_DIR}/gatk/GenomeAnalysisTK_3.7.jar
GATK38_FILE=${TOOL_DIR}/gatk/GenomeAnalysisTK_3.8.jar



# if the variables are empty, throw an error
if [ ! -x $(command -v plink) ]
then
	log "Error: plink not found." 1
	exit 6
fi



if [ ! -x $(command -v Rscript) ]
then
	log "Error: Rscript not found" 1
	exit 7
fi

R=$(command -v Rscript)



if [ ! -x $(command -v bwa) ]
then
	log "Error: bwa not found." 1
	exit 6
fi



if [ ! -x $(command -v samtools) ]
then
	log "Error: samtools not found." 1
	exit 6
fi



if [ ! -f ${GATK34_FILE} ] && [ ! -f ${GATK37_FILE} ] && [ ! -f ${GATK38_FILE} ]
then
	log "Error: gatk 3.4 and gatk 3.7 not found." 1
	exit 6
fi
#log "GATK was successfully found." 3



# GATK versions older than 3.8 cannot deal with spanning deletions, so default
# to this, and use older versions if necessary. 


if [ ! -f ${GATK_FILE} ]
then
	GATK_FILE=${GATK38_FILE}
fi


if [ ! -f ${GATK_FILE} ]
then
	GATK_FILE=${GATK37_FILE}
fi


if [ ! -f ${GATK_FILE} ]
then
	GATK_FILE=${GATK34_FILE}
fi






if [[ $( echo ${BUILD} | grep -E "19" ) ]]
then
	GRCH38=${REF_DIR}/${BUILD}/ucsc.hg19.fasta
	DBSNP146=${REF_DIR}/${BUILD}/dbsnp_138.hg19.vcf.gz
	DBSNP150=${REF_DIR}/${BUILD}/dbsnp_150.hg19.vcf.gz
	INDELS=${REF_DIR}/${BUILD}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
	THOUSANDG=${REF_DIR}/${BUILD}/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
	OMNI=${REF_DIR}/${BUILD}/1000G_omni2.5.hg19.sites.vcf.gz
	HAPMAP=${REF_DIR}/${BUILD}/hapmap_3.3.hg19.sites.vcf.gz
else
	GRCH38=${REF_DIR}/${BUILD}/GRCh38_full_analysis_set_plus_decoy_hla.fa
	DBSNP146=${REF_DIR}/${BUILD}/dbsnp_146.hg38.vcf.gz
	DBSNP150=${REF_DIR}/${BUILD}/dbsnp_150.hg38.chr.vcf.gz
	INDELS=${REF_DIR}/${BUILD}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
	#THOUSANDG=${REF_DIR}/${BUILD}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
	THOUSANDG=${REF_DIR}/${BUILD}/1000G_phase3_20170504.sites.GRCh38.vcf.gz
	OMNI=${REF_DIR}/${BUILD}/1000G_omni2.5.hg38.vcf.gz
	HAPMAP=${REF_DIR}/${BUILD}/hapmap_3.3.hg38.vcf.gz

fi
### Find the Resource Files


if [[ -z "${REF_FASTA}" ]]
then
	REF_FASTA=${GRCH38}
	log "Reference file $(basename ${REF_FASTA})" 4
elif [[ ! -f ${REF_FASTA} ]]
then
	REF_FASTA=${GRCH38}
	log "Specified reference file does not exist, using ${REF_FASTA}" 2
else
	log "Reference file $(basename ${REF_FASTA})" 4
fi



if [[ -z "${DBSNP}" ]]
then
	DBSNP=${DBSNP150}
	log "dbSNP file $(basename ${DBSNP})" 4
elif [[ ! -f ${DBSNP} ]]
then
	DBSNP=${DBSNP150}
	log "Specified dbSNP file does not exist, using ${DBSNP}" 2
else
	log "dbSNP file $(basename ${DBSNP})" 4
fi



if [[ ! -f ${INDELS} ]]
then
	log "The indel file does not exist at ${INDELS}" 1
	exit 1
fi

log "Indel reference file $(basename ${INDELS})" 4


if [[ ! -f ${THOUSANDG} ]]
then
	log "The 1000G reference VCF file does not exist at ${THOUSANDG}" 1
	exit 1
fi

log "1000G VCF reference file $(basename ${THOUSANDG})" 4


if [[ ! -f ${OMNI} ]]
then
	log "The Omni reference VCF file does not exist at ${OMNI}" 1
	exit 1
fi

log "Omni reference file $(basename ${OMNI})" 4


if [[ ! -f ${HAPMAP} ]]
then
	log "The HapMap reference VCF file does not exist at ${HAPMAP}" 1
	exit 1
fi

log "HapMap reference file $(basename ${HAPMAP})" 4

log " " 4






### Source name
# The base name of the raw data files.
SOURCE_FILE=$1



if [[ ! -f ${SOURCE_DIR}/${SOURCE_FILE} ]]
then
	log "Cannot find the input text file of gVCF files. Exiting ... " 1
	exit 1
fi


log "Found input files in: ${SOURCE_FILE}" 3


for FILE in `cat ${SOURCE_DIR}/${SOURCE_FILE}`
do
	if [[ ! -f ${FILE} ]]
	then
		log "Cannot find the gVCF file ${FILE}. Exiting ..." 1
		exit 1
	fi

	GVCF_BLOCKS=$( zcat ${FILE} | head -n 10000 | grep -v ^## | grep -v contig | grep NON_REF | wc -l )

	if [[ $GVCF_BLOCKS -eq 0 ]]
	then
		if [[ ${E} -le 1 || ${T} -eq 1 ]]
		then
			log "The file ${FILE} has no gVCF blocks, and joint genotyping is being called. Exiting ..." 1
			exit 1
		fi
#	else
#		if [[ ${S} -ge 1 || ${T} -gt 1 ]]
#		then
#			log "The file ${FILE} is a gVCF file, and has needs to be jointly genotyped. Exiting ..." 1
#			exit 1
#		fi
	fi
done


log "Outputting to file: ${OUTPUT_FILE}" 3


if [[ ${CLEAN} = true ]]
then 
	log "Cleaning input files once each stage is done!" 2
else
	log "Keeping all output files. " 3
fi

log " " 3

### Name the output files

GENO_VCF=${OUTPUT_FILE}.geno.vcf
GENO_VCF_GZ=${OUTPUT_FILE}.geno.vcf.gz

GENO_DP_TXT=${OUTPUT_FILE}.geno.dp.txt
DP_VCF=${OUTPUT_FILE}.geno.dp.vcf
DP_VCF_GZ=${OUTPUT_FILE}.geno.dp.vcf.gz

SNP_VCF=${OUTPUT_FILE}.geno.snp.vcf
SNP_VCF_GZ=${OUTPUT_FILE}.geno.snp.vcf.gz

INDEL_VCF=${OUTPUT_FILE}.geno.indel.vcf
INDEL_VCF_GZ=${OUTPUT_FILE}.geno.indel.vcf.gz

OTHER_VCF=${OUTPUT_FILE}.geno.other.vcf
OTHER_VCF_GZ=${OUTPUT_FILE}.geno.other.vcf.gz

SNP_RECAL_TXT=${OUTPUT_FILE}.geno.snp.recal.txt
SNP_TRANCHES=${OUTPUT_FILE}.geno.snp.tranches.txt
SNP_PLOTS=${OUTPUT_FILE}.geno.snp.plots.R
SNP_RECAL_VCF=${OUTPUT_FILE}.geno.snp.recal.vcf
SNP_RECAL_VCF_GZ=${OUTPUT_FILE}.geno.snp.recal.vcf.gz

INDEL_RECAL_TXT=${OUTPUT_FILE}.geno.indel.recal.txt
INDEL_TRANCHES=${OUTPUT_FILE}.geno.indel.tranches.txt
INDEL_PLOTS=${OUTPUT_FILE}.geno.indel.plots.R
INDEL_RECAL_VCF=${OUTPUT_FILE}.geno.indel.recal.vcf
INDEL_RECAL_VCF_GZ=${OUTPUT_FILE}.geno.indel.recal.vcf.gz

OTHER_HARD_VCF=${OUTPUT_FILE}.geno.other.hard.vcf
OTHER_HARD_VCF_GZ=${OUTPUT_FILE}.geno.other.hard.vcf.gz

RECAL_VCF=${OUTPUT_FILE}.recal.vcf
RECAL_VCF_GZ=${OUTPUT_FILE}.recal.vcf.gz
RECAL_TABLE=${OUTPUT_FILE}.recal.table

RENAME_VCF=${OUTPUT_FILE}.recal.rename.vcf
RENAME_VCF_GZ=${OUTPUT_FILE}.recal.rename.vcf.gz
RENAME_TABLE=${OUTPUT_FILE}.recal.rename.table



REPORT_TEX=${OUTPUT_FILE}.report.tex

FILTER_VCF=${OUTPUT_FILE}.recal.rename.filtered.vcf
FILTER_VCF_GZ=${OUTPUT_FILE}.recal.rename.filtered.vcf.gz


PLINK_FILE=${OUTPUT_FILE}.recal



