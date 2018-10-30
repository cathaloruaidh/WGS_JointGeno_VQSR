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




### Generate information files
# Run the vcftools command

vcftools \
	--gzvcf ${RESULTS_DIR}/${RENAME_VCF_GZ} \
	--out ${RESULTS_DIR}/${PLINK_FILE} \
	--temp ${TEMP_DIR} \
	--plink
	2> >(tee ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.freq.log >&2)

VCFT_RET=$?

log "Converted to Plink format. " 4


# If the mapping gave an error, return, else cleanup
if [ $VCFT_RET -ne 0 ] 
then
	exit 15
fi



# create the correct PED file manually from the output PED and the FAM files
if [[ ! -f ${FAM_FILE} ]]
then
	log "No FAM file was passed, converting PLINK files to binary files as they are." 3
	plink --file ${RESULTS_DIR}/${PLINK_FILE} --out ${RESULTS_DIR}/${PLINK_FILE} --make-bed 
else
	log "Editing PED file to include information from FAM file." 3
	log "Note: if the PED file has different number of columns for each line, the max will be used." 4


	rm -f ${RESULTS_DIR}/${PLINK_FILE}.family.ped
	touch ${RESULTS_DIR}/${PLINK_FILE}.family.ped

	log "Created family file: ${RESULTS_DIR}/${PLINK_FILE}.family.ped" 4

	# number of genotype columns is the max number of columns in the PED file minus 6
	COLS=$(( $(awk '{print NF}' ${RESULTS_DIR}/${PLINK_FILE}.ped | sort -nru | head -n 1) - 6 ))
	EMPTY=""

	log "Counted" 4

	# If we have no information, output "0" COLS times 
	EMPTY=$(yes "0" | head -n ${COLS} | tr '\n' '\t')
	log "Created EMPTY line (${COLS} values, length=${#EMPTY})" 4

	while IFS='' read -r LINE || [[ -n ${LINE}  ]]
	do
		# Get FID and IID from PED file
		FID=$(echo ${LINE} | awk '{print $1}')
		IID=$(echo ${LINE} | awk '{print $2}')

		# get the genotype information from PED file
		GENO_INFO=$(awk -v IID=${IID} '{if($2 == IID) print $0}' ${RESULTS_DIR}/${PLINK_FILE}.ped | cut -f 7- --output-delimiter=$'\t')

		# if the sample isn't in the PED file, use empty genotype information
		if [[ -z "${GENO_INFO}" ]]
		then
			#log "Individual ${IID} from family ${FID} not found! Empty genotype." 4
			GENO_INFO=$( echo -en ${EMPTY} )
		fi
		
		# add entry to family PED file
		echo -e "${LINE}\t${GENO_INFO}" >> ${RESULTS_DIR}/${PLINK_FILE}.family.ped

		log "Added individual ${IID}" 4
	done < ${FAM_FILE}
fi

cp ${RESULTS_DIR}/${PLINK_FILE}.map ${RESULTS_DIR}/${PLINK_FILE}.family.map


# Make binary plink files
log "Make binary files from the PED and MAP files. " 3

plink --file ${RESULTS_DIR}/${PLINK_FILE} --make-bed --out ${RESULTS_DIR}/${PLINK_FILE}
plink --file ${RESULTS_DIR}/${PLINK_FILE}.family --make-bed --out ${RESULTS_DIR}/${PLINK_FILE}.family 


#plink --bfile ${RESULTS_DIR}/${PLINK_FILE} --missing --het --freq --hardy --out ${RESULTS_DIR}/${PLINK_FILE}
#plink --bfile ${RESULTS_DIR}/${PLINK_FILE}.family --maf 0.05 --check-sex --genome rel-check full --out ${RESULTS_DIR}/${PLINK_FILE}


tabix -f ${RESULTS_DIR}/${RENAME_VCF_GZ}

peddy --plot -p ${NPROCS} --sites hg38 --prefix ${RESULTS_DIR}/${PLINK_FILE} ${RESULTS_DIR}/${RENAME_VCF_GZ} ${FAM_FILE} 2> >(tee ${LOG_DIR}/${OUTPUT_FILE}.JOINT.${1}.peddy.log >&2)
mv ${RESULTS_DIR}/${PLINK_NAME}*.png ${GRAPHICS_DIR}



# If the mapping gave an error, return, else cleanup
if [ $VCFT_RET -ne 0 ] 
then
	exit 15
fi

